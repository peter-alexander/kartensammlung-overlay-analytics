#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path
from urllib.request import urlretrieve

import geopandas as gpd
import pygeoops
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import linemerge, unary_union
from shapely import make_valid


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', default='config/centerline-widths.json')
	args = parser.parse_args()

	config = load_json(args.config)
	paths = config['paths']
	work_dir = Path(paths['work_dir'])
	output_dir = Path(paths['output_dir'])
	work_dir.mkdir(parents=True, exist_ok=True)
	output_dir.mkdir(parents=True, exist_ok=True)

	fmzk_raw = work_dir / 'fmzk_raw.geojson'
	sis_raw = work_dir / 'sis_raw.geojson'
	download(config['sources']['fmzk']['url'], fmzk_raw)
	download(config['sources']['sis']['url'], sis_raw)

	fmzk = load_prepare_source(fmzk_raw, config, config['sources']['fmzk'])
	sis = load_prepare_source(sis_raw, config, config['sources']['sis'])

	centerline_features = []
	measurement_features = []
	summary = {
		'method': config['method'],
		'fmzk_features': int(len(fmzk)),
		'sis_features': int(len(sis)),
		'centerline_ok': 0,
		'centerline_failed': 0,
		'centerline_parts': 0,
		'centerline_merged': 0,
		'centerline_post_simplified': 0,
		'centerline_straightened': 0,
		'centerline_straight_simplified': 0,
		'stations': 0,
		'measurements': 0,
		'output_segments': 0,
		'output_points_centerline': 0,
		'output_points_sis_midpoint': 0,
		'split_by_invalid_connection': 0,
		'dropped_single_point_groups': 0,
		'dropped_invalid_segments': 0,
		'fallback_unsimplified_segments': 0,
		'errors': {}
	}

	sis_union = unary_union(list(sis.geometry)) if len(sis) else None
	if sis_union is None or sis_union.is_empty:
		raise RuntimeError('SIS union is empty')

	for idx, row in fmzk.iterrows():
		geom = row.geometry
		props = row.drop(labels=['geometry']).to_dict()
		fmzk_id = props.get('FMZK_ID') or props.get('OBJECTID') or props.get('ID') or idx

		try:
			geom_clean = prepare_centerline_geom(geom, config)
			if geom_clean is None or geom_clean.is_empty:
				continue
			cl = pygeoops.centerline(
				geom_clean,
				densify_distance=config['centerline']['densify_distance'],
				min_branch_length=config['centerline']['min_branch_length'],
				simplifytolerance=config['centerline']['simplifytolerance'],
				extend=config['centerline']['extend']
			)
			lines = explode_lines(cl)
			lines = postprocess_centerlines(lines, config, summary)
			lines = [l for l in lines if l.length >= config['centerline']['min_line_length_m']]
			summary['centerline_ok'] += 1
		except Exception as exc:
			summary['centerline_failed'] += 1
			key = type(exc).__name__
			summary['errors'][key] = summary['errors'].get(key, 0) + 1
			continue

		summary['centerline_parts'] += len(lines)
		for part_index, line in enumerate(lines):
			centerline_features.append({
				'type': 'Feature',
				'geometry': mapping_line(line),
				'properties': {
					'fmzk_id': fmzk_id,
					'part_index': part_index,
					'length_m': round(line.length, 2),
					'source': 'pygeoops.centerline',
					'post_simplify_tolerance_m': config['centerline'].get('post_simplify_tolerance_m', 0)
				}
			})

			measurements = measure_line(line, sis_union, fmzk_id, part_index, config)
			summary['stations'] += measurements['stations']
			summary['measurements'] += len(measurements['points'])
			summary['output_points_centerline'] += measurements['output_points_centerline']
			summary['output_points_sis_midpoint'] += measurements['output_points_sis_midpoint']
			segments, group_stats = group_measurements(measurements['points'], config, fmzk_id, part_index, sis_union)
			summary['output_segments'] += len(segments)
			merge_summary(summary, group_stats)
			measurement_features.extend(segments)

	write_geojsonseq(paths['centerlines_metric'], centerline_features)
	write_geojsonseq(paths['measurements_metric'], measurement_features)
	Path(paths['summary']).parent.mkdir(parents=True, exist_ok=True)
	Path(paths['summary']).write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
	print(json.dumps(summary, ensure_ascii=False))


def load_json(path):
	return json.loads(Path(path).read_text(encoding='utf-8'))


def download(url, target):
	if str(url).startswith('http://') or str(url).startswith('https://'):
		urlretrieve(url, target)
	else:
		target.write_bytes(Path(url).read_bytes())


def load_prepare_source(path, config, source_cfg):
	gdf = gpd.read_file(path)
	if gdf.crs is None:
		gdf = gdf.set_crs(config['crs']['input'])
	gdf = gdf.to_crs(config['crs']['metric'])

	prop = source_cfg.get('filter_property')
	values = source_cfg.get('filter_values') or []
	if prop and values:
		gdf = gdf[gdf[prop].isin(values)].copy()

	gdf['geometry'] = gdf.geometry.apply(lambda g: make_valid(g) if g is not None and not g.is_empty else g)
	gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()
	gdf = gdf.explode(index_parts=False).reset_index(drop=True)
	return gdf


def prepare_centerline_geom(geom, config):
	if geom is None or geom.is_empty:
		return None
	g = make_valid(geom)
	if g.is_empty:
		return None
	min_area = config['centerline'].get('min_area_m2', 0)
	if hasattr(g, 'area') and g.area < min_area:
		return None
	tol = config['centerline'].get('simplify_tolerance_m', 0)
	if tol and tol > 0:
		g = g.simplify(tol, preserve_topology=True)
		g = make_valid(g)
	return g


def postprocess_centerlines(lines, config, summary):
	lines = merge_connected_centerlines(lines, config, summary)
	out = []
	for line in lines:
		if line is None or line.is_empty or line.length <= 0:
			continue

		line2 = straighten_if_nearly_linear(line, config, summary)
		line3 = straight_simplify_if_nearly_linear(line2, config, summary)
		line4 = simplify_centerline(line3, config, summary)
		out.extend(explode_lines(line4))
	return [line for line in out if line is not None and not line.is_empty and line.length > 0]


def merge_connected_centerlines(lines, config, summary):
	if not config['centerline'].get('merge_connected_lines', False):
		return lines
	if len(lines) <= 1:
		return lines
	try:
		merged = linemerge(unary_union(lines))
		merged_lines = explode_lines(merged)
		if merged_lines and len(merged_lines) < len(lines):
			summary['centerline_merged'] += len(lines) - len(merged_lines)
			return merged_lines
		return merged_lines or lines
	except Exception:
		return lines


def straighten_if_nearly_linear(line, config, summary):
	cfg = config['centerline']
	if not cfg.get('straighten_nearly_linear', False):
		return line
	if len(line.coords) < 3:
		return line
	if line.length < cfg.get('straighten_min_length_m', 8.0):
		return line

	coords = list(line.coords)
	start = coords[0]
	end = coords[-1]
	base = LineString([start, end])
	if base.length <= 1e-9:
		return line

	length_ratio = line.length / base.length
	if length_ratio > cfg.get('straighten_max_length_ratio', 2.5):
		return line

	max_dev = max(base.distance(Point(coord)) for coord in coords[1:-1])
	if max_dev <= cfg.get('straighten_max_deviation_m', 3.0):
		summary['centerline_straightened'] += 1
		return base
	return line


def straight_simplify_if_nearly_linear(line, config, summary):
	cfg = config['centerline']
	if not cfg.get('straight_simplify_nearly_linear', False):
		return line
	if len(line.coords) < 4:
		return line

	tol = cfg.get('straight_simplify_tolerance_m', 0)
	if tol <= 0:
		return line

	simplified = line.simplify(tol, preserve_topology=False)
	if simplified is None or simplified.is_empty or not isinstance(simplified, LineString):
		return line
	if len(simplified.coords) >= len(line.coords) - cfg.get('straight_simplify_min_reduction', 2):
		return line

	max_hausdorff = cfg.get('straight_simplify_max_hausdorff_m', tol)
	if line.hausdorff_distance(simplified) > max_hausdorff:
		return line

	summary['centerline_straight_simplified'] += 1
	return simplified


def simplify_centerline(line, config, summary):
	tol = config['centerline'].get('post_simplify_tolerance_m', 0)
	preserve = config['centerline'].get('post_simplify_preserve_topology', False)
	if not tol or tol <= 0:
		return line

	simplified = line.simplify(tol, preserve_topology=preserve)
	if simplified is None or simplified.is_empty:
		return line
	if len(explode_lines(simplified)) == 0:
		return line
	if isinstance(simplified, LineString) and len(simplified.coords) < len(line.coords):
		summary['centerline_post_simplified'] += 1
	return simplified


def explode_lines(geom):
	if geom is None or geom.is_empty:
		return []
	if isinstance(geom, LineString):
		return [geom]
	if isinstance(geom, MultiLineString):
		return list(geom.geoms)
	if hasattr(geom, 'geoms'):
		out = []
		for item in geom.geoms:
			out.extend(explode_lines(item))
		return out
	return []


def measure_line(line, sis_union, fmzk_id, part_index, config):
	step = config['measurement']['step_m']
	half = config['measurement']['cross_section_half_length_m']
	points = []
	stations = 0
	output_points_centerline = 0
	output_points_sis_midpoint = 0
	length = line.length
	if length <= 0:
		return {'stations': 0, 'points': [], 'output_points_centerline': 0, 'output_points_sis_midpoint': 0}

	s = 0.0
	while s <= length:
		p = line.interpolate(s)
		tangent = tangent_at(line, s, step)
		if tangent is not None:
			nx, ny = -tangent[1], tangent[0]
			cross = LineString([(p.x - nx * half, p.y - ny * half), (p.x + nx * half, p.y + ny * half)])
			inter = cross.intersection(sis_union)
			segments = explode_lines(inter)
			segments = [seg for seg in segments if seg.length >= config['measurement']['min_width_m']]
			segments.sort(key=lambda seg: seg.centroid.distance(p))
			if segments:
				seg = segments[0]
				mid = seg.interpolate(0.5, normalized=True)
				width = round(seg.length, 2)
				out_point, source = choose_output_point(p, mid, seg, config)
				if source == 'centerline':
					output_points_centerline += 1
				else:
					output_points_sis_midpoint += 1
				points.append({
					'point': (out_point.x, out_point.y),
					'centerline_point': (p.x, p.y),
					'sis_midpoint': (mid.x, mid.y),
					'output_point_source': source,
					'chain_m': s,
					'width_m': width,
					'width_class': width_class(width, config['width_classes']),
					'fmzk_id': fmzk_id,
					'part_index': part_index
				})
		stations += 1
		s += step
	return {
		'stations': stations,
		'points': points,
		'output_points_centerline': output_points_centerline,
		'output_points_sis_midpoint': output_points_sis_midpoint
	}


def choose_output_point(centerline_point, sis_midpoint, measured_segment, config):
	mode = config['measurement'].get('output_geometry', 'sis_midpoint')
	if mode == 'centerline':
		return centerline_point, 'centerline'
	if mode == 'centerline_when_inside_sis' and measured_segment.buffer(1e-6).covers(centerline_point):
		return centerline_point, 'centerline'
	return sis_midpoint, 'sis_midpoint'


def tangent_at(line, s, delta):
	length = line.length
	a = max(0.0, s - delta)
	b = min(length, s + delta)
	if b - a <= 1e-9:
		return None
	pa = line.interpolate(a)
	pb = line.interpolate(b)
	dx = pb.x - pa.x
	dy = pb.y - pa.y
	n = math.hypot(dx, dy)
	if n <= 1e-12:
		return None
	return (dx / n, dy / n)


def width_class(width, classes):
	for cls in classes:
		min_v = cls.get('min')
		max_v = cls.get('max')
		if (min_v is None or width >= min_v) and (max_v is None or width < max_v):
			return cls['id']
	return None


def group_measurements(points, config, fmzk_id, part_index, sis_union):
	stats = {
		'split_by_invalid_connection': 0,
		'dropped_single_point_groups': 0,
		'dropped_invalid_segments': 0,
		'fallback_unsimplified_segments': 0
	}
	if not points:
		return [], stats
	points.sort(key=lambda p: p['chain_m'])
	segments = []
	current = []
	last = None
	for p in points:
		new_group = False
		if last is None:
			new_group = True
		elif p['width_class'] != last['width_class']:
			new_group = True
		elif p['chain_m'] - last['chain_m'] > config['measurement']['max_station_gap_m']:
			new_group = True
		elif not connection_inside_surface(last['point'], p['point'], sis_union, config):
			new_group = True
			stats['split_by_invalid_connection'] += 1

		if new_group:
			if current:
				feature = make_segment_feature(current, config, fmzk_id, part_index, sis_union, stats)
				if feature:
					segments.append(feature)
			current = [p]
		else:
			current.append(p)
		last = p
	if current:
		feature = make_segment_feature(current, config, fmzk_id, part_index, sis_union, stats)
		if feature:
			segments.append(feature)
	return segments, stats


def connection_inside_surface(a, b, sis_union, config):
	line = LineString([a, b])
	if line.length <= 1e-9:
		return True
	sample_step = config['measurement'].get('output_validation_sample_step_m', 0.5)
	steps = max(1, math.ceil(line.length / sample_step))
	for i in range(0, steps + 1):
		p = line.interpolate(i / steps, normalized=True)
		if not sis_union.covers(p):
			return False
	return True


def make_segment_feature(items, config, fmzk_id, part_index, sis_union, stats):
	if len(items) < 2:
		stats['dropped_single_point_groups'] += 1
		return None
	coords = [item['point'] for item in items]
	line = LineString(coords)
	if not line_inside_surface(line, sis_union, config):
		stats['dropped_invalid_segments'] += 1
		return None

	output_line = line
	tol = config['measurement'].get('simplify_tolerance_m', 0)
	if tol and tol > 0:
		simplified = line.simplify(tol, preserve_topology=False)
		if simplified.length > 0 and line_inside_surface(simplified, sis_union, config):
			output_line = simplified
		else:
			stats['fallback_unsimplified_segments'] += 1

	widths = [item['width_m'] for item in items]
	point_sources = {item.get('output_point_source', 'unknown') for item in items}
	return {
		'type': 'Feature',
		'geometry': mapping_line(output_line),
		'properties': {
			'fmzk_id': fmzk_id,
			'part_index': part_index,
			'width_class': items[0]['width_class'],
			'min_width_m': round(min(widths), 2),
			'max_width_m': round(max(widths), 2),
			'avg_width_m': round(sum(widths) / len(widths), 2),
			'count': len(items),
			'from_m': round(items[0]['chain_m'], 2),
			'to_m': round(items[-1]['chain_m'], 2),
			'geometry_source': '+'.join(sorted(point_sources)),
			'method': 'centerline_widths_v1'
		}
	}


def line_inside_surface(line, sis_union, config):
	if line.length <= 1e-9:
		return False
	sample_step = config['measurement'].get('output_validation_sample_step_m', 0.5)
	steps = max(1, math.ceil(line.length / sample_step))
	for i in range(0, steps + 1):
		p = line.interpolate(i / steps, normalized=True)
		if not sis_union.covers(p):
			return False
	return True


def mapping_line(line):
	return {'type': 'LineString', 'coordinates': [(float(x), float(y)) for x, y in line.coords]}


def write_geojsonseq(path, features):
	Path(path).parent.mkdir(parents=True, exist_ok=True)
	with open(path, 'w', encoding='utf-8') as f:
		for feat in features:
			f.write(json.dumps(feat, ensure_ascii=False))
			f.write('\n')


def merge_summary(summary, stats):
	for key, value in stats.items():
		if isinstance(value, int):
			summary[key] = summary.get(key, 0) + value


if __name__ == '__main__':
	main()
