#!/usr/bin/env python3
import argparse
import json
import math
import os
from pathlib import Path
from urllib.request import urlretrieve

import geopandas as gpd
import pygeoops
from shapely.geometry import LineString, MultiLineString
from shapely.ops import unary_union
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
		'stations': 0,
		'measurements': 0,
		'output_segments': 0,
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
					'source': 'pygeoops.centerline'
				}
			})

			measurements = measure_line(line, sis_union, fmzk_id, part_index, config)
			summary['stations'] += measurements['stations']
			summary['measurements'] += len(measurements['points'])
			segments = group_measurements(measurements['points'], config, fmzk_id, part_index)
			summary['output_segments'] += len(segments)
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
	length = line.length
	if length <= 0:
		return {'stations': 0, 'points': []}

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
				width = seg.length
				points.append({
					'point': (mid.x, mid.y),
					'chain_m': s,
					'width_m': width,
					'width_class': width_class(width, config['width_classes']),
					'fmzk_id': fmzk_id,
					'part_index': part_index
				})
		stations += 1
		s += step
	return {'stations': stations, 'points': points}


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


def group_measurements(points, config, fmzk_id, part_index):
	if not points:
		return []
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

		if new_group:
			if current:
				feature = make_segment_feature(current, config, fmzk_id, part_index)
				if feature:
					segments.append(feature)
			current = [p]
		else:
			current.append(p)
		last = p
	if current:
		feature = make_segment_feature(current, config, fmzk_id, part_index)
		if feature:
			segments.append(feature)
	return segments


def make_segment_feature(items, config, fmzk_id, part_index):
	if len(items) < 2:
		return None
	coords = [item['point'] for item in items]
	line = LineString(coords)
	tol = config['measurement'].get('simplify_tolerance_m', 0)
	if tol and tol > 0:
		line = line.simplify(tol, preserve_topology=False)
	if line.length <= 0:
		return None
	widths = [item['width_m'] for item in items]
	return {
		'type': 'Feature',
		'geometry': mapping_line(line),
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
			'method': 'centerline_widths_v1'
		}
	}


def mapping_line(line):
	return {'type': 'LineString', 'coordinates': [(float(x), float(y)) for x, y in line.coords]}


def write_geojsonseq(path, features):
	Path(path).parent.mkdir(parents=True, exist_ok=True)
	with open(path, 'w', encoding='utf-8') as f:
		for feat in features:
			f.write(json.dumps(feat, ensure_ascii=False))
			f.write('\n')


if __name__ == '__main__':
	main()
