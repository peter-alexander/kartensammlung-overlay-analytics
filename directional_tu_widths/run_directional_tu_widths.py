#!/usr/bin/env python3
import argparse
import json
import math
import os
import subprocess
from pathlib import Path
from urllib.request import urlretrieve

import geopandas as gpd
from shapely import make_valid
from shapely.geometry import LineString, MultiLineString, mapping
from shapely.ops import nearest_points


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', default='config/directional-tu-widths.json')
	args = parser.parse_args()

	config = normalize_config(load_json(args.config))
	paths = config['paths']
	work_dir = Path(paths['work_dir'])
	output_dir = Path(paths['output_dir'])
	work_dir.mkdir(parents=True, exist_ok=True)
	output_dir.mkdir(parents=True, exist_ok=True)

	sis_raw = work_dir / 'sis_raw.geojson'
	download(config['sources']['sis']['url'], sis_raw)
	sis = load_prepare_surfaces(sis_raw, config, config['sources']['sis'])

	centerline_index = None
	if config['sources'].get('centerlines', {}).get('enabled', False):
		centerlines_raw = work_dir / 'centerlines_input.geojsonseq'
		download(config['sources']['centerlines']['url'], centerlines_raw)
		centerlines = load_prepare_lines(centerlines_raw, config, config['sources']['centerlines'], config['crs']['metric'])
		centerline_index = make_direction_index(
			centerlines,
			'centerline',
			config['direction']['centerline_search_m'],
			('centerline_id', 'fmzk_id', 'FMZK_ID', 'OBJECTID', 'ID')
		)
	else:
		centerlines = gpd.GeoDataFrame(geometry=[], crs=config['crs']['metric'])

	street_index = None
	if config['sources'].get('street_axis', {}).get('enabled', False):
		street_raw = work_dir / 'street_axis_input.fgb'
		download(config['sources']['street_axis']['url'], street_raw)
		streets = load_prepare_lines(street_raw, config, config['sources']['street_axis'], config['crs']['input'])
		streets = filter_street_axis(streets, config['sources']['street_axis'])
		street_index = make_direction_index(
			streets,
			'street_axis',
			config['direction']['street_search_m'],
			('GIP_OBJECTID', 'OBJECTID', 'ID', 'fid', 'FID')
		)
	else:
		streets = gpd.GeoDataFrame(geometry=[], crs=config['crs']['metric'])

	summary = {
		'method': config['method'],
		'sis_features': int(len(sis)),
		'centerline_features': int(len(centerlines)),
		'street_axis_features': int(len(streets)),
		'processed': 0,
		'skipped': 0,
		'direction_centerline': 0,
		'direction_street_axis': 0,
		'direction_polygon_axis': 0,
		'direction_longest_edge': 0,
		'width_below_min': 0,
		'width_above_max': 0,
		'width_below_min_by_source': source_count_summary(),
		'width_above_max_by_source': source_count_summary(),
		'errors': {}
	}

	rows = []
	for idx, row in sis.iterrows():
		geom = row.geometry
		if geom is None or geom.is_empty:
			summary['skipped'] += 1
			continue

		props = row_props(row)
		source_id = source_feature_id(props, idx)
		source_type = props.get('TYPE') or props.get('type')

		try:
			area_m2 = float(geom.area)
			perimeter_m = float(geom.length)
			polygon_axis_angle = minimum_rectangle_axis_angle(geom)
			direction = choose_direction(geom, centerline_index, street_index, config, polygon_axis_angle)
			coords = polygon_coords(geom)
			width_raw_m = directional_width(coords, direction['angle_rad'])
			width_m = round(width_raw_m, 2)
			width_cls = width_class(width_m, config['width_classes'])
		except Exception as exc:
			summary['skipped'] += 1
			add_error(summary, exc)
			continue

		if width_m < config['direction']['min_width_m']:
			summary['width_below_min'] += 1
			increment_source_count(summary['width_below_min_by_source'], direction['source'])
		if width_m > config['direction']['max_width_m']:
			summary['width_above_max'] += 1
			increment_source_count(summary['width_above_max_by_source'], direction['source'])

		summary_key = f"direction_{direction['source']}"
		if summary_key in summary:
			summary[summary_key] += 1

		rows.append({
			'geometry': geom,
			'method': config['method'],
			'area_m2': round(area_m2, 2),
			'perimeter_m': round(perimeter_m, 2),
			'polygon_axis_angle_deg': angle_deg_or_none(polygon_axis_angle, 2),
			'width_m': width_m,
			'width_class': width_cls,
			'direction_source': direction['source'],
			'direction_angle_deg': round(math.degrees(direction['angle_rad']), 2),
			'direction_distance_m': round_or_none(direction['distance_m'], 2),
			'street_axis_angle_deg': angle_deg_or_none(direction.get('street_axis_angle_rad'), 2),
			'angle_diff_deg': round_or_none(direction.get('angle_diff_deg'), 2),
			'source_type': source_type,
			'source_id': source_id
		})
		summary['processed'] += 1

	output = output_geodataframe(rows, config['crs']['metric'])
	write_geojsonseq_gdf(paths['surfaces_metric'], output)
	write_geojsonseq_gdf(paths['surfaces_wgs84'], output.to_crs(config['crs']['pmtiles']))

	Path(paths['summary']).parent.mkdir(parents=True, exist_ok=True)
	Path(paths['summary']).write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
	print(json.dumps(summary, ensure_ascii=False))


def normalize_config(config):
	config.setdefault('method', 'directional_tu_widths')
	config.setdefault('preprocess', {})
	config['preprocess'].setdefault('make_valid', True)
	config['preprocess'].setdefault('explode_multipart', True)
	config['preprocess'].setdefault('min_area_m2', 0.01)

	config.setdefault('direction', {})
	config['direction'].setdefault('centerline_search_m', 3.0)
	config['direction'].setdefault('street_search_m', 25.0)
	config['direction'].setdefault('min_width_m', 0.2)
	config['direction'].setdefault('max_width_m', 8.0)
	config['direction'].setdefault('tangent_delta_m', 2.0)

	config['paths'].setdefault('surfaces_wgs84', 'work/directional-tu-widths/directional_surfaces_4326.geojsonseq')
	return config


def load_json(path):
	return json.loads(Path(path).read_text(encoding='utf-8'))


def download(url, target):
	target = Path(target)
	target.parent.mkdir(parents=True, exist_ok=True)
	url = substitute_env(str(url))

	if url.startswith('http://') or url.startswith('https://'):
		urlretrieve(url, target)
	else:
		target.write_bytes(Path(url).read_bytes())


def substitute_env(value):
	for key, env_value in os.environ.items():
		value = value.replace('${' + key + '}', env_value)
	return os.path.expandvars(value)


def load_prepare_surfaces(path, config, source_cfg):
	gdf = read_geodataframe(path, config, source_cfg, config['crs']['input'])
	gdf = apply_property_filter(gdf, source_cfg)

	if config['preprocess'].get('make_valid', True):
		gdf['geometry'] = gdf.geometry.apply(valid_or_empty)

	gdf = gdf[(gdf.geometry.notna()) & (~gdf.geometry.is_empty)].copy()

	if config['preprocess'].get('explode_multipart', True):
		gdf = gdf.explode(index_parts=False).reset_index(drop=True)

	gdf = gdf[gdf.geometry.geom_type.isin(['Polygon', 'MultiPolygon'])].copy()

	min_area = float(config['preprocess'].get('min_area_m2', 0.0))
	if min_area > 0:
		gdf = gdf[gdf.geometry.area >= min_area].copy()

	return gdf.reset_index(drop=True)


def load_prepare_lines(path, config, source_cfg, assumed_crs):
	gdf = read_geodataframe(path, config, source_cfg, assumed_crs)
	gdf = apply_property_filter(gdf, source_cfg)
	gdf = gdf[(gdf.geometry.notna()) & (~gdf.geometry.is_empty)].copy()
	gdf = explode_lines_gdf(gdf, config['crs']['metric'])
	return gdf.reset_index(drop=True)


def read_geodataframe(path, config, source_cfg, assumed_crs):
	try:
		gdf = gpd.read_file(path)
	except Exception as exc:
		if Path(path).suffix.lower() != '.fgb':
			raise
		converted = Path(path).with_suffix('.geojsonseq')
		try:
			run_ogr2ogr(path, converted)
		except FileNotFoundError as ogr_exc:
			raise RuntimeError('ogr2ogr is required when GeoPandas cannot read the FGB source directly') from ogr_exc
		except subprocess.CalledProcessError as ogr_exc:
			raise RuntimeError('ogr2ogr failed while converting the FGB source') from ogr_exc
		try:
			gdf = gpd.read_file(converted)
		except Exception:
			raise exc

	source_crs = source_cfg.get('crs') or assumed_crs
	if source_cfg.get('crs'):
		gdf = gdf.set_crs(source_crs, allow_override=True)
	elif gdf.crs is None:
		gdf = gdf.set_crs(source_crs)

	return gdf.to_crs(config['crs']['metric'])


def run_ogr2ogr(source, target):
	subprocess.run(
		['ogr2ogr', '-f', 'GeoJSONSeq', str(target), str(source)],
		check=True
	)


def apply_property_filter(gdf, source_cfg):
	prop = source_cfg.get('filter_property')
	values = source_cfg.get('filter_values') or []
	if prop and values and prop in gdf.columns:
		return gdf[gdf[prop].isin(values)].copy()
	return gdf


def filter_street_axis(gdf, source_cfg):
	filter_cfg = source_cfg.get('filter', {})

	baustatus_property = filter_cfg.get('baustatus_property')
	baustatus_values = filter_cfg.get('baustatus_values')
	if baustatus_property and baustatus_values and baustatus_property in gdf.columns:
		gdf = gdf[gdf[baustatus_property].isin(baustatus_values)].copy()

	frc_property = filter_cfg.get('frc_property')
	frc_include = filter_cfg.get('frc_include')
	if frc_property and frc_include and frc_property in gdf.columns:
		gdf = gdf[gdf[frc_property].isin(frc_include)].copy()

	fow_property = filter_cfg.get('fow_property')
	fow_exclude = filter_cfg.get('fow_exclude')
	if fow_property and fow_exclude and fow_property in gdf.columns:
		gdf = gdf[~gdf[fow_property].isin(fow_exclude)].copy()

	return gdf.reset_index(drop=True)


def valid_or_empty(geom):
	if geom is None or geom.is_empty:
		return geom
	return make_valid(geom)


def explode_lines_gdf(gdf, crs):
	rows = []
	for idx, row in gdf.iterrows():
		props = row_props(row)
		for line in explode_lines(row.geometry):
			if line.length <= 0:
				continue
			item = dict(props)
			item['geometry'] = line
			item['_source_index'] = idx
			rows.append(item)

	if not rows:
		return gpd.GeoDataFrame(geometry=[], crs=crs)

	return gpd.GeoDataFrame(rows, geometry='geometry', crs=crs)


def explode_lines(geom):
	if geom is None or geom.is_empty:
		return []
	if isinstance(geom, LineString):
		return [geom]
	if isinstance(geom, MultiLineString):
		return list(geom.geoms)
	if hasattr(geom, 'geoms'):
		out = []
		for part in geom.geoms:
			out.extend(explode_lines(part))
		return out
	return []


def make_direction_index(gdf, source_name, search_m, id_keys):
	if gdf is None or len(gdf) == 0:
		return None

	gdf = gdf.reset_index(drop=True)
	return {
		'gdf': gdf,
		'sindex': gdf.sindex,
		'source_name': source_name,
		'search_m': float(search_m),
		'id_keys': id_keys
	}


def choose_direction(geom, centerline_index, street_index, config, polygon_axis_angle):
	if centerline_index is not None:
		match = nearest_direction(geom, centerline_index, config, polygon_axis_angle)
		if match is not None:
			return match

	if street_index is not None:
		match = nearest_direction(geom, street_index, config, polygon_axis_angle)
		if match is not None:
			return match

	return fallback_direction(geom, polygon_axis_angle)


def nearest_direction(geom, direction_index, config, polygon_axis_angle=None):
	search_geom = geom.buffer(direction_index['search_m'])
	candidate_indices = query_spatial_index(direction_index['sindex'], search_geom)
	if len(candidate_indices) == 0:
		return None

	best = None
	for raw_idx in candidate_indices:
		idx = int(raw_idx)
		row = direction_index['gdf'].iloc[idx]
		line = row.geometry
		distance = geom.distance(line)
		if distance > direction_index['search_m']:
			continue

		angle = line_direction_near(line, geom, config['direction']['tangent_delta_m'])
		if angle is None:
			continue

		if best is None or distance < best['distance_m']:
			best = {
				'source': direction_index['source_name'],
				'angle_rad': angle,
				'distance_m': distance,
				'source_id': direction_source_id(row_props(row), idx, direction_index['id_keys'])
			}
			if direction_index['source_name'] == 'street_axis':
				best['street_axis_angle_rad'] = angle
				if polygon_axis_angle is not None:
					best['angle_diff_deg'] = axis_angle_diff_deg(angle, polygon_axis_angle)

	return best


def query_spatial_index(sindex, geom):
	try:
		return list(sindex.query(geom, predicate='intersects'))
	except TypeError:
		return list(sindex.query(geom))


def line_direction_near(line, geom, delta):
	if line.length <= 0:
		return None

	_, point_on_line = nearest_points(geom.representative_point(), line)
	s = line.project(point_on_line)
	tangent = tangent_at(line, s, delta)
	if tangent is None:
		return None

	return normalize_axis_angle(math.atan2(tangent[1], tangent[0]))


def tangent_at(line, s, delta):
	length = line.length
	if length <= 0:
		return None

	a = max(0.0, s - delta)
	b = min(length, s + delta)
	if b - a <= 1e-9:
		a = 0.0
		b = length
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


def fallback_direction(geom, polygon_axis_angle=None):
	angle = polygon_axis_angle
	if angle is None:
		angle = minimum_rectangle_axis_angle(geom)
	if angle is not None:
		return {
			'source': 'polygon_axis',
			'angle_rad': angle,
			'distance_m': 0.0,
			'source_id': None
		}

	angle = longest_edge_angle(geom)
	if angle is None:
		angle = 0.0

	return {
		'source': 'longest_edge',
		'angle_rad': angle,
		'distance_m': 0.0,
		'source_id': None
	}


def minimum_rectangle_axis_angle(geom):
	try:
		rect = geom.minimum_rotated_rectangle
	except Exception:
		return None

	if rect is None or rect.is_empty:
		return None

	if isinstance(rect, LineString):
		return line_angle(rect)

	if rect.geom_type != 'Polygon':
		return None

	best = None
	for a, b in pairwise(list(rect.exterior.coords)):
		length = math.hypot(b[0] - a[0], b[1] - a[1])
		if length <= 1e-9:
			continue
		if best is None or length > best[0]:
			best = (length, math.atan2(b[1] - a[1], b[0] - a[0]))

	if best is None:
		return None

	return normalize_axis_angle(best[1])


def longest_edge_angle(geom):
	best = None
	for poly in iter_polygons(geom):
		for a, b in pairwise(list(poly.exterior.coords)):
			length = math.hypot(b[0] - a[0], b[1] - a[1])
			if length <= 1e-9:
				continue
			if best is None or length > best[0]:
				best = (length, math.atan2(b[1] - a[1], b[0] - a[0]))

	if best is None:
		return None

	return normalize_axis_angle(best[1])


def line_angle(line):
	coords = list(line.coords)
	if len(coords) < 2:
		return None
	a = coords[0]
	b = coords[-1]
	if math.hypot(b[0] - a[0], b[1] - a[1]) <= 1e-9:
		return None
	return normalize_axis_angle(math.atan2(b[1] - a[1], b[0] - a[0]))


def directional_width(coords, direction_angle_rad):
	normal_angle = direction_angle_rad + math.pi / 2
	nx = math.cos(normal_angle)
	ny = math.sin(normal_angle)

	min_projection = math.inf
	max_projection = -math.inf

	for coord in coords:
		x = float(coord[0])
		y = float(coord[1])
		projection = x * nx + y * ny
		if projection < min_projection:
			min_projection = projection
		if projection > max_projection:
			max_projection = projection

	return max_projection - min_projection


def normalize_axis_angle(angle):
	angle = angle % math.pi
	if angle < 0:
		angle += math.pi
	return angle


def axis_angle_diff_deg(a, b):
	diff = abs(normalize_axis_angle(a) - normalize_axis_angle(b))
	if diff > math.pi / 2:
		diff = math.pi - diff
	return math.degrees(diff)


def angle_deg_or_none(value, digits):
	if value is None:
		return None
	return round(math.degrees(float(value)), digits)


def polygon_coords(geom):
	coords = []
	for poly in iter_polygons(geom):
		coords.extend(list(poly.exterior.coords))
	if not coords:
		raise ValueError('polygon has no coordinates')
	return coords


def iter_polygons(geom):
	if geom is None or geom.is_empty:
		return
	if geom.geom_type == 'Polygon':
		yield geom
	elif geom.geom_type == 'MultiPolygon':
		for part in geom.geoms:
			yield part
	elif hasattr(geom, 'geoms'):
		for part in geom.geoms:
			yield from iter_polygons(part)


def pairwise(coords):
	for i in range(len(coords) - 1):
		yield coords[i], coords[i + 1]


def width_class(width, classes):
	for cls in classes:
		min_v = cls.get('min')
		max_v = cls.get('max')
		if (min_v is None or width >= min_v) and (max_v is None or width < max_v):
			return cls['id']
	return None


def source_count_summary():
	return {
		'centerline': 0,
		'street_axis': 0,
		'polygon_axis': 0,
		'longest_edge': 0
	}


def increment_source_count(counts, source):
	counts[source] = counts.get(source, 0) + 1


def source_feature_id(props, fallback):
	for key in ('SIS_ID', 'OBJECTID', 'ID', 'fid', 'FID'):
		value = props.get(key)
		if value is not None:
			return value
	return fallback


def direction_source_id(props, fallback, keys):
	for key in keys:
		value = props.get(key)
		if value is not None:
			part_index = props.get('part_index')
			if part_index is not None and key in ('centerline_id', 'fmzk_id', 'FMZK_ID'):
				return f"{value}:{part_index}"
			return value
	return fallback


def row_props(row):
	return {
		key: json_value(value)
		for key, value in row.drop(labels=['geometry']).to_dict().items()
		if not key.startswith('_')
	}


def json_value(value):
	if value is None:
		return None

	try:
		if value != value:
			return None
	except Exception:
		pass

	if hasattr(value, 'item'):
		try:
			value = value.item()
		except Exception:
			pass

	if isinstance(value, (str, bool, int)):
		return value

	if isinstance(value, float):
		if not math.isfinite(value):
			return None
		return value

	return str(value)


def round_or_none(value, digits):
	if value is None:
		return None
	return round(float(value), digits)


def add_error(summary, exc):
	key = type(exc).__name__
	summary['errors'][key] = summary['errors'].get(key, 0) + 1


def output_geodataframe(rows, crs):
	if rows:
		return gpd.GeoDataFrame(rows, geometry='geometry', crs=crs)

	return gpd.GeoDataFrame(
		columns=[
			'method',
			'area_m2',
			'perimeter_m',
			'polygon_axis_angle_deg',
			'width_m',
			'width_class',
			'direction_source',
			'direction_angle_deg',
			'direction_distance_m',
			'street_axis_angle_deg',
			'angle_diff_deg',
			'source_type',
			'source_id',
			'geometry'
		],
		geometry='geometry',
		crs=crs
	)


def write_geojsonseq_gdf(path, gdf):
	Path(path).parent.mkdir(parents=True, exist_ok=True)

	with open(path, 'w', encoding='utf-8') as f:
		for _, row in gdf.iterrows():
			geom = row.geometry
			if geom is None or geom.is_empty:
				continue
			props = {
				key: json_value(row[key])
				for key in gdf.columns
				if key != 'geometry'
			}
			feature = {
				'type': 'Feature',
				'geometry': mapping(geom),
				'properties': props
			}
			f.write(json.dumps(feature, ensure_ascii=False))
			f.write('\n')


if __name__ == '__main__':
	main()
