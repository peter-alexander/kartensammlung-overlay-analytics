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
			direction = choose_direction(geom, centerline_index, street_index, config)
			coords = polygon_coords(geom)
			width_m = directional_width(coords, direction['angle_rad'])
			width_cls = width_class(width_m, config['width_classes'])
		except Exception as exc:
			summary['skipped'] += 1
			add_error(summary, exc)
			continue

		if width_m < config['direction']['min_width_m']:
			summary['width_below_min'] += 1
		if width_m > config['direction']['max_width_m']:
			summary['width_above_max'] += 1

		summary_key = f"direction_{direction['source']}"
		if summary_key in summary:
			summary[summary_key] += 1

		rows.append({
			'geometry': geom,
			'method': config['method'],
			'width_m': round(width_m, 2),
			'width_class': width_cls,
			'direction_source': direction['source'],
			'direction_angle_deg': round(math.degrees(direction['angle_rad']), 2),
			'direction_distance_m': round_or_none(direction['distance_m'], 2),
			'source_type': source_type,
			'source_id': source_id
		})
		summary['processed'] += 1

	if rows:
		output = gpd.GeoDataFrame(rows, geometry='geometry', crs=config['crs']['metric'])
	else:
		output = gpd.GeoDataFrame(
			columns=[
				'method',
				'width_m',
				'width_class',
				'direction_source',
				'direction_angle_deg',
				'direction_distance_m',
				'source_type',
				'source_id',
				'geometry'
			],
			geometry='geometry',
			crs=config['crs']['metric']
		)
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
	except Exception:
		if Path(path).suffix.lower() != '.fgb':
			raise
		converted = Path(path).with_suffix('.geojsonseq')
		run_ogr2ogr(path, converted)
		gdf = gpd.read_file(converted)

	source_crs = source_cfg.get('crs') or assumed_crs
	if gdf.crs is None:
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
		gdf = gdfmÆÈ‹j◊ù