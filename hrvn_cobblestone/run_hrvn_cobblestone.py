#!/usr/bin/env python3
import argparse
import json
import math
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import geopandas as gpd
import requests
from requests.adapters import HTTPAdapter
from shapely import make_valid, union_all
from shapely.geometry import box, mapping, shape
from urllib3.util.retry import Retry


def main():
	parser = argparse.ArgumentParser(description='Find SIS cobblestone surfaces on Vienna main cycle routes.')
	parser.add_argument('--config', default='config/hrvn-cobblestone.json')
	parser.add_argument('--refresh', action='store_true', help='Redownload cached WFS tiles.')
	parser.add_argument('--limit-tiles', type=int, default=None, help='Process at most N tiles (debug/testing).')
	args = parser.parse_args()

	config = load_json(args.config)
	paths = config['paths']
	work_dir = Path(paths['work_dir'])
	tile_dir = work_dir / 'sis-tiles'
	output_dir = Path(paths['output_dir'])
	work_dir.mkdir(parents=True, exist_ok=True)
	tile_dir.mkdir(parents=True, exist_ok=True)
	output_dir.mkdir(parents=True, exist_ok=True)

	session = make_session(config['download'])
	geometry_property = discover_geometry_property(session, config)
	print(f'SIS geometry property: {geometry_property}')

	hrvn_path = work_dir / 'hrvn.geojson'
	download_file(session, config['sources']['hrvn']['url'], hrvn_path, refresh=args.refresh)
	hrvn = load_hrvn(hrvn_path, config)
	corridor = build_corridor(hrvn, config)

	tiles = build_tiles(corridor.bounds, config['download']['tile_size_m'])
	if args.limit_tiles is not None:
		tiles = tiles[:max(0, args.limit_tiles)]

	loaded_geometries = []
	tile_stats = []
	for index, tile in enumerate(tiles, start=1):
		tile_id = tile_id_for(tile)
		tile_path = tile_dir / f'{tile_id}.geojson'
		features, source = load_or_download_sis_tile(
			session=session,
			config=config,
			tile=tile,
			tile_path=tile_path,
			geometry_property=geometry_property,
			refresh=args.refresh
		)
		geometries = extract_cobblestone_geometries(features, config, tile)
		if geometries:
			loaded_geometries.append(union_all(geometries))
		tile_stats.append({
			'id': tile_id,
			'source': source,
			'features': len(features),
			'cobblestone_geometries': len(geometries)
		})
		print(f'[{index}/{len(tiles)}] {tile_id}: {len(features)} features, {len(geometries)} cobblestone geometries ({source})')

	if loaded_geometries:
		all_cobblestone = clean_polygonal(union_all(loaded_geometries))
	else:
		all_cobblestone = None

	if all_cobblestone is not None and not all_cobblestone.is_empty:
		hrvn_cobblestone = clean_polygonal(all_cobblestone.intersection(corridor))
	else:
		hrvn_cobblestone = None

	simplify_m = float(config['output']['simplify_tolerance_m'])
	debug_simplify_m = float(config['output'].get('debug_simplify_tolerance_m', simplify_m))
	all_cobblestone_simplified = simplify_polygonal(all_cobblestone, debug_simplify_m)
	hrvn_cobblestone_simplified = simplify_polygonal(hrvn_cobblestone, simplify_m)

	write_polygon_outputs(
		all_cobblestone_simplified,
		config['crs']['metric'],
		paths['debug_all_cobblestone_geojson'],
		paths.get('debug_all_cobblestone_fgb'),
		config['crs']['output'],
		'cobblestone_debug'
	)
	write_polygon_outputs(
		hrvn_cobblestone_simplified,
		config['crs']['metric'],
		paths['hrvn_cobblestone_geojson'],
		paths.get('hrvn_cobblestone_fgb'),
		config['crs']['output'],
		'hrvn_cobblestone'
	)
	write_polygon_outputs(
		corridor,
		config['crs']['metric'],
		paths['debug_hrvn_corridor_geojson'],
		None,
		config['crs']['output'],
		'hrvn_corridor'
	)

	summary = {
		'hrvn_features': int(len(hrvn)),
		'hrvn_ranks': config['sources']['hrvn']['filter_values'],
		'buffer_m_each_side': float(config['analysis']['buffer_m']),
		'sis_belags': config['sources']['sis']['filter_values'],
		'tile_size_m': float(config['download']['tile_size_m']),
		'tiles_processed': len(tiles),
		'corridor_area_m2': round(float(corridor.area), 2),
		'all_loaded_cobblestone_area_m2': area_or_zero(all_cobblestone_simplified),
		'hrvn_cobblestone_area_m2': area_or_zero(hrvn_cobblestone_simplified),
		'tiles': tile_stats
	}
	Path(paths['summary']).parent.mkdir(parents=True, exist_ok=True)
	Path(paths['summary']).write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
	print(json.dumps({key: value for key, value in summary.items() if key != 'tiles'}, ensure_ascii=False, indent=2))


def load_json(path):
	return json.loads(Path(path).read_text(encoding='utf-8'))


def make_session(download_cfg):
	retry = Retry(
		total=int(download_cfg.get('retries', 5)),
		connect=int(download_cfg.get('retries', 5)),
		read=int(download_cfg.get('retries', 5)),
		backoff_factor=float(download_cfg.get('retry_backoff_s', 1.0)),
		status_forcelist=(429, 500, 502, 503, 504),
		allowed_methods=frozenset(['GET'])
	)
	session = requests.Session()
	session.mount('https://', HTTPAdapter(max_retries=retry))
	session.headers.update({'User-Agent': 'kartensammlung-overlay-analytics/hrvn-cobblestone'})
	return session


def download_file(session, url, target, refresh=False):
	target = Path(target)
	if target.exists() and target.stat().st_size > 0 and not refresh:
		return 'cache'

	target.parent.mkdir(parents=True, exist_ok=True)
	response = session.get(url, timeout=120)
	response.raise_for_status()
	target.write_bytes(response.content)
	return 'download'


def load_hrvn(path, config):
	gdf = gpd.read_file(path)
	if gdf.crs is None:
		gdf = gdf.set_crs(config['crs']['input'])
	else:
		gdf = gdf.to_crs(config['crs']['input'])

	filter_property = config['sources']['hrvn']['filter_property']
	filter_values = {str(value) for value in config['sources']['hrvn']['filter_values']}
	if filter_property not in gdf.columns:
		raise RuntimeError(f'HRVN filter property {filter_property!r} is missing')

	gdf = gdf[gdf[filter_property].astype(str).isin(filter_values)].copy()
	gdf = gdf[(gdf.geometry.notna()) & (~gdf.geometry.is_empty)].copy()
	if gdf.empty:
		raise RuntimeError('HRVN filter returned no features')

	return gdf.to_crs(config['crs']['metric']).reset_index(drop=True)


def build_corridor(hrvn, config):
	buffer_m = float(config['analysis']['buffer_m'])
	buffered = [geom.buffer(buffer_m, cap_style='round', join_style='round') for geom in hrvn.geometry if geom is not None and not geom.is_empty]
	corridor = clean_polygonal(union_all(buffered))
	if corridor is None or corridor.is_empty:
		raise RuntimeError('HRVN corridor is empty')
	return corridor


def build_tiles(bounds, tile_size_m):
	minx, miny, maxx, maxy = bounds
	tile_size_m = float(tile_size_m)
	start_x = math.floor(minx / tile_size_m) * tile_size_m
	start_y = math.floor(miny / tile_size_m) * tile_size_m
	end_x = math.ceil(maxx / tile_size_m) * tile_size_m
	end_y = math.ceil(maxy / tile_size_m) * tile_size_m

	tiles = []
	x = start_x
	while x < end_x:
		y = start_y
		while y < end_y:
			tiles.append(box(x, y, min(x + tile_size_m, end_x), min(y + tile_size_m, end_y)))
			y += tile_size_m
		x += tile_size_m
	return tiles


def tile_id_for(tile):
	minx, miny, maxx, maxy = tile.bounds
	return f'{int(round(minx))}_{int(round(miny))}_{int(round(maxx))}_{int(round(maxy))}'


def discover_geometry_property(session, config):
	sis_cfg = config['sources']['sis']
	explicit = sis_cfg.get('geometry_property')
	if explicit:
		return explicit

	params = {
		'service': 'WFS',
		'version': sis_cfg.get('version', '1.1.0'),
		'request': 'DescribeFeatureType',
		'typeName': sis_cfg['type_name']
	}
	response = session.get(sis_cfg['base_url'], params=params, timeout=float(config['download'].get('timeout_s', 180)))
	response.raise_for_status()
	return geometry_property_from_xsd(response.content)


def geometry_property_from_xsd(xml_content):
	root = ET.fromstring(xml_content)
	candidates = []
	for element in root.iter():
		if not str(element.tag).endswith('element'):
			continue
		name = element.attrib.get('name')
		type_name = element.attrib.get('type', '')
		if not name:
			continue
		type_lower = type_name.lower()
		if 'gml:' in type_lower and any(token in type_lower for token in ('geometry', 'surface', 'polygon', 'curve', 'point', 'line')):
			candidates.append(name)

	if not candidates:
		raise RuntimeError('DescribeFeatureType did not expose a GML geometry property')
	return candidates[0]


def load_or_download_sis_tile(session, config, tile, tile_path, geometry_property, refresh=False):
	if tile_path.exists() and tile_path.stat().st_size > 0 and not refresh:
		return load_geojson_features(tile_path), 'cache'

	sis_cfg = config['sources']['sis']
	minx, miny, maxx, maxy = tile.bounds
	params = {
		'service': 'WFS',
		'version': sis_cfg.get('version', '1.1.0'),
		'request': 'GetFeature',
		'typeName': sis_cfg['type_name'],
		'srsName': config['crs']['metric'],
		'outputFormat': sis_cfg.get('output_format', 'json')
	}
	bbox_cql = (
		f"BBOX({geometry_property},{fmt(minx)},{fmt(miny)},{fmt(maxx)},{fmt(maxy)},"
		f"'{config['crs']['metric']}')"
	)
	if sis_cfg.get('server_filter', True):
		values = ','.join(f"'{escape_cql(value)}'" for value in sis_cfg['filter_values'])
		params['cql_filter'] = f"{bbox_cql} AND {sis_cfg['filter_property']} IN ({values})"
	else:
		params['cql_filter'] = bbox_cql

	response = session.get(sis_cfg['base_url'], params=params, timeout=float(config['download'].get('timeout_s', 180)))
	response.raise_for_status()
	max_bytes = int(float(config['download'].get('max_response_mb', 512)) * 1024 * 1024)
	if len(response.content) > max_bytes:
		raise RuntimeError(f'SIS tile {tile_id_for(tile)} exceeded max_response_mb ({len(response.content) / 1024 / 1024:.1f} MiB)')

	try:
		payload = response.json()
	except ValueError as exc:
		snippet = response.text[:500].replace('\n', ' ')
		raise RuntimeError(f'SIS tile {tile_id_for(tile)} did not return JSON: {snippet}') from exc

	if payload.get('type') != 'FeatureCollection' or not isinstance(payload.get('features'), list):
		raise RuntimeError(f'SIS tile {tile_id_for(tile)} returned unexpected payload')

	tile_path.parent.mkdir(parents=True, exist_ok=True)
	tile_path.write_text(json.dumps(payload, ensure_ascii=False), encoding='utf-8')
	delay_s = float(config['download'].get('delay_between_requests_s', 0.0))
	if delay_s > 0:
		time.sleep(delay_s)
	return payload['features'], 'download'


def load_geojson_features(path):
	payload = json.loads(Path(path).read_text(encoding='utf-8'))
	features = payload.get('features')
	if not isinstance(features, list):
		raise RuntimeError(f'Cached tile {path} has no feature list')
	return features


def extract_cobblestone_geometries(features, config, tile):
	sis_cfg = config['sources']['sis']
	property_name = sis_cfg['filter_property']
	allowed_values = {str(value) for value in sis_cfg['filter_values']}
	geometries = []

	for feature in features:
		properties = feature.get('properties') or {}
		if str(properties.get(property_name, '')) not in allowed_values:
			continue
		geometry_data = feature.get('geometry')
		if not geometry_data:
			continue
		try:
			geom = make_valid(shape(geometry_data))
		except Exception:
			continue
		geom = clean_polygonal(geom.intersection(tile))
		if geom is None or geom.is_empty:
			continue
		geometries.append(geom)

	return geometries


def clean_polygonal(geom):
	if geom is None or geom.is_empty:
		return None
	geom = make_valid(geom)
	if geom.geom_type in ('Polygon', 'MultiPolygon'):
		return geom
	if hasattr(geom, 'geoms'):
		parts = [part for part in geom.geoms if part.geom_type in ('Polygon', 'MultiPolygon') and not part.is_empty]
		if parts:
			return make_valid(union_all(parts))
	return None


def simplify_polygonal(geom, tolerance_m):
	if geom is None or geom.is_empty:
		return None
	if tolerance_m <= 0:
		return clean_polygonal(geom)
	return clean_polygonal(geom.simplify(float(tolerance_m), preserve_topology=True))


def write_polygon_outputs(geom, source_crs, geojson_path, fgb_path, output_crs, layer_name):
	gdf = geometry_to_gdf(geom, source_crs, layer_name)
	if output_crs and str(output_crs) != str(source_crs):
		gdf = gdf.to_crs(output_crs)

	geojson_path = Path(geojson_path)
	geojson_path.parent.mkdir(parents=True, exist_ok=True)
	write_geojson(geojson_path, gdf)

	if fgb_path:
		fgb_path = Path(fgb_path)
		fgb_path.parent.mkdir(parents=True, exist_ok=True)
		if fgb_path.exists():
			fgb_path.unlink()
		if not gdf.empty:
			gdf.to_file(fgb_path, driver='FlatGeobuf')


def geometry_to_gdf(geom, crs, layer_name):
	rows = []
	if geom is not None and not geom.is_empty:
		parts = list(geom.geoms) if geom.geom_type == 'MultiPolygon' else [geom]
		for idx, part in enumerate(parts, start=1):
			rows.append({
				'id': idx,
				'layer': layer_name,
				'area_m2': round(float(part.area), 2),
				'geometry': part
			})
	if rows:
		return gpd.GeoDataFrame(rows, geometry='geometry', crs=crs)
	return gpd.GeoDataFrame(columns=['id', 'layer', 'area_m2', 'geometry'], geometry='geometry', crs=crs)


def write_geojson(path, gdf):
	features = []
	for _, row in gdf.iterrows():
		features.append({
			'type': 'Feature',
			'geometry': mapping(row.geometry),
			'properties': {
				'id': int(row['id']),
				'layer': str(row['layer']),
				'area_m2': float(row['area_m2'])
			}
		})
	payload = {'type': 'FeatureCollection', 'features': features}
	path.write_text(json.dumps(payload, ensure_ascii=False, separators=(',', ':')), encoding='utf-8')


def area_or_zero(geom):
	if geom is None or geom.is_empty:
		return 0.0
	return round(float(geom.area), 2)


def fmt(value):
	return f'{float(value):.3f}'.rstrip('0').rstrip('.')


def escape_cql(value):
	return str(value).replace("'", "''")


if __name__ == '__main__':
	main()
