#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path
from urllib.request import urlretrieve

import geopandas as gpd
from shapely import make_valid
from shapely.geometry import Point


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', default='config/tu-method-widths.json')
	args = parser.parse_args()

	config = load_json(args.config)
	paths = config['paths']
	work_dir = Path(paths['work_dir'])
	output_dir = Path(paths['output_dir'])
	work_dir.mkdir(parents=True, exist_ok=True)
	output_dir.mkdir(parents=True, exist_ok=True)

	sis_raw = work_dir / 'sis_raw.geojson'
	download(config['sources']['sis']['url'], sis_raw)
	sis = load_prepare_source(sis_raw, config, config['sources']['sis'])

	features = []
	summary = {
		'method': config['method'],
		'sis_features': int(len(sis)),
		'processed': 0,
		'skipped': 0,
		'method_1_ok': 0,
		'method_2_ok': 0,
		'method_3_ok': 0,
		'errors': {}
	}

	for idx, row in sis.iterrows():
		geom = row.geometry
		if geom is None or geom.is_empty:
			summary['skipped'] += 1
			continue
		props = row.drop(labels=['geometry']).to_dict()
		area_m2 = float(geom.area)
		perimeter_m = float(geom.length)
		if area_m2 < config['preprocess'].get('min_area_m2', 0.0):
			summary['skipped'] += 1
			continue

		try:
			w1 = method_1_circle_from_area(area_m2)
			summary['method_1_ok'] += 1
		except Exception as exc:
			w1 = None
			add_error(summary, exc)

		try:
			w2, m2_point = method_2_largest_inscribed_circle_approx(geom, config)
			summary['method_2_ok'] += 1
		except Exception as exc:
			w2 = None
			m2_point = None
			add_error(summary, exc)

		features.append({
			'type': 'Feature',
			'geometry': geom.__geo_interface__,
			'properties': {
				'sis_id': props.get('SIS_ID') or props.get('OBJECTID') or props.get('ID') or idx,
				'type': props.get('TYPE'),
				'type_txt': props.get('TYPE_TXT'),
				'area_m2': round(area_m2, 2),
				'perimeter_m': round(perimeter_m, 2),
				'w_m1': round_or_none(w1, 2),
				'w_m2': round_or_none(w2, 2),
				'width_class_m1': width_class(w1, config['width_classes']) if w1 is not None else None,
				'width_class_m2': width_class(w2, config['width_classes']) if w2 is not None else None,
				'method': config['method']
			}
		})
		summary['processed'] += 1

	write_geojsonseq(paths['surfaces_metric'], features)
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

	if config['preprocess'].get('make_valid', True):
		gdf['geometry'] = gdf.geometry.apply(lambda g: make_valid(g) if g is not None and not g.is_empty else g)

	gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()
	if config['preprocess'].get('explode_multipart', True):
		gdf = gdf.explode(index_parts=False).reset_index(drop=True)
	return gdf


def method_1_circle_from_area(area_m2):
	# TU method 1 approximation: sidewalk width as diameter of an equal-area circle.
	return 2.0 * math.sqrt(area_m2 / math.pi)


def method_2_largest_inscribed_circle_approx(geom, config):
	# Approximation of largest inscribed circle: iterative grid search for the point
	# with maximum distance to the polygon boundary while remaining inside the polygon.
	minx, miny, maxx, maxy = geom.bounds
	grid_nodes = int(config['method_2'].get('grid_nodes', 21))
	iterations = int(config['method_2'].get('iterations', 8))
	shrink = float(config['method_2'].get('shrink_factor', math.sqrt(2)))

	best_point = geom.representative_point()
	best_dist = best_point.distance(geom.boundary)
	x0, y0, x1, y1 = minx, miny, maxx, maxy

	for _ in range(iterations):
		if x1 <= x0 or y1 <= y0:
			break
		for i in range(grid_nodes):
			x = x0 + (x1 - x0) * i / max(1, grid_nodes - 1)
			for j in range(grid_nodes):
				y = y0 + (y1 - y0) * j / max(1, grid_nodes - 1)
				p = Point(x, y)
				if not geom.covers(p):
					continue
				d = p.distance(geom.boundary)
				if d > best_dist:
					best_dist = d
					best_point = p

		w = (x1 - x0) / shrink
		h = (y1 - y0) / shrink
		x0 = best_point.x - w / 2
		x1 = best_point.x + w / 2
		y0 = best_point.y - h / 2
		y1 = best_point.y + h / 2

	return 2.0 * best_dist, best_point


def width_class(width, classes):
	for cls in classes:
		min_v = cls.get('min')
		max_v = cls.get('max')
		if (min_v is None or width >= min_v) and (max_v is None or width < max_v):
			return cls['id']
	return None


def round_or_none(value, digits):
	if value is None:
		return None
	return round(float(value), digits)


def add_error(summary, exc):
	key = type(exc).__name__
	summary['errors'][key] = summary['errors'].get(key, 0) + 1


def write_geojsonseq(path, features):
	Path(path).parent.mkdir(parents=True, exist_ok=True)
	with open(path, 'w', encoding='utf-8') as f:
		for feat in features:
			f.write(json.dumps(feat, ensure_ascii=False))
			f.write('\n')


if __name__ == '__main__':
	main()
