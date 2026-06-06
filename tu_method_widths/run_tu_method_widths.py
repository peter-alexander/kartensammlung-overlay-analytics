#!/usr/bin/env python3
import argparse
import json
import math
import random
from pathlib import Path
from urllib.request import urlretrieve

import geopandas as gpd
from shapely import make_valid
from shapely.geometry import Point
from shapely.ops import unary_union


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

	surface_features = []
	street_features = []
	summary = {
		'method': config['method'],
		'sis_features': int(len(sis)),
		'processed': 0,
		'skipped': 0,
		'method_1_ok': 0,
		'method_2_ok': 0,
		'method_3_ok': 0,
		'street_axis_loaded': 0,
		'street_axis_filtered': 0,
		'street_axis_processed': 0,
		'street_axis_skipped': 0,
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

		surface_type = props.get('TYPE')
		
		try:
			d_c = minimum_enclosing_circle_diameter(geom)
			w1 = tu_method_1_width(surface_type, area_m2, d_c)
			summary['method_1_ok'] += 1
		except Exception as exc:
			d_c = None
			w1 = None
			add_error(summary, exc)
		
		try:
			d_i, m2_point = method_2_largest_inscribed_circle_approx(geom, config)
			w2 = tu_method_2_width(surface_type, area_m2, d_i)
			summary['method_2_ok'] += 1
		except Exception as exc:
			d_i = None
			w2 = None
			m2_point = None
			add_error(summary, exc)

		w_tu = w2
		width_class_m1 = width_class(w1, config['width_classes']) if w1 is not None else None
		width_class_m2 = width_class(w2, config['width_classes']) if w2 is not None else None
		width_class_tu = width_class_m2
		
		if w1 is not None and w2 is not None:
			m1_m2_diff_m = abs(w1 - w2)
		else:
			m1_m2_diff_m = None
		
		if width_class_m1 is not None and width_class_m2 is not None:
			m1_m2_class_diff = abs(width_class_m1 - width_class_m2)
		else:
			m1_m2_class_diff = None
		
		if m1_m2_class_diff is None:
			tu_quality = 'unknown'
		elif m1_m2_class_diff == 0:
			tu_quality = 'ok'
		elif m1_m2_class_diff == 1:
			tu_quality = 'uncertain'
		else:
			tu_quality = 'conflict'

		surface_features.append({
			'type': 'Feature',
			'geometry': geom.__geo_interface__,
			'properties': {
				'sis_id': props.get('SIS_ID') or props.get('OBJECTID') or props.get('ID') or idx,
				'type': props.get('TYPE'),
				'type_txt': props.get('TYPE_TXT'),
				'area_m2': round(area_m2, 2),
				'perimeter_m': round(perimeter_m, 2),
				'd_c': round_or_none(d_c, 2),
				'd_i': round_or_none(d_i, 2),
				'd_c': round_or_none(d_c, 2),
				'd_i': round_or_none(d_i, 2),
				
				'w_m1': round_or_none(w1, 2),
				'w_m2': round_or_none(w2, 2),
				
				'width_class_m1': width_class_m1,
				'width_class_m2': width_class_m2,
				
				'w_tu': round_or_none(w_tu, 2),
				'width_class_tu': width_class_tu,
				'tu_method': 'm2',
				
				'm1_m2_diff_m': round_or_none(m1_m2_diff_m, 2),
				'm1_m2_class_diff': m1_m2_class_diff,
				'tu_quality': tu_quality,
				
				'method': config['method']
			}
		})
		summary['processed'] += 1

	if config.get('method_3', {}).get('enabled', False) and config['sources'].get('street_axis', {}).get('enabled', False):
		try:
			street_features = build_method_3_street_axis_features(sis, config, work_dir, summary)
		except Exception as exc:
			add_error(summary, exc)

	write_geojsonseq(paths['surfaces_metric'], surface_features)
	write_geojsonseq(paths['streets_metric'], street_features)

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


def load_street_axis(path, config, source_cfg):
	gdf = gpd.read_file(path)

	if gdf.crs is None:
		gdf = gdf.set_crs(config['crs']['input'])

	gdf = gdf.to_crs(config['crs']['metric'])
	gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()

	filter_cfg = source_cfg.get('filter', {})

	baustatus_property = filter_cfg.get('baustatus_property')
	baustatus_values = filter_cfg.get('baustatus_values')
	if baustatus_property and baustatus_values and baustatus_property in gdf.columns:
		gdf = gdf[gdf[baustatus_property].isin(baustatus_values)].copy()

	frc_property = filter_cfg.get('frc_property', 'FRC')
	frc_include = filter_cfg.get('frc_include')
	if frc_property and frc_include and frc_property in gdf.columns:
		gdf = gdf[gdf[frc_property].isin(frc_include)].copy()

	fow_property = filter_cfg.get('fow_property', 'FOW')
	fow_exclude = filter_cfg.get('fow_exclude')
	if fow_property and fow_exclude and fow_property in gdf.columns:
		gdf = gdf[~gdf[fow_property].isin(fow_exclude)].copy()

	return gdf.reset_index(drop=True)


def build_method_3_street_axis_features(sis, config, work_dir, summary):
	source_cfg = config['sources']['street_axis']
	street_raw = work_dir / 'street_axis_input.fgb'
	download(source_cfg['url'], street_raw)

	streets = load_street_axis(street_raw, config, source_cfg)
	summary['street_axis_loaded'] = int(len(streets))
	summary['street_axis_filtered'] = int(len(streets))

	if len(streets) == 0 or len(sis) == 0:
		return []

	sis_union = unary_union(list(sis.geometry))
	buffer_m = float(config.get('method_3', {}).get('buffer_m', 20.0))

	features = []

	for idx, row in streets.iterrows():
		geom = row.geometry
		if geom is None or geom.is_empty:
			summary['street_axis_skipped'] += 1
			continue

		length_m = float(geom.length)
		if length_m <= 0:
			summary['street_axis_skipped'] += 1
			continue

		buffer_geom = geom.buffer(buffer_m)
		intersection = sis_union.intersection(buffer_geom)
		area_sidewalk_m2 = float(intersection.area) if not intersection.is_empty else 0.0

		w3 = 0.5 * area_sidewalk_m2 / (length_m + 2.0 * buffer_m)

		props = row.drop(labels=['geometry']).to_dict()

		features.append({
			'type': 'Feature',
			'geometry': geom.__geo_interface__,
			'properties': {
				'street_id': props.get('GIP_OBJECTID') or props.get('OBJECTID') or props.get('ID') or idx,
				'name': props.get('FEATURENAME') or props.get('MAINNAMETEXT') or props.get('NAME'),
				'frc': props.get('FRC'),
				'fow': props.get('FOW'),
				'baustatus': props.get('BAUSTATUS'),
				'axis_length_m': round(length_m, 2),
				'buffer_m': round(buffer_m, 2),
				'area_sidewalk_m2': round(area_sidewalk_m2, 2),
				'w_m3': round(w3, 2),
				'width_class_m3': width_class(w3, config['width_classes']),
				'method': config['method']
			}
		})

		summary['method_3_ok'] += 1
		summary['street_axis_processed'] += 1

	return features


def minimum_enclosing_circle_diameter(geom):
	points = []

	for x, y in geom.convex_hull.exterior.coords:
		points.append((float(x), float(y)))

	if not points:
		return None

	circle = make_minimum_enclosing_circle(points)

	if circle is None:
		return None

	return 2.0 * circle[2]


def make_minimum_enclosing_circle(points):
	points = list(dict.fromkeys(points))
	random.seed(1)
	random.shuffle(points)

	circle = None

	for i, p in enumerate(points):
		if circle is None or not point_in_circle(p, circle):
			circle = make_circle_one_point(points[:i + 1], p)

	return circle


def make_circle_one_point(points, p):
	circle = (p[0], p[1], 0.0)

	for i, q in enumerate(points):
		if not point_in_circle(q, circle):
			if circle[2] == 0.0:
				circle = make_diameter_circle(p, q)
			else:
				circle = make_circle_two_points(points[:i + 1], p, q)

	return circle


def make_circle_two_points(points, p, q):
	circle = make_diameter_circle(p, q)
	left = None
	right = None

	px, py = p
	qx, qy = q

	for r in points:
		if point_in_circle(r, circle):
			continue

		cross = cross_product(px, py, qx, qy, r[0], r[1])
		c = make_circumcircle(p, q, r)

		if c is None:
			continue

		cross_center = cross_product(px, py, qx, qy, c[0], c[1])

		if cross > 0 and (left is None or cross_center > cross_product(px, py, qx, qy, left[0], left[1])):
			left = c
		elif cross < 0 and (right is None or cross_center < cross_product(px, py, qx, qy, right[0], right[1])):
			right = c

	if left is None and right is None:
		return circle
	if left is None:
		return right
	if right is None:
		return left

	return left if left[2] <= right[2] else right


def make_diameter_circle(p, q):
	cx = (p[0] + q[0]) / 2.0
	cy = (p[1] + q[1]) / 2.0
	r = math.hypot(cx - p[0], cy - p[1])
	return (cx, cy, r)


def make_circumcircle(p, q, r):
	ax, ay = p
	bx, by = q
	cx, cy = r

	d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))

	if abs(d) < 1e-12:
		return None

	ux = (
		(ax * ax + ay * ay) * (by - cy)
		+ (bx * bx + by * by) * (cy - ay)
		+ (cx * cx + cy * cy) * (ay - by)
	) / d

	uy = (
		(ax * ax + ay * ay) * (cx - bx)
		+ (bx * bx + by * by) * (ax - cx)
		+ (cx * cx + cy * cy) * (bx - ax)
	) / d

	radius = math.hypot(ux - ax, uy - ay)
	return (ux, uy, radius)


def point_in_circle(p, circle):
	eps = 1e-9
	return math.hypot(p[0] - circle[0], p[1] - circle[1]) <= circle[2] + eps


def cross_product(ax, ay, bx, by, cx, cy):
	return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def tu_method_1_width(surface_type, area_m2, d_c):
	if d_c is None or d_c <= 0:
		return None

	if surface_type == 'GG':
		return area_m2 / d_c

	if surface_type in ('EE', 'HH'):
		return d_c

	return area_m2 / d_c


def tu_method_2_width(surface_type, area_m2, d_i):
	if d_i is None or d_i <= 0:
		return None

	if surface_type == 'GG':
		return d_i

	if surface_type in ('EE', 'HH'):
		return area_m2 / d_i

	return d_i

def method_2_largest_inscribed_circle_approx(geom, config):
	# Approximation of largest inscribed circle:
	# iterative grid search for the point with maximum distance to polygon boundary.
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
