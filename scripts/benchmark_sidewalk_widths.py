#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import resource
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import requests
from shapely import LineString, Point, box, line_interpolate_point, make_valid
from shapely.geometry import GeometryCollection, LineString as GeoLineString
from shapely.geometry import MultiLineString, MultiPolygon, Polygon, mapping, shape
from shapely.prepared import prep

WFS_BASE = "https://data.wien.gv.at/daten/geo"
EPSG_METRIC = "EPSG:31256"
TIMEOUT = 180


@dataclass(slots=True)
class AoiDef:
	name: str
	bbox: tuple[float, float, float, float]
	clip_geom: Polygon | MultiPolygon


@dataclass(slots=True)
class MeasureConfig:
	step_m: float
	span_m: float
	mark_half_m: float
	vertex_margin_m: float
	tangent_window_m: float
	inner_probe_m: float
	min_width_m: float
	only_outer_rings: bool
	debug_feature_limit: int


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser()
	parser.add_argument("--config", required=True)
	parser.add_argument("--output", required=True)
	return parser.parse_args()


def load_config(path: str | Path) -> dict[str, Any]:
	with open(path, "r", encoding="utf-8") as f:
		return json.load(f)


def fetch_geojson(session: requests.Session, type_name: str, cql_filter: str | None = None) -> dict[str, Any]:
	params = {
		"service": "WFS",
		"version": "1.1.0",
		"request": "GetFeature",
		"typeName": type_name,
		"outputFormat": "json",
		"SRSNAME": EPSG_METRIC,
	}
	if cql_filter:
		params["CQL_FILTER"] = cql_filter

	res = session.get(WFS_BASE, params=params, timeout=TIMEOUT)
	res.raise_for_status()
	return res.json()


def feature_geometries(fc: dict[str, Any]) -> list[Any]:
	geoms: list[Any] = []
	for feat in fc.get("features", []):
		geom = feat.get("geometry")
		if not geom:
			continue
		geoms.append(shape(geom))
	return geoms


def extract_polygons(geom: Any) -> list[Polygon]:
	if geom is None or geom.is_empty:
		return []
	if isinstance(geom, Polygon):
		return [geom]
	if isinstance(geom, MultiPolygon):
		return [g for g in geom.geoms if not g.is_empty]
	if isinstance(geom, GeometryCollection):
		out: list[Polygon] = []
		for g in geom.geoms:
			out.extend(extract_polygons(g))
		return out
	return []


def extract_lines(geom: Any) -> list[LineString]:
	if geom is None or geom.is_empty:
		return []
	if isinstance(geom, LineString):
		return [geom]
	if isinstance(geom, MultiLineString):
		return [g for g in geom.geoms if not g.is_empty]
	if isinstance(geom, GeometryCollection):
		out: list[LineString] = []
		for g in geom.geoms:
			out.extend(extract_lines(g))
		return out
	return []


def segment_length(a: tuple[float, float], b: tuple[float, float]) -> float:
	dx = b[0] - a[0]
	dy = b[1] - a[1]
	return math.hypot(dx, dy)


def unit_vector(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
	length = segment_length(a, b)
	if length == 0:
		return (0.0, 0.0)
	return ((b[0] - a[0]) / length, (b[1] - a[1]) / length)


def ring_chain_distances(coords: list[tuple[float, float]]) -> list[float]:
	out = [0.0]
	for i in range(1, len(coords)):
		out.append(out[-1] + segment_length(coords[i - 1], coords[i]))
	return out


def tangent_at_distance(ring_line: LineString, chain_m: float, window_m: float, fallback: tuple[float, float]) -> tuple[float, float]:
	half = max(0.25, window_m * 0.5)
	length = ring_line.length
	d0 = max(0.0, chain_m - half)
	d1 = min(length, chain_m + half)
	if d1 - d0 < 0.05:
		return fallback

	p0 = line_interpolate_point(ring_line, d0)
	p1 = line_interpolate_point(ring_line, d1)
	dx = p1.x - p0.x
	dy = p1.y - p0.y
	vec_len = math.hypot(dx, dy)
	if vec_len < 1e-9:
		return fallback
	return (dx / vec_len, dy / vec_len)


def choose_inward_normal(prepared_poly: Any, px: float, py: float, tx: float, ty: float, probe_m: float) -> tuple[float, float] | None:
	nl = (-ty, tx)
	nr = (ty, -tx)

	for scale in (1.0, 2.0, 4.0):
		probe = probe_m * scale
		left_inside = prepared_poly.contains(Point(px + nl[0] * probe, py + nl[1] * probe))
		right_inside = prepared_poly.contains(Point(px + nr[0] * probe, py + nr[1] * probe))

		if left_inside and not right_inside:
			return nl
		if right_inside and not left_inside:
			return nr
		if left_inside and right_inside:
			return nl

	return None


def best_cross_segment(intersection_geom: Any, center: Point, min_width_m: float) -> LineString | None:
	segments = []
	for seg in extract_lines(intersection_geom):
		if seg.length >= min_width_m:
			segments.append(seg)

	if not segments:
		return None

	containing = [seg for seg in segments if seg.distance(center) <= 1e-7]
	if containing:
		return max(containing, key=lambda g: g.length)

	return min(segments, key=lambda g: g.distance(center))


def build_output_line(mx: float, my: float, tx: float, ty: float, half_len: float) -> LineString:
	return GeoLineString([
		(mx - tx * half_len, my - ty * half_len),
		(mx + tx * half_len, my + ty * half_len),
	])


def resolve_aois(session: requests.Session, config: dict[str, Any]) -> list[AoiDef]:
	out: list[AoiDef] = []
	for aoi in config.get("aois", []):
		mode = aoi["mode"]
		name = aoi["name"]

		if mode == "vienna_boundary":
			fc = fetch_geojson(session, "ogdwien:LANDESGRENZEOGD")
			geoms = feature_geometries(fc)
			polys: list[Polygon] = []
			for geom in geoms:
				valid = make_valid(geom) if not geom.is_valid else geom
				polys.extend(extract_polygons(valid))
			if not polys:
				raise RuntimeError("Keine Wien-Grenze gefunden")
			clip_geom = polys[0]
			for poly in polys[1:]:
				clip_geom = clip_geom.union(poly)
			bbox = clip_geom.bounds
			out.append(AoiDef(name=name, bbox=bbox, clip_geom=clip_geom))
			continue

		if mode == "bbox":
			bbox_list = aoi["bbox"]
			if len(bbox_list) != 4:
				raise ValueError(f"AOI {name}: bbox muss 4 Werte haben")
			bbox = tuple(float(v) for v in bbox_list)
			out.append(AoiDef(name=name, bbox=bbox, clip_geom=box(*bbox)))
			continue

		raise ValueError(f"Unbekannter AOI-Modus: {mode}")

	if not out:
		raise ValueError("Keine AOIs in der Config")
	return out


def fetch_sidewalk_geometries(session: requests.Session, aoi: AoiDef) -> list[Any]:
	minx, miny, maxx, maxy = aoi.bbox
	bbox_sql = f"BBOX(SHAPE, {minx},{miny},{maxx},{maxy}, '{EPSG_METRIC}')"
	cql = f"{bbox_sql} AND LAYER = 'Gehsteig'"
	fc = fetch_geojson(session, "ogdwien:FMZKVERKEHR2OGD", cql_filter=cql)
	return feature_geometries(fc)


def measure_sidewalks(geoms: Iterable[Any], aoi: AoiDef, cfg: MeasureConfig) -> tuple[list[dict[str, Any]], dict[str, Any]]:
	stats: dict[str, Any] = {
		"input_features": 0,
		"repaired_features": 0,
		"polygons_after_prepare": 0,
		"stations_attempted": 0,
		"stations_written": 0,
		"short_segments_skipped": 0,
		"failed_normal": 0,
		"empty_intersections": 0,
		"intersection_errors": 0,
	}
	out_features: list[dict[str, Any]] = []
	prepared_clip = prep(aoi.clip_geom)

	for geom in geoms:
		stats["input_features"] += 1
		work = geom

		if not work.is_valid:
			stats["repaired_features"] += 1
			work = make_valid(work)

		if work.is_empty:
			continue

		if not prepared_clip.intersects(work):
			continue

		work = work.intersection(aoi.clip_geom)
		if work.is_empty:
			continue

		for poly in extract_polygons(work):
			if poly.is_empty or poly.area <= 0:
				continue
			stats["polygons_after_prepare"] += 1

			prepared_poly = prep(poly)
			rings = [("outer", list(poly.exterior.coords))]
			if not cfg.only_outer_rings:
				for ring in poly.interiors:
					rings.append(("inner", list(ring.coords)))

			for ring_type, coords in rings:
				if len(coords) < 4:
					continue

				ring_line = GeoLineString(coords)
				chain = ring_chain_distances(coords)

				for i in range(1, len(coords)):
					a = coords[i - 1]
					b = coords[i]
					seg_len = segment_length(a, b)
					if seg_len < max(0.5, cfg.step_m * 0.6):
						stats["short_segments_skipped"] += 1
						continue

					local_tx, local_ty = unit_vector(a, b)
					start = max(cfg.vertex_margin_m, cfg.step_m * 0.5)
					if seg_len <= start * 2:
						continue

					s = start
					while s < seg_len - start:
						stats["stations_attempted"] += 1
						px = a[0] + local_tx * s
						py = a[1] + local_ty * s
						chain_m = chain[i - 1] + s

						tx, ty = tangent_at_distance(ring_line, chain_m, cfg.tangent_window_m, (local_tx, local_ty))
						normal = choose_inward_normal(prepared_poly, px, py, tx, ty, cfg.inner_probe_m)
						if normal is None:
							stats["failed_normal"] += 1
							s += cfg.step_m
							continue

						nx, ny = normal
						cx = px + nx * cfg.inner_probe_m
						cy = py + ny * cfg.inner_probe_m
						center = Point(cx, cy)
						cross = GeoLineString([
							(cx - nx * cfg.span_m, cy - ny * cfg.span_m),
							(cx + nx * cfg.span_m, cy + ny * cfg.span_m),
						])

						try:
							inter = poly.intersection(cross)
						except Exception:
							stats["intersection_errors"] += 1
							s += cfg.step_m
							continue

						best = best_cross_segment(inter, center, cfg.min_width_m)
						if best is None:
							stats["empty_intersections"] += 1
							s += cfg.step_m
							continue

						mid = best.interpolate(0.5, normalized=True)
						out_line = build_output_line(mid.x, mid.y, tx, ty, cfg.mark_half_m)
						out_features.append({
							"type": "Feature",
							"geometry": mapping(out_line),
							"properties": {
								"width_m": round(best.length, 2),
								"step_m": cfg.step_m,
								"ring_type": ring_type,
								"station_chain_m": round(chain_m, 2),
								"method": "smoothed_ring_cross_v1",
							},
						})
						stats["stations_written"] += 1
						s += cfg.step_m

	return out_features, stats


def write_geojson(path: Path, features: list[dict[str, Any]]) -> None:
	fc = {
		"type": "FeatureCollection",
		"features": features,
	}
	path.write_text(json.dumps(fc, ensure_ascii=False), encoding="utf-8")


def max_rss_kb() -> int:
	return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def main() -> None:
	args = parse_args()
	config = load_config(args.config)
	out_dir = Path(args.output)
	out_dir.mkdir(parents=True, exist_ok=True)

	steps = [float(x) for x in config.get("steps_m", [])]
	if not steps:
		raise ValueError("steps_m fehlt oder ist leer")

	session = requests.Session()
	resolved_aois = resolve_aois(session, config)

	summary: dict[str, Any] = {
		"created_at_unix": int(time.time()),
		"epsg_metric": EPSG_METRIC,
		"config": config,
		"runs": [],
	}

	for aoi in resolved_aois:
		fetch_t0 = time.perf_counter()
		sidewalk_geoms = fetch_sidewalk_geometries(session, aoi)
		fetch_seconds = round(time.perf_counter() - fetch_t0, 3)

		for step_m in steps:
			cfg = MeasureConfig(
				step_m=step_m,
				span_m=float(config.get("span_m", 30.0)),
				mark_half_m=float(config.get("mark_half_m", 1.2)),
				vertex_margin_m=float(config.get("vertex_margin_m", 0.5)),
				tangent_window_m=float(config.get("tangent_window_m", 6.0)),
				inner_probe_m=float(config.get("inner_probe_m", 0.08)),
				min_width_m=float(config.get("min_width_m", 0.2)),
				only_outer_rings=bool(config.get("only_outer_rings", True)),
				debug_feature_limit=int(config.get("debug_feature_limit", 3000)),
			)

			measure_t0 = time.perf_counter()
			features, stats = measure_sidewalks(sidewalk_geoms, aoi, cfg)
			measure_seconds = round(time.perf_counter() - measure_t0, 3)

			debug_features = features[: cfg.debug_feature_limit]
			stem = f"{aoi.name}_step-{str(step_m).replace('.', '_')}"
			write_geojson(out_dir / f"{stem}.geojson", debug_features)

			run_summary = {
				"aoi": aoi.name,
				"bbox": [round(v, 3) for v in aoi.bbox],
				"step_m": step_m,
				"fetch_seconds": fetch_seconds,
				"measure_seconds": measure_seconds,
				"debug_features_written": len(debug_features),
				"total_features_generated": len(features),
				"max_rss_kb": max_rss_kb(),
				"stats": stats,
			}
			summary["runs"].append(run_summary)
			print(json.dumps(run_summary, ensure_ascii=False))

	(out_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
	main()
