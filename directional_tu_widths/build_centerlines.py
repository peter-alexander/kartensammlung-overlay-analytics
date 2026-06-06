#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path

import pygeoops

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
	sys.path.insert(0, str(REPO_ROOT))

from centerline_widths.run_centerline_widths import (  # noqa: E402
	download,
	explode_lines,
	load_prepare_source,
	mapping_line,
	postprocess_centerlines,
	prepare_centerline_geom,
	write_geojsonseq,
)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', default='config/directional-tu-widths.json')
	args = parser.parse_args()

	config = load_json(args.config)
	paths = config['paths']
	work_dir = Path(paths['work_dir'])
	output_dir = Path(paths['output_dir'])
	work_dir.mkdir(parents=True, exist_ok=True)
	output_dir.mkdir(parents=True, exist_ok=True)

	build_cfg = config['centerline_build']
	centerline_config = {
		'crs': config['crs'],
		'centerline': build_cfg['centerline']
	}
	fmzk_source = config['sources']['fmzk']
	centerline_output = Path(config['sources']['centerlines']['url'])
	if str(centerline_output).startswith(('http://', 'https://')):
		raise RuntimeError('sources.centerlines.url must be a local output path for generated centerlines')

	fmzk_raw = work_dir / 'fmzk_raw.geojson'
	download(fmzk_source['url'], fmzk_raw)
	fmzk = load_prepare_source(fmzk_raw, centerline_config, fmzk_source)

	summary = {
		'method': 'directional_centerline_build',
		'fmzk_features': int(len(fmzk)),
		'centerline_ok': 0,
		'centerline_failed': 0,
		'centerline_parts': 0,
		'centerline_merged': 0,
		'centerline_post_simplified': 0,
		'centerline_straightened': 0,
		'centerline_straight_simplified': 0,
		'errors': {}
	}
	features = []

	for idx, row in fmzk.iterrows():
		geom = row.geometry
		props = row.drop(labels=['geometry']).to_dict()
		fmzk_id = props.get('FMZK_ID') or props.get('OBJECTID') or props.get('ID') or idx

		try:
			geom_clean = prepare_centerline_geom(geom, centerline_config)
			if geom_clean is None or geom_clean.is_empty:
				continue
			centerline = pygeoops.centerline(
				geom_clean,
				densify_distance=centerline_config['centerline']['densify_distance'],
				min_branch_length=centerline_config['centerline']['min_branch_length'],
				simplifytolerance=centerline_config['centerline']['simplifytolerance'],
				extend=centerline_config['centerline']['extend']
			)
			lines = explode_lines(centerline)
			lines = postprocess_centerlines(lines, centerline_config, summary)
			lines = [line for line in lines if line.length >= centerline_config['centerline']['min_line_length_m']]
			summary['centerline_ok'] += 1
		except Exception as exc:
			summary['centerline_failed'] += 1
			key = type(exc).__name__
			summary['errors'][key] = summary['errors'].get(key, 0) + 1
			continue

		summary['centerline_parts'] += len(lines)
		for part_index, line in enumerate(lines):
			features.append({
				'type': 'Feature',
				'geometry': mapping_line(line),
				'properties': {
					'fmzk_id': fmzk_id,
					'part_index': part_index,
					'length_m': round(line.length, 2),
					'source': 'pygeoops.centerline',
					'post_simplify_tolerance_m': centerline_config['centerline'].get('post_simplify_tolerance_m', 0)
				}
			})

	write_geojsonseq(centerline_output, features)
	Path(build_cfg['summary']).parent.mkdir(parents=True, exist_ok=True)
	Path(build_cfg['summary']).write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
	print(json.dumps(summary, ensure_ascii=False))


def load_json(path):
	return json.loads(Path(path).read_text(encoding='utf-8'))


if __name__ == '__main__':
	main()
