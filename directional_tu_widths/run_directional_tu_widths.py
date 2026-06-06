#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path
from urllib.request import urlretrieve

import geopandas as gpd
from shapely import make_valid
from shapely.geometry import LineString


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', default='config/directional-tu-widths.json')
	args = parser.parse_args()

	config = load_json(args.config)
	paths = config['paths']
	work_dir = Path(paths['work_dir'])
	Path(paths['output_dir']).mkdir(parents