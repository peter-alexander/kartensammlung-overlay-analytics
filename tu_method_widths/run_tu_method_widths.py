#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path
from urllib.request import urlretrieve

import geopandas as gpd
from shapely import make_valid
from shapely.geometry import Point
from shapely.ops import unary_union


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', default='config/tu-method-widths.json')