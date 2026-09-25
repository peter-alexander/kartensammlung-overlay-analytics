import importlib.util
import unittest
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, box


MODULE_PATH = Path(__file__).resolve().parents[1] / 'hrvn_cobblestone' / 'run_hrvn_cobblestone.py'
spec = importlib.util.spec_from_file_location('hrvn_cobblestone_runner', MODULE_PATH)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class HrvnCobblestoneTests(unittest.TestCase):
	def test_build_tiles_uses_requested_grid_size(self):
		tiles = runner.build_tiles((0, 0, 4000, 2000), 2000)
		self.assertEqual(len(tiles), 2)
		self.assertEqual(tuples(t.bounds for t in tiles), [(0.0, 0.0, 2000.0, 2000.0), (2000.0, 0.0, 4000.0, 2000.0)])

	def test_extract_filters_types_and_clips_to_tile(self):
		features = [
			feature('GM', [[[-10, 0], [10, 0], [10, 10], [-10, 10], [-10, 0]]]),
			feature('AB', [[[0, 0], [5, 0], [5, 5], [0, 5], [0, 0]]])
		]
		config = {'sources': {'sis': {'filter_property': 'TYPE', 'filter_values': ['GM', 'GO', 'KL']}}}
		geometries = runner.extract_cobblestone_geometries(features, config, box(0, 0, 20, 20))
		self.assertEqual(len(geometries), 1)
		self.assertAlmostEqual(geometries[0].area, 100.0)

	def test_geometry_property_from_xsd(self):
		xml = b'''<?xml version="1.0"?>
		<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:gml="http://www.opengis.net/gml">
			<xsd:complexType name="SISBELAGOGDType">
				<xsd:sequence>
					<xsd:element name="TYPE" type="xsd:string"/>
					<xsd:element name="SHAPE" type="gml:MultiSurfacePropertyType"/>
				</xsd:sequence>
			</xsd:complexType>
		</xsd:schema>'''
		self.assertEqual(runner.geometry_property_from_xsd(xml), 'SHAPE')

	def test_clip_hrvn_to_cobblestone_uses_unbuffered_lines(self):
		hrvn = gpd.GeoDataFrame(
			[
				{'M18_RANG_SUB': 'B', 'geometry': LineString([(-5, 5), (15, 5)])},
				{'M18_RANG_SUB': 'G', 'geometry': LineString([(-5, 20), (15, 20)])}
			],
			geometry='geometry',
			crs='EPSG:31256'
		)
		result = runner.clip_hrvn_to_cobblestone(hrvn, box(0, 0, 10, 10), 'M18_RANG_SUB')
		self.assertEqual(len(result), 1)
		self.assertEqual(result.iloc[0]['rank'], 'B')
		self.assertAlmostEqual(result.iloc[0]['length_m'], 10.0)

	def test_simplify_preserves_polygonal_geometry(self):
		geom = box(0, 0, 10, 10)
		result = runner.simplify_polygonal(geom, 0.2)
		self.assertEqual(result.geom_type, 'Polygon')
		self.assertAlmostEqual(result.area, 100.0)


def feature(surface_type, coordinates):
	return {
		'type': 'Feature',
		'properties': {'TYPE': surface_type},
		'geometry': {'type': 'Polygon', 'coordinates': coordinates}
	}


def tuples(values):
	return [tuple(value) for value in values]


if __name__ == '__main__':
	unittest.main()
