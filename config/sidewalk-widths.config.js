const TESTDATA_BASE_URL = `https://${process.env.EASYNAME_HTTP_HOST || 'tiles.radlobby.at'}/Gehsteigbreiten/Testdaten`;

const WORK_DIR = 'work/sidewalk-widths';

const WIDTH_CLASSES = [
	{ id: 0, label: '< 1.2 m', min: null, max: 1.2 },
	{ id: 1, label: '>= 1.2 m - < 1.5 m', min: 1.2, max: 1.5 },
	{ id: 2, label: '>= 1.5 m - < 2.0 m', min: 1.5, max: 2.0 },
	{ id: 3, label: '>= 2.0 m - < 3.0 m', min: 2.0, max: 3.0 },
	{ id: 4, label: '>= 3.0 m', min: 3.0, max: null }
];

export default {
	method: 'fmzk_guided_sis_multiband_cross_v1',

	crs: {
		input: 'EPSG:4326',
		metric: 'EPSG:31256',
		pmtiles: 'EPSG:4326'
	},

	paths: {
		workDir: WORK_DIR,
		sisIndex: `${WORK_DIR}/sis-index.sqlite`,
		metricOutput: `${WORK_DIR}/sidewalk_widths_metric.geojsonseq`,
		wgs84Output: `${WORK_DIR}/sidewalk_widths_4326.geojsonseq`,
		pmtilesOutput: 'dist/sidewalk_widths.pmtiles',
		summary: 'dist/sidewalk_widths_summary.json'
	},

	sources: {
		fmzk: {
			name: 'FMZK Gehsteig Testausschnitt',
			url: `${TESTDATA_BASE_URL}/FMZKVERKEHR2OGD.json`,
			crs: 'EPSG:4326',
			filter: {
				property: 'LAYER',
				values: ['Gehsteig']
			}
		},

		sis: {
			name: 'SIS Belagsflächen Testausschnitt',
			url: `${TESTDATA_BASE_URL}/SISBELAGOGD.json`,
			crs: 'EPSG:4326',
			typeProperty: 'TYPE',
			allowedTypes: ['GG', 'EE', 'HH']
		}
	},

	measurement: {
		stationStepM: 1.0,
		crossSectionHalfLengthM: 30.0,
		innerProbeM: 0.05,
		minWidthM: 0.2,

		fmzkSisSearchBufferM: 35.0,

		bandTrackingMaxDistanceM: 2.0,
		maxStationGapM: 2.5,

		outputSimplifyToleranceM: 0.2,
		outputValidationSampleStepM: 0.5,

		directionSmoothingWindowStations: 2,
		directionSpreadWindowStations: 2,
		directionOutlierDeg: 45,
		directionHardOutlierDeg: 70,
		directionStableSpreadDeg: 25
	},

	widthClasses: WIDTH_CLASSES,

	pmtiles: {
		layerName: 'sidewalk_widths',
		minZoom: 12,
		maxZoom: 17,
		baseZoom: 14
	}
};
