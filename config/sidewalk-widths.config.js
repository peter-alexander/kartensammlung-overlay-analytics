export default {
	method: 'fmzk_guided_sis_multiband_cross_v1',
	crs: {
		metric: 'EPSG:31256',
		pmtiles: 'EPSG:4326'
	},
	paths: {
		workDir: 'work/sidewalk-widths',
		sisIndex: 'work/sidewalk-widths/sis-index.sqlite',
		metricOutput: 'work/sidewalk-widths/sidewalk_widths_metric.geojsonseq',
		wgs84Output: 'work/sidewalk-widths/sidewalk_widths_4326.geojsonseq',
		pmtilesOutput: 'dist/sidewalk_widths.pmtiles',
		summary: 'dist/sidewalk_widths_summary.json'
	},
	sources: {
		fmzk: {
			name: 'FMZK Gehsteig Testausschnitt',
			kind: 'geojson',
			url: 'https://fahrrad.lima-city.de/Gehsteigbreiten-Testdaten/FMZKVERKEHR2OGD.json',
			filter: {
				property: 'LAYER',
				values: ['Gehsteig']
			}
		},
		sis: {
			name: 'SIS Belagsflächen Testausschnitt',
			kind: 'geojson',
			url: 'https://fahrrad.lima-city.de/Gehsteigbreiten-Testdaten/SISBELAGOGD.json',
			typeProperty: 'TYPE',
			allowedTypes: ['GG', 'EE', 'HH']
		}
	},
	measurement: {
		stepM: 1.0,
		crossSectionHalfLengthM: 30.0,
		innerProbeM: 0.08,
		stationMarginM: 0.5,
		minWidthM: 0.2,
		fmzkSisSearchBufferM: 1.0,
		directionWindowStations: 3,
		directionOutlierDeg: 45,
		directionHardOutlierDeg: 70,
		directionStableSpreadDeg: 25,
		simplifyToleranceM: 0.2,
		maxStationGapM: 2.5,
		outputValidationSampleStepM: 0.5,
		bandTrackingMaxDistanceM: 2.0
	},
	widthClasses: [
		{ id: 0, label: '< 1.2 m', min: null, max: 1.2 },
		{ id: 1, label: '1.2 m - 1.5 m', min: 1.2, max: 1.5 },
		{ id: 2, label: '1.5 m - 2.0 m', min: 1.5, max: 2.0 },
		{ id: 3, label: '2.0 m - 3.0 m', min: 2.0, max: 3.0 },
		{ id: 4, label: '> 3.0 m', min: 3.0, max: null }
	],
	pmtiles: {
		layerName: 'sidewalk_widths',
		minZoom: 12,
		maxZoom: 17,
		baseZoom: 14
	}
};
