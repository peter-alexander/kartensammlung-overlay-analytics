const TESTDATA_BASE_URL = `https://${process.env.EASYNAME_HTTP_HOST || 'tiles.radlobby.at'}/Gehsteigbreiten/Testdaten`;

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
			url: `${TESTDATA_BASE_URL}/FMZKVERKEHR2OGD.json`,
			filter: {
				property: 'LAYER',
	