import fs from 'node:fs';
import path from 'node:path';
import { downloadToFile, ensureDir } from '../io/download.js';
import { iterateGeoJsonFeatures, writeGeoJsonSeqFeature } from '../io/geojsonFeatures.js';
import { buildSisIndex, openSisIndex } from '../io/sisIndex.js';
import { measureFmzkFeature } from './measure.js';

export async function runSidewalkWidthsPipeline(config) {
	await ensureDir(config.paths.workDir);
	await ensureDir(path.dirname(config.paths.metricOutput));
	await ensureDir(path.dirname(config.paths.pmtilesOutput));

	const fmzkPath = path.join(config.paths.workDir, 'fmzk.geojson');
	const sisPath = path.join(config.paths.workDir, 'sis.geojson');

	console.log(`Download FMZK: ${config.sources.fmzk.url}`);
	await downloadToFile(config.sources.fmzk.url, fmzkPath);

	console.log(`Download SIS: ${config.sources.sis.url}`);
	await downloadToFile(config.sources.sis.url, sisPath);

	console.log('Baue SIS SQLite/RTree Index');
	const sisIndexStats = await buildSisIndex({
		sourcePath: sisPath,
		dbPath: config.paths.sisIndex,
		typeProperty: config.sources.sis.typeProperty,
		allowedTypes: config.sources.sis.allowedTypes
	});
	console.log(JSON.stringify({ sisIndexStats }));

	const sisDb = openSisIndex(config.paths.sisIndex);
	const output = fs.createWriteStream(config.paths.metricOutput, { encoding: 'utf8' });

	const summary = {
		method: config.method,
		created_at: new Date().toISOString(),
		crs_metric: config.crs.metric,
		width_classes: config.widthClasses,
		sources: {
			fmzk: config.sources.fmzk.name,
			sis: config.sources.sis.name
		},
		sis_index: sisIndexStats,
		fmzk_features: 0,
		fmzk_features_filtered_out: 0,
		fmzk_features_processed: 0,
		output_features: 0,
		fmzk_diagnostics: {
			geometryTypeCounts: {},
			samplePropertyKeys: []
		},
		stats: {
			fmzkPolygons: 0,
			sisCandidates: 0,
			stations: 0,
			crossSections: 0,
			measurements: 0
		}
	};

	try {
		for await (const feature of iterateGeoJsonFeatures(fmzkPath)) {
			summary.fmzk_features++;
			collectFmzkDiagnostics(summary.fmzk_diagnostics, feature);

			if (!passesFeatureFilter(feature, config.sources.fmzk.filter)) {
				summary.fmzk_features_filtered_out++;
				continue;
			}

			summary.fmzk_features_processed++;
			const result = measureFmzkFeature({ feature, sisDb, config });

			for (const outFeature of result.features) {
				writeGeoJsonSeqFeature(output, outFeature);
				summary.output_features++;
			}

			mergeStats(summary.stats, result.stats);

			if (summary.fmzk_features % 100 === 0) {
				console.log(JSON.stringify({
					fmzk_features: summary.fmzk_features,
					fmzk_features_processed: summary.fmzk_features_processed,
					output_features: summary.output_features,
					stations: summary.stats.stations,
					measurements: summary.stats.measurements
				}));
			}
		}
	} finally {
		await new Promise((resolve) => output.end(resolve));
		sisDb.close();
	}

	await fs.promises.writeFile(config.paths.summary, JSON.stringify(summary, null, 2), 'utf8');
	console.log(JSON.stringify({ done: true, summary }));
	return summary;
}

function passesFeatureFilter(feature, filter) {
	if (!filter) return true;
	const value = feature.properties?.[filter.property];
	const allowed = new Set(filter.values || []);
	return allowed.has(value);
}

function collectFmzkDiagnostics(target, feature) {
	const geomType = feature.geometry?.type || 'null';
	target.geometryTypeCounts[geomType] = (target.geometryTypeCounts[geomType] || 0) + 1;

	const properties = feature.properties || {};
	const keys = new Set(target.samplePropertyKeys);
	for (const key of Object.keys(properties).slice(0, 40)) keys.add(key);
	target.samplePropertyKeys = [...keys];
}

function mergeStats(target, source) {
	for (const [key, value] of Object.entries(source)) {
		if (key === 'fmzkFeatures') continue;
		if (typeof value === 'number') target[key] = (target[key] || 0) + value;
	}
}
