#!/usr/bin/env node
import fs from 'node:fs';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { loadConfig } from '../src/config/loadConfig.js';

export async function buildPmtiles(config) {
	await fs.promises.mkdir(path.dirname(config.paths.pmtilesOutput), { recursive: true });
	await fs.promises.mkdir(path.dirname(config.paths.wgs84Output), { recursive: true });

	if (!fs.existsSync(config.paths.metricOutput)) {
		console.warn(`Metrischer Output fehlt, PMTiles-Build wird übersprungen: ${config.paths.metricOutput}`);
		return;
	}

	const stat = await fs.promises.stat(config.paths.metricOutput);
	if (stat.size === 0) {
		console.warn(`Metrischer Output ist leer, PMTiles-Build wird übersprungen: ${config.paths.metricOutput}`);
		return;
	}

	await transformMetricOutputToWgs84(config);

	const args = [
		'-o', config.paths.pmtilesOutput,
		'--force',
		'--layer', config.pmtiles.layerName,
		'--minimum-zoom', String(config.pmtiles.minZoom),
		'--maximum-zoom', String(config.pmtiles.maxZoom),
		'--base-zoom', String(config.pmtiles.baseZoom),
		'--drop-densest-as-needed',
		config.paths.wgs84Output
	];

	await runCommand('tippecanoe', args);
}

async function transformMetricOutputToWgs84(config) {
	await fs.promises.rm(config.paths.wgs84Output, { force: true });

	await runCommand('ogr2ogr', [
		'-f', 'GeoJSONSeq',
		'-s_srs', config.crs.metric,
		'-t_srs', config.crs.pmtiles,
		config.paths.wgs84Output,
		config.paths.metricOutput
	]);
}

function runCommand(command, args) {
	return new Promise((resolve, reject) => {
		const child = spawn(command, args, { stdio: 'inherit' });
		child.on('error', reject);
		child.on('close', (code) => {
			if (code === 0) resolve();
			else reject(new Error(`${command} beendet mit Exit-Code ${code}`));
		});
	});
}

if (import.meta.url === `file://${process.argv[1]}`) {
	const idx = process.argv.indexOf('--config');
	const configPath = idx === -1 ? 'config/sidewalk-widths.config.js' : process.argv[idx + 1];
	const config = await loadConfig(configPath);
	await buildPmtiles(config);
}
