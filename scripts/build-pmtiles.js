#!/usr/bin/env node
import fs from 'node:fs';
import path from 'node:path';
import { spawn } from 'node:child_process';
import { loadConfig } from '../src/config/loadConfig.js';

export async function buildPmtiles(config) {
	await fs.promises.mkdir(path.dirname(config.paths.pmtilesOutput), { recursive: true });

	if (!fs.existsSync(config.paths.wgs84Output)) {
		console.warn(`WGS84-Output fehlt, verwende metrischen Output direkt: ${config.paths.metricOutput}`);
	}

	const input = fs.existsSync(config.paths.wgs84Output)
		? config.paths.wgs84Output
		: config.paths.metricOutput;

	const args = [
		'-o', config.paths.pmtilesOutput,
		'--force',
		'--layer', config.pmtiles.layerName,
		'--minimum-zoom', String(config.pmtiles.minZoom),
		'--maximum-zoom', String(config.pmtiles.maxZoom),
		'--base-zoom', String(config.pmtiles.baseZoom),
		'--drop-densest-as-needed',
		input
	];

	await runCommand('tippecanoe', args);
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
