#!/usr/bin/env node
import { loadConfig } from '../src/config/loadConfig.js';
import { runSidewalkWidthsPipeline } from '../src/sidewalk/runPipeline.js';
import { buildPmtiles } from './build-pmtiles.js';

function argValue(name, fallback = null) {
	const idx = process.argv.indexOf(name);
	if (idx === -1 || idx + 1 >= process.argv.length) return fallback;
	return process.argv[idx + 1];
}

const configPath = argValue('--config', 'config/sidewalk-widths.config.js');
const skipPmtiles = process.argv.includes('--skip-pmtiles');

const config = await loadConfig(configPath);
await runSidewalkWidthsPipeline(config);

if (!skipPmtiles) {
	await buildPmtiles(config);
}
