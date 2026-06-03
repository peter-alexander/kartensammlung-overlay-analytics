#!/usr/bin/env node
import { loadConfig } from '../src/config/loadConfig.js';
import { downloadToFile } from '../src/io/download.js';
import { buildSisIndex } from '../src/io/sisIndex.js';

function argValue(name, fallback = null) {
	const idx = process.argv.indexOf(name);
	if (idx === -1 || idx + 1 >= process.argv.length) return fallback;
	return process.argv[idx + 1];
}

const configPath = argValue('--config', 'config/sidewalk-widths.config.js');
const config = await loadConfig(configPath);
const sourcePath = `${config.paths.workDir}/sis.geojson`;

await downloadToFile(config.sources.sis.url, sourcePath);
const stats = await buildSisIndex({
	sourcePath,
	dbPath: config.paths.sisIndex,
	typeProperty: config.sources.sis.typeProperty,
	allowedTypes: config.sources.sis.allowedTypes
});

console.log(JSON.stringify(stats, null, 2));
