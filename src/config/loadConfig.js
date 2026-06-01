import { pathToFileURL } from 'node:url';
import path from 'node:path';

export async function loadConfig(configPath) {
	const abs = path.resolve(configPath);
	const mod = await import(pathToFileURL(abs).href + `?t=${Date.now()}`);
	const config = mod.default;

	if (!config || typeof config !== 'object') {
		throw new Error(`Config ${configPath} exportiert kein default-Objekt.`);
	}

	return config;
}
