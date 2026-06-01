import fs from 'node:fs';
import path from 'node:path';
import { pipeline } from 'node:stream/promises';

export async function ensureDir(dir) {
	await fs.promises.mkdir(dir, { recursive: true });
}

export async function downloadToFile(url, targetPath) {
	await ensureDir(path.dirname(targetPath));

	if (url.startsWith('file://')) {
		const sourcePath = new URL(url);
		await fs.promises.copyFile(sourcePath, targetPath);
		return targetPath;
	}

	if (!/^https?:\/\//i.test(url)) {
		await fs.promises.copyFile(url, targetPath);
		return targetPath;
	}

	const res = await fetch(url);
	if (!res.ok) {
		throw new Error(`Download fehlgeschlagen ${res.status} ${res.statusText}: ${url}`);
	}
	if (!res.body) {
		throw new Error(`Download ohne Body: ${url}`);
	}

	await pipeline(res.body, fs.createWriteStream(targetPath));
	return targetPath;
}
