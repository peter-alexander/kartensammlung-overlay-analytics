import fs from 'node:fs';
import readline from 'node:readline';
import { chain } from 'stream-chain';
import { parser } from 'stream-json';
import { pick } from 'stream-json/filters/Pick.js';
import { streamArray } from 'stream-json/streamers/StreamArray.js';

export async function* iterateGeoJsonFeatures(filePath) {
	if (/(\.geojsonseq|\.ndjson|\.jsonl)$/i.test(filePath)) {
		yield* iterateGeoJsonSeq(filePath);
		return;
	}

	const pipeline = chain([
		fs.createReadStream(filePath),
		parser(),
		pick({ filter: 'features' }),
		streamArray()
	]);

	for await (const item of pipeline) {
		yield item.value;
	}
}

async function* iterateGeoJsonSeq(filePath) {
	const rl = readline.createInterface({
		input: fs.createReadStream(filePath),
		crlfDelay: Infinity
	});

	for await (const line of rl) {
		const trimmed = line.trim();
		if (!trimmed || trimmed === '\u001e') continue;
		yield JSON.parse(trimmed.replace(/^\u001e/, ''));
	}
}

export function writeGeoJsonSeqFeature(stream, feature) {
	stream.write(JSON.stringify(feature));
	stream.write('\n');
}
