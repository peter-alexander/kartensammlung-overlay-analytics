import fs from 'node:fs';
import readline from 'node:readline';
import chainPkg from 'stream-chain';
import streamJsonPkg from 'stream-json';
import pickPkg from 'stream-json/filters/Pick.js';
import streamArrayPkg from 'stream-json/streamers/StreamArray.js';

const { chain } = chainPkg;
const { parser } = streamJsonPkg;
const { pick } = pickPkg;
const { streamArray } = streamArrayPkg;

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
