import fs from 'node:fs';
import path from 'node:path';
import Database from 'better-sqlite3';
import { iterateGeoJsonFeatures } from './geojsonFeatures.js';
import { bboxOfGeometry } from '../geom/bbox.js';

export function openSisIndex(dbPath) {
	return new Database(dbPath);
}

export async function buildSisIndex({ sourcePath, dbPath, typeProperty, allowedTypes }) {
	await fs.promises.mkdir(path.dirname(dbPath), { recursive: true });
	await fs.promises.rm(dbPath, { force: true });

	const db = new Database(dbPath);
	db.pragma('journal_mode = WAL');
	db.pragma('synchronous = NORMAL');

	db.exec(`
		CREATE TABLE sis_features (
			id INTEGER PRIMARY KEY,
			type TEXT,
			min_x REAL NOT NULL,
			min_y REAL NOT NULL,
			max_x REAL NOT NULL,
			max_y REAL NOT NULL,
			geometry_json TEXT NOT NULL,
			properties_json TEXT NOT NULL
		);
		CREATE VIRTUAL TABLE sis_rtree USING rtree(
			id,
			min_x,
			max_x,
			min_y,
			max_y
		);
	`);

	const insertFeature = db.prepare(`
		INSERT INTO sis_features (type, min_x, min_y, max_x, max_y, geometry_json, properties_json)
		VALUES (?, ?, ?, ?, ?, ?, ?)
	`);
	const insertRtree = db.prepare(`
		INSERT INTO sis_rtree (id, min_x, max_x, min_y, max_y)
		VALUES (?, ?, ?, ?, ?)
	`);

	let read = 0;
	let indexed = 0;
	const allowed = new Set(allowedTypes || []);

	const batch = db.transaction((items) => {
		for (const feature of items) {
			const properties = feature.properties || {};
			const type = String(properties[typeProperty] ?? '').trim();
			if (allowed.size && !allowed.has(type)) continue;
			if (!feature.geometry) continue;

			const bbox = bboxOfGeometry(feature.geometry);
			if (!bbox) continue;

			const result = insertFeature.run(
				type,
				bbox.minX,
				bbox.minY,
				bbox.maxX,
				bbox.maxY,
				JSON.stringify(feature.geometry),
				JSON.stringify(properties)
			);
			insertRtree.run(result.lastInsertRowid, bbox.minX, bbox.maxX, bbox.minY, bbox.maxY);
			indexed++;
		}
	});

	let buffer = [];
	for await (const feature of iterateGeoJsonFeatures(sourcePath)) {
		read++;
		buffer.push(feature);
		if (buffer.length >= 1000) {
			batch(buffer);
			buffer = [];
		}
	}
	if (buffer.length) batch(buffer);

	db.close();
	return { read, indexed };
}

export function querySisIndex(db, bbox) {
	const stmt = db.prepare(`
		SELECT f.id, f.type, f.geometry_json, f.properties_json
		FROM sis_rtree r
		JOIN sis_features f ON f.id = r.id
		WHERE r.min_x <= @maxX
			AND r.max_x >= @minX
			AND r.min_y <= @maxY
			AND r.max_y >= @minY
	`);

	return stmt.all(bbox).map((row) => ({
		id: row.id,
		type: row.type,
		geometry: JSON.parse(row.geometry_json),
		properties: JSON.parse(row.properties_json)
	}));
}
