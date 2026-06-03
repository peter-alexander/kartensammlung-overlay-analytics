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
	let skippedByType = 0;
	let skippedNoGeometry = 0;
	let skippedNoBbox = 0;
	const allowed = new Set(allowedTypes || []);
	const typeCounts = new Map();
	const geometryTypeCounts = new Map();
	const samplePropertyKeys = new Set();

	const prepared = [];
	for await (const feature of iterateGeoJsonFeatures(sourcePath)) {
		read++;
		const properties = feature.properties || {};
		const type = String(properties[typeProperty] ?? '').trim();
		const geomType = feature.geometry?.type || 'null';
		typeCounts.set(type || '(empty)', (typeCounts.get(type || '(empty)') || 0) + 1);
		geometryTypeCounts.set(geomType, (geometryTypeCounts.get(geomType) || 0) + 1);
		for (const key of Object.keys(properties).slice(0, 30)) samplePropertyKeys.add(key);
		prepared.push({ feature, type });
	}

	const hasAnyAllowedType = allowed.size === 0 || prepared.some((item) => allowed.has(item.type));
	const applyTypeFilter = allowed.size > 0 && hasAnyAllowedType;

	const batch = db.transaction((items) => {
		for (const item of items) {
			const feature = item.feature;
			const type = item.type;
			const properties = feature.properties || {};

			if (applyTypeFilter && !allowed.has(type)) {
				skippedByType++;
				continue;
			}
			if (!feature.geometry) {
				skippedNoGeometry++;
				continue;
			}

			const bbox = bboxOfGeometry(feature.geometry);
			if (!bbox) {
				skippedNoBbox++;
				continue;
			}

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

	for (let i = 0; i < prepared.length; i += 1000) {
		batch(prepared.slice(i, i + 1000));
	}

	db.close();
	return {
		read,
		indexed,
		skippedByType,
		skippedNoGeometry,
		skippedNoBbox,
		typeProperty,
		typeFilterApplied: applyTypeFilter,
		allowedTypes: [...allowed],
		typeCounts: Object.fromEntries(typeCounts),
		geometryTypeCounts: Object.fromEntries(geometryTypeCounts),
		samplePropertyKeys: [...samplePropertyKeys]
	};
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
