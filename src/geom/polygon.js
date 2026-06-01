import { ensureClosedRing } from './basic.js';

export function extractPolygons(geometry) {
	if (!geometry) return [];
	if (geometry.type === 'Polygon') return [geometry.coordinates];
	if (geometry.type === 'MultiPolygon') return geometry.coordinates;
	if (geometry.type === 'GeometryCollection') {
		return (geometry.geometries || []).flatMap(extractPolygons);
	}
	return [];
}

export function pointInRing(point, ring) {
	const x = point[0];
	const y = point[1];
	let inside = false;
	const r = ensureClosedRing(ring);

	for (let i = 0, j = r.length - 1; i < r.length; j = i++) {
		const xi = r[i][0];
		const yi = r[i][1];
		const xj = r[j][0];
		const yj = r[j][1];
		const crosses = ((yi > y) !== (yj > y)) &&
			(x < ((xj - xi) * (y - yi)) / ((yj - yi) || 1e-30) + xi);
		if (crosses) inside = !inside;
	}

	return inside;
}

export function pointInPolygonCoords(point, polygonCoords) {
	if (!polygonCoords?.length) return false;
	if (!pointInRing(point, polygonCoords[0])) return false;
	for (let i = 1; i < polygonCoords.length; i++) {
		if (pointInRing(point, polygonCoords[i])) return false;
	}
	return true;
}
