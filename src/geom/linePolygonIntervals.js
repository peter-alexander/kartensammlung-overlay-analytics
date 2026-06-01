import { dist, lerp } from './basic.js';
import { extractPolygons, pointInPolygonCoords } from './polygon.js';

const EPS = 1e-9;

function segmentLineT(lineA, lineB, segA, segB) {
	const x1 = lineA[0];
	const y1 = lineA[1];
	const x2 = lineB[0];
	const y2 = lineB[1];
	const x3 = segA[0];
	const y3 = segA[1];
	const x4 = segB[0];
	const y4 = segB[1];

	const den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	if (Math.abs(den) < EPS) return null;

	const t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
	const u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den;

	if (t < -EPS || t > 1 + EPS || u < -EPS || u > 1 + EPS) return null;
	return Math.max(0, Math.min(1, t));
}

function collectRingIntersections(lineA, lineB, ring, out) {
	for (let i = 1; i < ring.length; i++) {
		const t = segmentLineT(lineA, lineB, ring[i - 1], ring[i]);
		if (t !== null) out.push(t);
	}
}

export function lineGeometryInsideIntervals(lineA, lineB, geometry) {
	const lineLen = dist(lineA, lineB);
	if (lineLen <= EPS) return [];

	const intervals = [];
	for (const polygon of extractPolygons(geometry)) {
		const ts = [0, 1];
		for (const ring of polygon) collectRingIntersections(lineA, lineB, ring, ts);
		const sorted = [...new Set(ts.map((t) => Number(t.toFixed(12))))].sort((a, b) => a - b);

		for (let i = 1; i < sorted.length; i++) {
			const t0 = sorted[i - 1];
			const t1 = sorted[i];
			if (t1 - t0 <= EPS) continue;
			const mid = lerp(lineA, lineB, (t0 + t1) * 0.5);
			if (pointInPolygonCoords(mid, polygon)) {
				intervals.push({ t0, t1 });
			}
		}
	}

	return mergeIntervals(intervals);
}

export function mergeIntervals(intervals) {
	if (!intervals.length) return [];
	const sorted = intervals
		.map((i) => ({ t0: Math.max(0, i.t0), t1: Math.min(1, i.t1) }))
		.filter((i) => i.t1 - i.t0 > EPS)
		.sort((a, b) => a.t0 - b.t0);

	const out = [];
	for (const item of sorted) {
		const last = out[out.length - 1];
		if (last && item.t0 <= last.t1 + EPS) {
			last.t1 = Math.max(last.t1, item.t1);
		} else {
			out.push({ ...item });
		}
	}
	return out;
}
