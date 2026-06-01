import { dist } from './basic.js';

export function simplifyLine(points, tolerance) {
	if (points.length <= 2 || tolerance <= 0) return removeDuplicatePoints(points);
	const deduped = removeDuplicatePoints(points);
	if (deduped.length <= 2) return deduped;
	return douglasPeucker(deduped, tolerance);
}

function removeDuplicatePoints(points) {
	const out = [];
	for (const p of points) {
		const last = out[out.length - 1];
		if (!last || dist(last, p) > 1e-9) out.push(p);
	}
	return out;
}

function douglasPeucker(points, tolerance) {
	let maxDist = 0;
	let index = 0;
	const start = points[0];
	const end = points[points.length - 1];

	for (let i = 1; i < points.length - 1; i++) {
		const d = perpendicularDistance(points[i], start, end);
		if (d > maxDist) {
			maxDist = d;
			index = i;
		}
	}

	if (maxDist > tolerance) {
		const left = douglasPeucker(points.slice(0, index + 1), tolerance);
		const right = douglasPeucker(points.slice(index), tolerance);
		return left.slice(0, -1).concat(right);
	}

	return [start, end];
}

function perpendicularDistance(p, a, b) {
	const dx = b[0] - a[0];
	const dy = b[1] - a[1];
	const len2 = dx * dx + dy * dy;
	if (len2 <= 1e-12) return dist(p, a);
	let t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2;
	t = Math.max(0, Math.min(1, t));
	const q = [a[0] + dx * t, a[1] + dy * t];
	return dist(p, q);
}
