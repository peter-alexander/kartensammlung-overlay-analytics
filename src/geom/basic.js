export function dist(a, b) {
	return Math.hypot(b[0] - a[0], b[1] - a[1]);
}

export function lerp(a, b, t) {
	return [
		a[0] + (b[0] - a[0]) * t,
		a[1] + (b[1] - a[1]) * t
	];
}

export function normalize(dx, dy) {
	const len = Math.hypot(dx, dy);
	if (len <= 1e-12) return [0, 0];
	return [dx / len, dy / len];
}

export function ringSignedArea(ring) {
	let sum = 0;
	for (let i = 0; i < ring.length - 1; i++) {
		const a = ring[i];
		const b = ring[i + 1];
		sum += (a[0] * b[1]) - (b[0] * a[1]);
	}
	return sum * 0.5;
}

export function ensureClosedRing(ring) {
	if (!ring.length) return ring;
	const first = ring[0];
	const last = ring[ring.length - 1];
	if (first[0] === last[0] && first[1] === last[1]) return ring;
	return [...ring, first];
}
