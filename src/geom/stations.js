import { dist, lerp, normalize, ringSignedArea } from './basic.js';
import { axisAngleDeg, classifyDirection, localAxisSpread, smoothAxisAngles } from '../sidewalk/directionSmoothing.js';
import { pointInPolygonCoords } from './polygon.js';

export function buildRingStations(ring, polygonCoords, cfg) {
	const stations = [];
	const isCcw = ringSignedArea(ring) > 0;

	for (let i = 1; i < ring.length; i++) {
		const a = ring[i - 1];
		const b = ring[i];
		const segLen = dist(a, b);
		if (segLen < Math.max(0.5, cfg.stepM * 0.6)) continue;

		const [tx, ty] = normalize(b[0] - a[0], b[1] - a[1]);
		const start = Math.max(cfg.stationMarginM, cfg.stepM * 0.5);
		if (segLen <= start * 2) continue;

		for (let s = start; s < segLen - start; s += cfg.stepM) {
			const point = lerp(a, b, s / segLen);
			stations.push({
				point,
				chainM: null,
				rawTx: tx,
				rawTy: ty,
				rawDeg: axisAngleDeg(tx, ty),
				segIndex: i - 1,
				segOffsetM: s,
				isCcw
			});
		}
	}

	let chain = 0;
	let cursor = 0;
	for (let i = 1; i < ring.length; i++) {
		const a = ring[i - 1];
		const b = ring[i];
		const segLen = dist(a, b);
		while (cursor < stations.length && stations[cursor].segIndex === i - 1) {
			stations[cursor].chainM = chain + stations[cursor].segOffsetM;
			cursor++;
		}
		chain += segLen;
	}

	const rawAngles = stations.map((s) => s.rawDeg);
	const usedAngles = smoothAxisAngles(rawAngles, cfg.directionWindowStations);

	for (let i = 0; i < stations.length; i++) {
		const theta = usedAngles[i] * Math.PI / 180;
		const tx = Math.cos(theta);
		const ty = Math.sin(theta);
		let normal = chooseInwardNormal(stations[i].point, tx, ty, polygonCoords, cfg.innerProbeM);

		if (!normal) {
			normal = stations[i].isCcw ? [-ty, tx] : [ty, -tx];
		}

		const spread = localAxisSpread(rawAngles, i, cfg.directionWindowStations);
		const directionQuality = classifyDirection(rawAngles[i], usedAngles[i], spread, cfg);

		stations[i].usedTx = tx;
		stations[i].usedTy = ty;
		stations[i].normalX = normal[0];
		stations[i].normalY = normal[1];
		stations[i].usedDeg = usedAngles[i];
		stations[i].directionSpreadDeg = spread;
		stations[i].directionQuality = directionQuality.quality;
		stations[i].directionDeltaDeg = directionQuality.deltaDeg;
	}

	return stations;
}

function chooseInwardNormal(point, tx, ty, polygonCoords, probeM) {
	const left = [-ty, tx];
	const right = [ty, -tx];

	for (const scale of [1, 2, 4]) {
		const probe = probeM * scale;
		const lp = [point[0] + left[0] * probe, point[1] + left[1] * probe];
		const rp = [point[0] + right[0] * probe, point[1] + right[1] * probe];
		const leftInside = pointInPolygonCoords(lp, polygonCoords);
		const rightInside = pointInPolygonCoords(rp, polygonCoords);

		if (leftInside && !rightInside) return left;
		if (rightInside && !leftInside) return right;
	}

	return null;
}
