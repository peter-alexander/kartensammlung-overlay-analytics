const DEG = 180 / Math.PI;
const RAD = Math.PI / 180;

export function axisAngleDeg(tx, ty) {
	let deg = Math.atan2(ty, tx) * DEG;
	deg = ((deg % 180) + 180) % 180;
	return deg;
}

export function axisDeltaDeg(a, b) {
	const d = Math.abs(a - b) % 180;
	return Math.min(d, 180 - d);
}

export function smoothAxisAngles(rawAnglesDeg, windowStations) {
	const out = [];

	for (let i = 0; i < rawAnglesDeg.length; i++) {
		let sx = 0;
		let sy = 0;
		let weightSum = 0;

		for (let j = Math.max(0, i - windowStations); j <= Math.min(rawAnglesDeg.length - 1, i + windowStations); j++) {
			const d = Math.abs(i - j);
			const weight = windowStations + 1 - d;
			const theta = rawAnglesDeg[j] * RAD * 2;
			sx += Math.cos(theta) * weight;
			sy += Math.sin(theta) * weight;
			weightSum += weight;
		}

		if (weightSum <= 0 || (Math.abs(sx) < 1e-12 && Math.abs(sy) < 1e-12)) {
			out.push(rawAnglesDeg[i]);
			continue;
		}

		let smoothed = 0.5 * Math.atan2(sy, sx) * DEG;
		smoothed = ((smoothed % 180) + 180) % 180;
		out.push(smoothed);
	}

	return out;
}

export function classifyDirection(rawDeg, usedDeg, neighborSpreadDeg, cfg) {
	const delta = axisDeltaDeg(rawDeg, usedDeg);
	if (neighborSpreadDeg > 45) return { quality: 'uncertain_direction', deltaDeg: delta };
	if (delta >= cfg.directionHardOutlierDeg && neighborSpreadDeg <= cfg.directionStableSpreadDeg) {
		return { quality: 'direction_hard_smoothed_outlier', deltaDeg: delta };
	}
	if (delta >= cfg.directionOutlierDeg && neighborSpreadDeg <= cfg.directionStableSpreadDeg) {
		return { quality: 'direction_smoothed_outlier', deltaDeg: delta };
	}
	if (delta >= cfg.directionOutlierDeg) {
		return { quality: 'direction_smoothed', deltaDeg: delta };
	}
	return { quality: 'ok', deltaDeg: delta };
}

export function localAxisSpread(rawAnglesDeg, index, windowStations) {
	const values = [];
	for (let j = Math.max(0, index - windowStations); j <= Math.min(rawAnglesDeg.length - 1, index + windowStations); j++) {
		values.push(rawAnglesDeg[j]);
	}
	if (values.length <= 1) return 0;

	let best = Infinity;
	for (const candidate of values) {
		const deltas = values.map((v) => axisDeltaDeg(v, candidate)).sort((a, b) => a - b);
		const p90 = deltas[Math.floor((deltas.length - 1) * 0.9)];
		best = Math.min(best, p90);
	}
	return best;
}
