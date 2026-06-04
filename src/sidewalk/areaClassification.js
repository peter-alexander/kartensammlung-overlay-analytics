export function classifyMeasurementContext(measurement) {
	const spreadComponent = ramp(measurement.directionSpreadDeg, 35, 90);
	const deltaComponent = ramp(measurement.directionDeltaDeg, 35, 90);
	const widthComponent = widthRamp(measurement.widthM, 4.0, 10.0);
	const partsComponent = Math.min(measurement.sisPartsCount, 6) / 6;

	const cornerScore = clamp01(
		deltaComponent * 0.65 +
		spreadComponent * deltaComponent * 0.30 +
		spreadComponent * 0.05
	);

	const plazaScore = clamp01(
		widthComponent * 0.65 +
		spreadComponent * 0.20 +
		partsComponent * 0.15
	);

	const complexityScore = clamp01(
		cornerScore * 0.35 +
		plazaScore * 0.35 +
		partsComponent * 0.20 +
		spreadComponent * 0.10
	);

	return {
		cornerScore: round(cornerScore, 3),
		plazaScore: round(plazaScore, 3),
		complexityScore: round(complexityScore, 3),
		areaClass: classifyArea({ cornerScore, plazaScore, complexityScore })
	};
}

export function aggregateSegmentClassification(group) {
	const count = Math.max(1, group.count);
	const avgCornerScore = group.sumCornerScore / count;
	const avgPlazaScore = group.sumPlazaScore / count;
	const avgComplexityScore = group.sumComplexityScore / count;

	const cornerScore = Math.max(avgCornerScore, group.maxCornerScore * 0.35);
	const plazaScore = Math.max(avgPlazaScore, group.maxPlazaScore * 0.65);
	const complexityScore = Math.max(avgComplexityScore, group.maxComplexityScore * 0.50);

	return {
		areaClass: classifyArea({ cornerScore, plazaScore, complexityScore }),
		cornerScore: round(cornerScore, 3),
		plazaScore: round(plazaScore, 3),
		complexityScore: round(complexityScore, 3)
	};
}

function classifyArea({ cornerScore, plazaScore, complexityScore }) {
	if (plazaScore >= 0.62) return 'plaza_candidate';
	if (cornerScore >= 0.62) return 'corner_candidate';
	if (complexityScore >= 0.58) return 'complex_area';
	if (complexityScore >= 0.38 || cornerScore >= 0.45 || plazaScore >= 0.45) return 'uncertain';
	return 'normal_sidewalk';
}

function ramp(value, min, max) {
	if (value <= min) return 0;
	if (value >= max) return 1;
	return (value - min) / (max - min);
}

function widthRamp(widthM, minM, maxM) {
	if (widthM <= minM) return 0;
	if (widthM >= maxM) return 1;
	return (widthM - minM) / (maxM - minM);
}

function clamp01(value) {
	return Math.max(0, Math.min(1, Number.isFinite(value) ? value : 0));
}

function round(value, digits) {
	const factor = 10 ** digits;
	return Math.round(value * factor) / factor;
}
