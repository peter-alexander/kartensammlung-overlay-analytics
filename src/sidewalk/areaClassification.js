export function classifyMeasurementContext(measurement) {
	const cornerScore = clamp01(
		measurement.directionSpreadDeg / 60 * 0.7 +
		measurement.directionDeltaDeg / 90 * 0.3
	);

	const plazaScore = clamp01(
		widthRamp(measurement.widthM, 4.0, 10.0) * 0.65 +
		(measurement.directionSpreadDeg / 60) * 0.25 +
		Math.min(measurement.sisPartsCount, 6) / 6 * 0.10
	);

	const complexityScore = clamp01(
		cornerScore * 0.45 +
		plazaScore * 0.30 +
		Math.min(measurement.sisPartsCount, 8) / 8 * 0.25
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

	const cornerScore = Math.max(group.maxCornerScore, avgCornerScore * 0.75);
	const plazaScore = Math.max(group.maxPlazaScore, avgPlazaScore * 0.85);
	const complexityScore = Math.max(group.maxComplexityScore, avgComplexityScore * 0.85);

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
