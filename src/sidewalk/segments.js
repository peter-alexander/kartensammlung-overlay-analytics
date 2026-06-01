import { simplifyLine } from '../geom/simplify.js';

export function groupMeasurements({ measurements, config, fmzkProperties }) {
	if (!measurements.length) return [];

	const byBand = new Map();
	for (const m of measurements) {
		if (m.widthClass === null) continue;
		const key = String(m.bandIndex);
		if (!byBand.has(key)) byBand.set(key, []);
		byBand.get(key).push(m);
	}

	const features = [];
	for (const [bandKey, items] of byBand) {
		items.sort((a, b) => a.chainM - b.chainM);
		let current = null;

		for (const item of items) {
			const startsNew = !current ||
				current.widthClass !== item.widthClass ||
				item.chainM - current.lastChainM > config.measurement.maxStationGapM;

			if (startsNew) {
				if (current) features.push(makeFeature(current, config, fmzkProperties));
				current = startGroup(item, Number(bandKey));
			} else {
				appendGroup(current, item);
			}
		}

		if (current) features.push(makeFeature(current, config, fmzkProperties));
	}

	return features;
}

function startGroup(item, bandIndex) {
	return {
		bandIndex,
		widthClass: item.widthClass,
		points: [item.point],
		count: 1,
		minWidthM: item.widthM,
		maxWidthM: item.widthM,
		sumWidthM: item.widthM,
		firstChainM: item.chainM,
		lastChainM: item.chainM,
		directionQualities: new Set([item.directionQuality]),
		maxDirectionDeltaDeg: item.directionDeltaDeg,
		maxDirectionSpreadDeg: item.directionSpreadDeg,
		maxSisPartsCount: item.sisPartsCount
	};
}

function appendGroup(group, item) {
	group.points.push(item.point);
	group.count++;
	group.minWidthM = Math.min(group.minWidthM, item.widthM);
	group.maxWidthM = Math.max(group.maxWidthM, item.widthM);
	group.sumWidthM += item.widthM;
	group.lastChainM = item.chainM;
	group.directionQualities.add(item.directionQuality);
	group.maxDirectionDeltaDeg = Math.max(group.maxDirectionDeltaDeg, item.directionDeltaDeg);
	group.maxDirectionSpreadDeg = Math.max(group.maxDirectionSpreadDeg, item.directionSpreadDeg);
	group.maxSisPartsCount = Math.max(group.maxSisPartsCount, item.sisPartsCount);
}

function makeFeature(group, config, fmzkProperties) {
	const coords = simplifyLine(group.points, config.measurement.simplifyToleranceM);
	if (coords.length < 2 && group.points.length === 1) {
		const p = group.points[0];
		coords.push([p[0] + 0.001, p[1] + 0.001]);
	}

	return {
		type: 'Feature',
		geometry: {
			type: 'LineString',
			coordinates: coords
		},
		properties: {
			method: config.method,
			width_class: group.widthClass,
			min_width_m: round(group.minWidthM, 2),
			max_width_m: round(group.maxWidthM, 2),
			avg_width_m: round(group.sumWidthM / group.count, 2),
			count: group.count,
			band_index: group.bandIndex,
			from_m: round(group.firstChainM, 2),
			to_m: round(group.lastChainM, 2),
			quality: summarizeQuality(group.directionQualities),
			max_direction_delta_deg: round(group.maxDirectionDeltaDeg, 2),
			max_direction_spread_deg: round(group.maxDirectionSpreadDeg, 2),
			max_sis_parts_count: group.maxSisPartsCount,
			fmzk_id: fmzkProperties.OBJECTID ?? fmzkProperties.ID ?? fmzkProperties.FID ?? null
		}
	};
}

function summarizeQuality(qualities) {
	if (qualities.has('uncertain_direction')) return 'uncertain_direction';
	if (qualities.has('direction_hard_smoothed_outlier')) return 'direction_hard_smoothed_outlier';
	if (qualities.has('direction_smoothed_outlier')) return 'direction_smoothed_outlier';
	if (qualities.has('direction_smoothed')) return 'direction_smoothed';
	return 'ok';
}

function round(value, digits) {
	const f = 10 ** digits;
	return Math.round(value * f) / f;
}
