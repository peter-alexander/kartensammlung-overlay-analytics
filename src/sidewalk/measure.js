import { bboxOfGeometry, expandBbox } from '../geom/bbox.js';
import { extractPolygons, pointInPolygonCoords } from '../geom/polygon.js';
import { buildRingStations } from '../geom/stations.js';
import { dist, lerp, ensureClosedRing } from '../geom/basic.js';
import { lineGeometryInsideIntervals, mergeIntervals } from '../geom/linePolygonIntervals.js';
import { querySisIndex } from '../io/sisIndex.js';
import { widthClass } from './classes.js';
import { groupMeasurements, linePassesPointValidator } from './segments.js';

export function measureFmzkFeature({ feature, sisDb, config }) {
	const geometry = feature.geometry;
	if (!geometry) return { features: [], stats: emptyStats() };

	const stats = emptyStats();
	stats.fmzkFeatures = 1;

	const bbox = bboxOfGeometry(geometry);
	if (!bbox) return { features: [], stats };

	const sisCandidates = querySisIndex(
		sisDb,
		expandBbox(bbox, config.measurement.fmzkSisSearchBufferM)
	);
	stats.sisCandidates = sisCandidates.length;

	if (!sisCandidates.length) return { features: [], stats };

	const out = [];
	const polygons = extractPolygons(geometry);
	stats.fmzkPolygons = polygons.length;

	const pointValidator = (point) => pointInsideAnySis(point, sisCandidates);
	const connectionValidator = (a, b) => linePassesPointValidator(
		a,
		b,
		pointValidator,
		config.measurement.outputValidationSampleStepM || 0.5
	);

	for (const polygonCoords of polygons) {
		if (!polygonCoords?.length) continue;
		const outerRing = ensureClosedRing(polygonCoords[0]);
		const stations = buildRingStations(outerRing, polygonCoords, config.measurement);
		stats.stations += stations.length;

		const measurements = [];
		let previousBands = [];
		for (const station of stations) {
			const stationMeasurementsRaw = measureStation({ station, sisCandidates, config });
			const tracked = trackBands({
				measurements: stationMeasurementsRaw,
				previousBands,
				maxDistanceM: config.measurement.bandTrackingMaxDistanceM || 2.0,
				stats
			});
			previousBands = tracked.nextBands;

			stats.crossSections++;
			stats.measurements += tracked.measurements.length;
			measurements.push(...tracked.measurements);
		}

		const grouped = groupMeasurements({
			measurements,
			config,
			fmzkProperties: feature.properties || {},
			connectionValidator
		});
		out.push(...grouped.features);
		mergeStats(stats, grouped.stats);
	}

	return { features: out, stats };
}

function measureStation({ station, sisCandidates, config }) {
	const cfg = config.measurement;
	const half = cfg.crossSectionHalfLengthM;
	const cx = station.point[0] + station.normalX * cfg.innerProbeM;
	const cy = station.point[1] + station.normalY * cfg.innerProbeM;

	const lineA = [cx - station.normalX * half, cy - station.normalY * half];
	const lineB = [cx + station.normalX * half, cy + station.normalY * half];
	const lineLen = dist(lineA, lineB);

	let intervals = [];
	let parts = 0;
	for (const candidate of sisCandidates) {
		const candidateIntervals = lineGeometryInsideIntervals(lineA, lineB, candidate.geometry);
		if (!candidateIntervals.length) continue;
		parts++;
		intervals.push(...candidateIntervals);
	}

	intervals = mergeIntervals(intervals);
	if (!intervals.length) return [];

	const measured = [];
	for (const interval of intervals) {
		const widthM = (interval.t1 - interval.t0) * lineLen;
		if (widthM < cfg.minWidthM) continue;

		const midT = (interval.t0 + interval.t1) * 0.5;
		const mid = lerp(lineA, lineB, midT);
		const dx = mid[0] - station.point[0];
		const dy = mid[1] - station.point[1];
		const inwardDistanceM = dx * station.normalX + dy * station.normalY;

		measured.push({
			point: mid,
			chainM: station.chainM,
			widthM,
			widthClass: widthClass(widthM, config.widthClasses),
			inwardDistanceM,
			directionQuality: station.directionQuality,
			directionRawDeg: round(station.rawDeg, 2),
			directionUsedDeg: round(station.usedDeg, 2),
			directionDeltaDeg: round(station.directionDeltaDeg, 2),
			directionSpreadDeg: round(station.directionSpreadDeg, 2),
			sisPartsCount: parts
		});
	}

	measured.sort((a, b) => a.inwardDistanceM - b.inwardDistanceM);
	return measured;
}

function trackBands({ measurements, previousBands, maxDistanceM, stats }) {
	if (!measurements.length) return { measurements: [], nextBands: [] };

	const assigned = [];
	const usedPrevious = new Set();
	let nextBandId = previousBands.reduce((max, band) => Math.max(max, band.bandIndex + 1), 0);

	for (const measurement of measurements) {
		let best = null;
		let bestDistance = Infinity;

		for (let i = 0; i < previousBands.length; i++) {
			if (usedPrevious.has(i)) continue;
			const d = dist(measurement.point, previousBands[i].point);
			if (d < bestDistance) {
				bestDistance = d;
				best = i;
			}
		}

		if (best !== null && bestDistance <= maxDistanceM) {
			const bandIndex = previousBands[best].bandIndex;
			usedPrevious.add(best);
			assigned.push({ ...measurement, bandIndex });
			stats.bandTrackingMatched++;
		} else {
			assigned.push({ ...measurement, bandIndex: nextBandId++ });
			stats.bandTrackingNewBands++;
		}
	}

	assigned.sort((a, b) => a.inwardDistanceM - b.inwardDistanceM);
	return {
		measurements: assigned,
		nextBands: assigned.map((m) => ({ bandIndex: m.bandIndex, point: m.point }))
	};
}

function pointInsideAnySis(point, sisCandidates) {
	for (const candidate of sisCandidates) {
		for (const polygonCoords of extractPolygons(candidate.geometry)) {
			if (pointInPolygonCoords(point, polygonCoords)) return true;
		}
	}
	return false;
}

function emptyStats() {
	return {
		fmzkFeatures: 0,
		fmzkPolygons: 0,
		sisCandidates: 0,
		stations: 0,
		crossSections: 0,
		measurements: 0,
		bandTrackingMatched: 0,
		bandTrackingNewBands: 0,
		groupsStarted: 0,
		outputSegmentsWritten: 0,
		droppedSinglePointGroups: 0,
		splitByWidthClass: 0,
		splitByStationGap: 0,
		splitByInvalidConnection: 0
	};
}

function mergeStats(target, source) {
	for (const [key, value] of Object.entries(source)) {
		if (typeof value === 'number') target[key] = (target[key] || 0) + value;
	}
}

function round(value, digits) {
	const f = 10 ** digits;
	return Math.round(value * f) / f;
}
