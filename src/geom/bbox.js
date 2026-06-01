export function bboxOfGeometry(geometry) {
	if (!geometry) return null;

	let minX = Infinity;
	let minY = Infinity;
	let maxX = -Infinity;
	let maxY = -Infinity;

	function visitCoord(coord) {
		const x = Number(coord[0]);
		const y = Number(coord[1]);
		if (!Number.isFinite(x) || !Number.isFinite(y)) return;
		minX = Math.min(minX, x);
		minY = Math.min(minY, y);
		maxX = Math.max(maxX, x);
		maxY = Math.max(maxY, y);
	}

	function visitCoords(coords) {
		if (!Array.isArray(coords)) return;
		if (typeof coords[0] === 'number') {
			visitCoord(coords);
			return;
		}
		for (const item of coords) visitCoords(item);
	}

	if (geometry.type === 'GeometryCollection') {
		for (const g of geometry.geometries || []) {
			const b = bboxOfGeometry(g);
			if (!b) continue;
			minX = Math.min(minX, b.minX);
			minY = Math.min(minY, b.minY);
			maxX = Math.max(maxX, b.maxX);
			maxY = Math.max(maxY, b.maxY);
		}
	} else {
		visitCoords(geometry.coordinates);
	}

	if (!Number.isFinite(minX)) return null;
	return { minX, minY, maxX, maxY };
}

export function expandBbox(bbox, byM) {
	return {
		minX: bbox.minX - byM,
		minY: bbox.minY - byM,
		maxX: bbox.maxX + byM,
		maxY: bbox.maxY + byM
	};
}
