export function widthClass(widthM, classes) {
	for (const cls of classes) {
		const minOk = cls.min === null || widthM >= cls.min;
		const maxOk = cls.max === null || widthM < cls.max;
		if (minOk && maxOk) return cls.id;
	}
	return null;
}
