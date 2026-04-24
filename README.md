# kartensammlung-overlay-analytics

Erste Benchmark-Struktur für Gehsteigbreiten.

## Start lokal

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/benchmark_sidewalk_widths.py --config config/sidewalk_benchmarks.json --output outputs/local-run
```

## Ziel der ersten Runde

- Wien-Daten direkt vom WFS holen
- nur outer rings messen
- geglättete lokale Richtung verwenden
- je Messstation genau ein Hauptsegment wählen
- JSON-Metriken und Debug-GeoJSON schreiben

## Nächster Schritt danach

Wenn die Benchmark-Runde sauber läuft und die Ausgabe optisch passt, kann die Struktur in `src/`-Module aufgeteilt und um `pyogrio` / `GeoPackage` erweitert werden.
