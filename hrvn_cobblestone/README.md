# Kopfsteinpflaster auf dem Hauptradverkehrsnetz

Diese Analyse findet SIS-Kopfsteinpflasterflächen im Korridor des Wiener Hauptradverkehrsnetzes (HRVN).

## Datenquellen

- HRVN: ogdwien:RADNETZOGD
- SIS-Belagsflächen: ogdwien:SISBELAGOGD

HRVN wird wie im Kartensammlungs-Overlay hrvn auf M18_RANG_SUB = B, G, E gefiltert.

Als gesuchtes Kopfsteinpflaster gelten die SIS-Belagsarten:

- GM – Großsteinpflaster mit Fugenverguss
- GO – Großsteinpflaster ohne Fugenverguss

## Methode

1. HRVN laden und auf B, G, E filtern.
2. In EPSG:31256 einen konfigurierbaren Korridor um das Netz bilden. Standard: 15 m je Seite.
3. Den HRVN-Gesamtbereich in maximal 2×2-km-Kacheln teilen.
4. SIS je Kachel per WFS laden. cql_filter beschränkt die Serverantwort bereits auf BELAG = GM oder GO.
5. Die Flächen jeder Kachel an der Kachelgrenze schneiden und vereinigen; dadurch entstehen an Kachelgrenzen keine doppelten Flächen.
6. Alle geladenen Kopfsteinpflasterflächen vereinigen.
7. Die Gesamtfläche mit dem HRVN-Korridor verschneiden.
8. Zusätzlich die zusammengeführten GM/GO-Flächen mit dem ungepufferten HRVN schneiden. Dadurch entsteht ein Linienlayer nur auf den tatsächlichen B/G/E-Trassen.
9. Debug- und Ergebnisgeometrien topologieerhaltend vereinfachen und nach EPSG:4326 schreiben.

Die Stadt Wien empfiehlt für SIS ausdrücklich kleine WFS-Ausschnitte, beispielsweise 2×2 km. Deshalb wird kein Wien-Gesamtdownload verwendet.

## Installation

    python -m pip install -r requirements-hrvn-cobblestone.txt

## Start

    python hrvn_cobblestone/run_hrvn_cobblestone.py

Vorhandene SIS-Kacheln werden wiederverwendet. Für einen kompletten Neuabruf:

    python hrvn_cobblestone/run_hrvn_cobblestone.py --refresh

Für einen kurzen technischen Test können nur die ersten Kacheln verarbeitet werden:

    python hrvn_cobblestone/run_hrvn_cobblestone.py --limit-tiles 2

## Outputs

- dist/hrvn-cobblestone/hrvn_cobblestone.geojson – vereinfachte Kopfsteinpflasterflächen innerhalb des HRVN-Korridors
- dist/hrvn-cobblestone/hrvn_cobblestone.fgb – derselbe Ergebnisdatensatz als FlatGeobuf
- dist/hrvn-cobblestone/debug_sis_cobblestone_merged.geojson – alle in den SIS-Kacheln gefundenen und zusammengeführten Kopfsteinpflasterflächen vor dem HRVN-Schnitt
- dist/hrvn-cobblestone/debug_sis_cobblestone_merged.fgb – derselbe Debug-Datensatz als FlatGeobuf
- dist/hrvn-cobblestone/debug_hrvn_corridor.geojson – verwendeter HRVN-Puffer
- dist/hrvn-cobblestone/debug_sis_cobblestone_merged_hrvn_unbuffered.geojson – GM/GO-Flächen mit dem ungepufferten HRVN B/G/E verschnitten; Linienlayer ohne seitlichen Straßenpuffer
- dist/hrvn-cobblestone/summary.json – Parameter, Flächenstatistik und Kachelstatistik

## Wichtige Parameter

In config/hrvn-cobblestone.json:

- analysis.buffer_m: Straßenkorridor je Seite; Standard 15 m
- download.tile_size_m: SIS-Kachelgröße; Standard 2000 m
- output.simplify_tolerance_m: Vereinfachung des Endergebnisses; Standard 0,2 m
- sources.sis.filter_property / filter_values: aktuell BELAG = GM, GO

Der Debug-Datensatz umfasst den rechteckigen Gesamtbereich des gefilterten HRVN. Er wird bewusst vor dem HRVN-Schnitt geschrieben, damit sich fehlende oder fälschlich erfasste Pflasterflächen gegenüber dem Korridor kontrollieren lassen.

Der ungepufferte Debug-Layer enthält pro Liniensegment `rank` (B/G/E) und `length_m` in Metern. Er eignet sich insbesondere, um Gehsteige, Parkspuren und andere seitlich zum HRVN liegende GM/GO-Flächen aus der Betrachtung auszuschließen.
