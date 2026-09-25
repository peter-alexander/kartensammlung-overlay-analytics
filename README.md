# kartensammlung-overlay-analytics

Neue JavaScript-Pipeline für Gehsteigbreiten der Kartensammlung.

## Methode

Die Pipeline verwendet FMZK-Gehsteigflächen als Richtungs- und Strukturquelle und SIS-Belagsflächen als tatsächliche Messfläche.

- FMZK wird featureweise gelesen und verarbeitet.
- SIS wird vorher lokal in eine SQLite/RTree-Datei indiziert.
- Pro FMZK-Feature werden Messstationen entlang des Außenrings erzeugt.
- Die lokale Gehrichtung wird geglättet und Ausreißer werden korrigiert.
- Quer zur lokalen Gehrichtung wird in den zugeordneten SIS-Flächen gemessen.
- Mehrere parallele Gehflächen pro Messstation bleiben erhalten.
- Benachbarte Messungen gleicher Breitenklasse werden sofort zu Liniensegmenten zusammengefasst.
- Am Ende wird ein einziges PMTiles-File gebaut.

## Breitenklassen

- `0`: `< 1.2 m`
- `1`: `>= 1.2 m - < 1.5 m`
- `2`: `>= 1.5 m - < 2.0 m`
- `3`: `>= 2.0 m - < 3.0 m`
- `4`: `>= 3.0 m`

## Lokaler Start

```bash
npm ci
npm run build:sidewalk-widths -- --config config/sidewalk-widths.config.js
```

## Datenquellen

Die finalen URLs sind in `config/sidewalk-widths.config.js` vorbereitet, aber auskommentiert.
Aktiv sind Test-URLs für kleine vorbereitete Ausschnitte auf Lima-City.


## Kopfsteinpflaster auf dem Hauptradverkehrsnetz

Die Analyse unter hrvn_cobblestone/ verschneidet das Wiener Hauptradverkehrsnetz (B/G/E) mit den SIS-Kopfsteinpflastertypen GM, GO und KL. SIS wird entsprechend der Empfehlung der Stadt Wien kachelweise geladen und bereits serverseitig per CQL auf die relevanten Belagsarten eingeschränkt.

Details, Parameter und Outputs: hrvn_cobblestone/README.md
