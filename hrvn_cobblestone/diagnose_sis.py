#!/usr/bin/env python3
import json
import requests

BASE = "https://data.wien.gv.at/daten/geo"
TYPE_NAME = "ogdwien:SISBELAGOGD"
BBOX = "0,340000,2000,342000,EPSG:31256"

def run(label, extra):
	params = {
		"service": "WFS",
		"version": "1.1.0",
		"request": "GetFeature",
		"typeName": TYPE_NAME,
		"srsName": "EPSG:31256",
		"outputFormat": "json",
		"maxFeatures": "5",
	}
	params.update(extra)
	r = requests.get(BASE, params=params, timeout=180)
	print("\n===", label, "===")
	print("URL:", r.url)
	print("status:", r.status_code)
	print("content-type:", r.headers.get("content-type"))
	print("bytes:", len(r.content))
	try:
		data = r.json()
	except Exception:
		print(r.text[:2000])
		return
	features = data.get("features", [])
	print("features:", len(features))
	for feature in features[:5]:
		print(json.dumps(feature.get("properties", {}), ensure_ascii=False, sort_keys=True))

run("BBOX only", {"BBOX": BBOX})
run("TYPE GM/GO/KL only", {"cql_filter": "TYPE IN ('GM','GO','KL')"})
run("CQL spatial only", {"cql_filter": "BBOX(SHAPE,0,340000,2000,342000,'EPSG:31256')"})
run("CQL spatial + TYPE", {"cql_filter": "BBOX(SHAPE,0,340000,2000,342000,'EPSG:31256') AND TYPE IN ('GM','GO','KL')"})
