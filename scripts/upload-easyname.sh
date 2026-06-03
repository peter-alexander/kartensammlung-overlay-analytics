#!/usr/bin/env bash
set -euo pipefail

: "${EASYNAME_FTP_HOST:?EASYNAME_FTP_HOST fehlt}"
: "${EASYNAME_FTP_USER:?EASYNAME_FTP_USER fehlt}"
: "${EASYNAME_FTP_PASSWORD:?EASYNAME_FTP_PASSWORD fehlt}"

REMOTE_DIR="${1:-Gehsteigbreiten}"
PMTILES_FILE="dist/sidewalk_widths.pmtiles"
SUMMARY_FILE="dist/sidewalk_widths_summary.json"

test -f "$PMTILES_FILE"
test -f "$SUMMARY_FILE"

lftp -u "$EASYNAME_FTP_USER","$EASYNAME_FTP_PASSWORD" "$EASYNAME_FTP_HOST" <<EOF
set ftp:ssl-allow no
mkdir -p "$REMOTE_DIR"
cd "$REMOTE_DIR"
put "$PMTILES_FILE" -o sidewalk_widths.pmtiles
put "$SUMMARY_FILE" -o sidewalk_widths_summary.json
bye
EOF
