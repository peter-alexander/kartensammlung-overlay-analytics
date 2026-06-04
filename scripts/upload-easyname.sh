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

if [[ "$EASYNAME_FTP_HOST" == ftp://* || "$EASYNAME_FTP_HOST" == ftps://* ]]; then
	LFTP_URL="$EASYNAME_FTP_HOST"
elif [[ "$EASYNAME_FTP_HOST" == *:990 ]]; then
	LFTP_URL="ftps://${EASYNAME_FTP_HOST}"
else
	LFTP_URL="ftp://${EASYNAME_FTP_HOST}"
fi

echo "Upload-Ziel: ${LFTP_URL}/${REMOTE_DIR}"
echo "PMTiles-Größe: $(du -h "$PMTILES_FILE" | cut -f1)"
echo "Summary-Größe: $(du -h "$SUMMARY_FILE" | cut -f1)"

timeout 180s lftp -d -u "$EASYNAME_FTP_USER","$EASYNAME_FTP_PASSWORD" "$LFTP_URL" <<EOF
set net:timeout 15
set net:max-retries 2
set net:reconnect-interval-base 5
set net:reconnect-interval-max 10
set ftp:passive-mode on
set ftp:ssl-allow yes
set ftp:ssl-force no
set ftp:ssl-protect-data yes
set ssl:verify-certificate no
set cmd:fail-exit yes
mkdir -p "$REMOTE_DIR"
cd "$REMOTE_DIR"
put "$PMTILES_FILE" -o sidewalk_widths.pmtiles
put "$SUMMARY_FILE" -o sidewalk_widths_summary.json
bye
EOF
