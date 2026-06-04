#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <local-path> <remote-path>" >&2
	exit 1
fi

: "${EASYNAME_FTP_HOST:?EASYNAME_FTP_HOST is required}"
: "${EASYNAME_FTP_USER:?EASYNAME_FTP_USER is required}"
: "${EASYNAME_FTP_PASSWORD:?EASYNAME_FTP_PASSWORD is required}"

LOCAL_PATH="$1"
REMOTE_PATH="$2"
REMOTE_PATH="${REMOTE_PATH#/}"

if [ ! -e "$LOCAL_PATH" ]; then
	echo "Local path does not exist: $LOCAL_PATH" >&2
	exit 1
fi

if ! command -v lftp >/dev/null 2>&1; then
	echo "lftp ist nicht installiert." >&2
	exit 1
fi

SSL_VERIFY="${LFTP_SSL_VERIFY:-no}"
PARALLEL="${LFTP_PARALLEL:-4}"
TIMEOUT="${LFTP_TIMEOUT:-20}"
MAX_RETRIES="${LFTP_MAX_RETRIES:-2}"

LFTP_CMDS=$(mktemp)
trap 'rm -f "$LFTP_CMDS"' EXIT

{
	echo "set net:timeout $TIMEOUT"
	echo "set net:max-retries $MAX_RETRIES"
	echo "set net:reconnect-interval-base 5"
	echo "set net:reconnect-interval-max 10"
	echo "set ftp:ssl-force true"
	echo "set ftp:passive-mode true"
	echo "set ssl:verify-certificate $SSL_VERIFY"
	echo "set cmd:fail-exit true"
	echo "open \"$EASYNAME_FTP_HOST\""
	echo "user \"$EASYNAME_FTP_USER\" \"$EASYNAME_FTP_PASSWORD\""

	if [ -d "$LOCAL_PATH" ]; then
		echo "mirror -R --delete --verbose --parallel=$PARALLEL \"$LOCAL_PATH\" \"$REMOTE_PATH\""
	else
		REMOTE_DIR="$(dirname -- "$REMOTE_PATH")"
		REMOTE_FILE="$(basename -- "$REMOTE_PATH")"
		if [ "$REMOTE_DIR" != "." ]; then
			echo "cd \"$REMOTE_DIR\""
		fi
		echo "put \"$LOCAL_PATH\" -o \"$REMOTE_FILE\""
	fi

	echo 'bye'
} > "$LFTP_CMDS"

lftp -f "$LFTP_CMDS"
