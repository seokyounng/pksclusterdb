#!/usr/bin/env bash
set -euo pipefail

# Download and extract MiBIG JSON/GBK archives.

command -v wget >/dev/null || { echo "wget is not available on PATH" >&2; exit 1; }
command -v tar >/dev/null || { echo "tar is not available on PATH" >&2; exit 1; }

MIBIG_DIR="${MIBIG_DIR:-data/mibig}"
MIBIG_VERSION="${MIBIG_VERSION:-4.0}"
MIBIG_URL="${MIBIG_URL:-https://dl.secondarymetabolites.org/mibig}"

mkdir -p "$MIBIG_DIR"

for kind in json gbk; do
  archive="mibig_${kind}_${MIBIG_VERSION}.tar.gz"
  wget -c -P "$MIBIG_DIR" "$MIBIG_URL/$archive"
  tar -xzf "$MIBIG_DIR/$archive" -C "$MIBIG_DIR"
done

echo "MiBIG archives extracted under $MIBIG_DIR"
