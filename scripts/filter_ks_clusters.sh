#!/usr/bin/env bash
set -euo pipefail

# Filter every BLAST outfmt 6 result under outputs/tblastn and merge candidate
# KS-cluster loci into outputs/ks_filter/clusters_all.txt.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

TBLASTN_DIR="${TBLASTN_DIR:-$REPO_DIR/outputs/tblastn}"
OUT_DIR="${OUT_DIR:-$REPO_DIR/outputs/ks_filter}"
MERGED_OUT="${MERGED_OUT:-$OUT_DIR/clusters_all.txt}"
KS_FILTER="${KS_FILTER:-$REPO_DIR/ks_filter.py}"

mkdir -p "$OUT_DIR"
: > "$MERGED_OUT"

find "$TBLASTN_DIR" -name "*.out" -print | sort | while read -r result; do
  base=${result#"$TBLASTN_DIR"/}
  base=${base%.out}
  base=${base//\//_}
  out_file="$OUT_DIR/${base}_clusters.txt"
  python "$KS_FILTER" "$result" > "$out_file"
  cat "$out_file" >> "$MERGED_OUT"
done

echo "Merged KS-filtered clusters: $MERGED_OUT"
