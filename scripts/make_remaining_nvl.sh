#!/usr/bin/env bash
set -euo pipefail

# Build a new .nvl alias file that excludes WGS/SRA accessions already present
# in the current working directory.

if [ "$#" -ne 2 ]; then
  echo "Usage: bash scripts/make_remaining_nvl.sh INPUT.nvl OUTPUT.nvl" >&2
  exit 2
fi

input_nvl="$1"
output_nvl="$2"
tmp_dir="${TMPDIR:-/tmp}/pksclusterdb_nvl_$$"

if ! grep -q '^VDBLIST ' "$input_nvl"; then
  echo "No VDBLIST line found in $input_nvl" >&2
  exit 1
fi

mkdir -p "$tmp_dir"
trap 'rm -rf "$tmp_dir"' EXIT

grep 'VDBLIST' "$input_nvl" | sed 's/^VDBLIST //' | tr ' ' '\n' | sed '/^$/d' | sort -u > "$tmp_dir/all_entries.txt"
ls -1 | grep -E '^[A-Z0-9]{5,}$' | sort -u > "$tmp_dir/downloaded_entries.txt" || true
comm -23 "$tmp_dir/all_entries.txt" "$tmp_dir/downloaded_entries.txt" > "$tmp_dir/remaining_entries.txt"

{
  grep -v '^VDBLIST ' "$input_nvl" || true
  echo -n "VDBLIST "
  tr '\n' ' ' < "$tmp_dir/remaining_entries.txt"
  echo
} > "$output_nvl"

echo "Wrote $output_nvl with $(wc -l < "$tmp_dir/remaining_entries.txt") remaining accessions"
