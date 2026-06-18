#!/usr/bin/env bash
set -euo pipefail

# Split a WGS .nvl alias file into smaller alias files for parallel prefetch.

if [ "$#" -ne 1 ]; then
  echo "Usage: CHUNK_SIZE=45000 GROUP=bacteria bash scripts/split_nvl.sh INPUT.nvl" >&2
  exit 2
fi

alias_file="$1"
chunk_size="${CHUNK_SIZE:-45000}"
group="${GROUP:-$(basename "$alias_file" .nvl)}"
tmp_dir="${TMPDIR:-/tmp}/pksclusterdb_split_nvl_$$"

if ! grep -q '^VDBLIST ' "$alias_file"; then
  echo "No VDBLIST line found in $alias_file" >&2
  exit 1
fi

mkdir -p "$tmp_dir"
trap 'rm -rf "$tmp_dir"' EXIT

header_end_line=$(grep -n '^VDBLIST' "$alias_file" | cut -d: -f1)
head -n $((header_end_line - 1)) "$alias_file" > "$tmp_dir/header.txt"
sed -n "${header_end_line}p" "$alias_file" | sed 's/^VDBLIST //' | tr ' ' '\n' | sed '/^$/d' > "$tmp_dir/accessions.txt"

split -d -a 2 -l "$chunk_size" "$tmp_dir/accessions.txt" "${group}_chunk_"

for chunk_file in "${group}"_chunk_*; do
  chunk_num=$(printf "%s" "$chunk_file" | grep -o '[0-9]\+$')
  {
    cat "$tmp_dir/header.txt"
    echo -n "VDBLIST "
    tr '\n' ' ' < "$chunk_file"
    echo
  } > "${alias_file%.nvl}-${chunk_num}.nvl"
  rm "$chunk_file"
done

echo "Wrote split alias files matching ${alias_file%.nvl}-*.nvl"
