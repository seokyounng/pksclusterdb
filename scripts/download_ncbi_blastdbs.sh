#!/usr/bin/env bash
set -euo pipefail

# Download and extract preformatted NCBI BLAST nucleotide databases.

command -v wget >/dev/null || { echo "wget is not available on PATH" >&2; exit 1; }
command -v tar >/dev/null || { echo "tar is not available on PATH" >&2; exit 1; }

DB_ROOT="${DB_ROOT:-data/ncbi_blastdbs}"
NCBI_BLASTDB_URL="${NCBI_BLASTDB_URL:-https://ftp.ncbi.nih.gov/blast/db}"

if [ "$#" -eq 0 ]; then
  set -- ref_prok_rep_genomes:00-24 ref_euk_rep_genomes:000-161 env_nt patnt:00-07 tsa_nt:00-03 nt:000-223
fi

mkdir -p "$DB_ROOT"

download_archive() {
  archive=$1
  archive_dir=$2
  wget -c -P "$archive_dir" "${NCBI_BLASTDB_URL%/}/$archive"
}

for spec in "$@"; do
  db=${spec%%:*}
  range=""
  if [ "$spec" != "$db" ]; then
    range=${spec#*:}
  fi

  db_dir="$DB_ROOT/$db"
  archive_dir="$DB_ROOT/zips_$db"
  mkdir -p "$db_dir" "$archive_dir"

  echo "Downloading $db archives into $archive_dir"
  if [ -n "$range" ]; then
    start=${range%-*}
    end=${range#*-}
    width=${#start}
    for ((n=10#$start; n<=10#$end; n++)); do
      printf -v i "%0${width}d" "$n"
      download_archive "${db}.${i}.tar.gz" "$archive_dir"
    done
  else
    download_archive "${db}.tar.gz" "$archive_dir"
  fi

  echo "Extracting $db into $db_dir"
  for archive in "$archive_dir"/"${db}"*.tar.gz; do
    [ -e "$archive" ] || continue
    tar -xzf "$archive" -C "$db_dir"
  done

  echo "Done: $db_dir"
done
