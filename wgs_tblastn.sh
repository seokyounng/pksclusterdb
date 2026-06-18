#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Run SRA/WGS prefetch and tblastn_vdb for one taxonomic group.

if [ -n "${BLAST_BIN:-}" ]; then
  export PATH="$BLAST_BIN:$PATH"
fi
if [ -n "${SRA_TOOLS_BIN:-}" ]; then
  export PATH="$SRA_TOOLS_BIN:$PATH"
fi
if [ -n "${MODULE_SETUP:-}" ]; then
  eval "$MODULE_SETUP"
fi

GROUP="${GROUP:-bacteria}"
case "$GROUP" in
  bacteria) default_taxid=2 ;;
  archaea) default_taxid=2157 ;;
  eukaryota) default_taxid=2759 ;;
  *) default_taxid="" ;;
esac

TAXID="${TAXID:-$default_taxid}"
if [ -z "$TAXID" ]; then
  echo "Set TAXID for GROUP=$GROUP" >&2
  exit 2
fi

QUERY="${QUERY:-queries/KSSignatureConsensusPKSDB.fasta}"
EVALUE="${EVALUE:-0.0001}"
THREADS="${THREADS:-4}"
PARALLEL_JOBS="${PARALLEL_JOBS:-4}"
PREFETCH_PARALLEL="${PREFETCH_PARALLEL:-4}"
CHUNK_SIZE="${CHUNK_SIZE:-5000}"
ALIAS_FILE="${ALIAS_FILE:-${GROUP}-wgs}"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/data/wgs_${GROUP}}"
TBLASTN_DIR="${TBLASTN_DIR:-$SCRIPT_DIR/outputs/tblastn/wgs_${GROUP}}"
CREATE_ALIAS="${CREATE_ALIAS:-0}"
RUN_PREFETCH="${RUN_PREFETCH:-1}"
RUN_TBLASTN="${RUN_TBLASTN:-1}"
REDOWNLOAD_FAILED="${REDOWNLOAD_FAILED:-1}"
MERGE_OUTPUT="${MERGE_OUTPUT:-0}"

case "$QUERY" in
  /*) ;;
  *) QUERY="$SCRIPT_DIR/$QUERY" ;;
esac
case "$TBLASTN_DIR" in
  /*) ;;
  *) TBLASTN_DIR="$SCRIPT_DIR/$TBLASTN_DIR" ;;
esac

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

command -v parallel >/dev/null || { echo "GNU parallel is not available on PATH" >&2; exit 1; }
if [ "$RUN_PREFETCH" = "1" ] || [ "$REDOWNLOAD_FAILED" = "1" ]; then
  command -v prefetch >/dev/null || { echo "prefetch is not available on PATH" >&2; exit 1; }
fi
if [ "$RUN_TBLASTN" = "1" ]; then
  command -v tblastn_vdb >/dev/null || { echo "tblastn_vdb is not available on PATH" >&2; exit 1; }
fi

SECONDS=0
echo "Starting WGS pipeline at $(date)"
echo "Group=$GROUP TaxID=$TAXID Alias=${ALIAS_FILE}.nvl"

if [ "$CREATE_ALIAS" = "1" ]; then
  command -v perl >/dev/null || { echo "perl is not available on PATH" >&2; exit 1; }
  title="$(printf '%s' "$GROUP" | tr '[:lower:]' '[:upper:]') WGS"
  perl "$SCRIPT_DIR/scripts/taxid2wgs.pl" -title "$title" -alias_file "$ALIAS_FILE" "$TAXID"
fi

if [ "$RUN_PREFETCH" = "1" ]; then
  if [ ! -s "${ALIAS_FILE}.nvl" ]; then
    echo "Missing alias file: ${ALIAS_FILE}.nvl" >&2
    echo "Create one with CREATE_ALIAS=1 before running prefetch." >&2
    exit 1
  fi

  echo "Splitting ${ALIAS_FILE}.nvl into chunks of $CHUNK_SIZE accessions"
  accessions_file="${GROUP}_accessions.txt"
  grep '^VDBLIST ' "${ALIAS_FILE}.nvl" | sed 's/^VDBLIST //' | tr ' ' '\n' | sed '/^$/d' > "$accessions_file"
  split -d -a 2 -l "$CHUNK_SIZE" "$accessions_file" "${GROUP}_part_"

  for part in "${GROUP}"_part_*; do
    [ -e "$part" ] || continue
    {
      echo "# Alias chunk"
      echo "TITLE ${GROUP} WGS CHUNK"
      echo -n "VDBLIST "
      tr '\n' ' ' < "$part"
      echo
    } > "${part}.nvl"
  done
  rm -f "$accessions_file"

  parts=( "${GROUP}"_part_*.nvl )
  if [ ! -e "${parts[0]}" ]; then
    echo "No alias chunks were created from ${ALIAS_FILE}.nvl" >&2
    exit 1
  fi

  echo "Prefetching ${#parts[@]} chunks with $PREFETCH_PARALLEL concurrent jobs"
  counter=0
  for part in "${parts[@]}"; do
    prefetch --option-file "$part" &
    counter=$((counter + 1))
    if (( counter % PREFETCH_PARALLEL == 0 )); then
      wait
    fi
  done
  wait

  rm -f "${GROUP}"_part_* "${GROUP}"_part_*.nvl
fi

if [ "$RUN_TBLASTN" = "1" ]; then
  mkdir -p "$TBLASTN_DIR"

  comm -23 \
    <(ls -1 | grep -E '^[A-Z0-9]{5,}$' | sort) \
    <(ls -1 "$TBLASTN_DIR" 2>/dev/null | sed 's/\.tblastn\.out$//' | sort) \
    | sed '/^$/d' > accessions.txt

  missing_count=$(wc -l < accessions.txt | tr -d ' ')
  echo "Found $missing_count accessions missing tblastn output"

  if [ "$missing_count" -gt 0 ]; then
    redownload_sra() {
      acc=$1
      echo "Attempting redownload of SRA data for $acc"
      rm -rf "$acc"
      prefetch "$acc"
    }
    export -f redownload_sra

    parallel -j "$PARALLEL_JOBS" --halt soon,fail=0 '
      acc={}
      out="'"${TBLASTN_DIR}"'/${acc}.tblastn.out"
      if tblastn_vdb -db "$acc" -query "'"${QUERY}"'" -evalue "'"${EVALUE}"'" -out "$out" -outfmt 6 -num_threads "'"${THREADS}"'"; then
        (flock -x 200; echo "$acc" >> done_ids.txt) 200>done_ids.txt.lock
      elif [ "'"${REDOWNLOAD_FAILED}"'" = "1" ] && bash -c "redownload_sra $acc" &&
           tblastn_vdb -db "$acc" -query "'"${QUERY}"'" -evalue "'"${EVALUE}"'" -out "$out" -outfmt 6 -num_threads "'"${THREADS}"'"; then
        (flock -x 200; echo "$acc" >> done_ids.txt) 200>done_ids.txt.lock
      else
        : > "$out"
        (flock -x 200; echo "$acc" >> failed_ids.txt) 200>failed_ids.txt.lock
      fi
    ' :::: accessions.txt
  fi

  rm -f accessions.txt
fi

if [ "$MERGE_OUTPUT" = "1" ]; then
  merged="${GROUP}_wgs_tblastn_merged_${EVALUE}.out"
  find "$TBLASTN_DIR" -name "*.tblastn.out" -print0 | xargs -0 cat > "$merged"
  echo "Merged results written to: $merged"
fi

echo "Finished at: $(date)"
echo "Total elapsed time: $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m $((SECONDS % 60))s"
