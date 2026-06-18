#!/usr/bin/env bash
set -euo pipefail

if [ -n "${BLAST_BIN:-}" ]; then
  export PATH="$BLAST_BIN:$PATH"
fi

command -v tblastn >/dev/null || { echo "tblastn is not available on PATH" >&2; exit 1; }

THREADS="${THREADS:-8}"
QUERY="${QUERY:-queries/KSSignatureConsensusPKSDB.fasta}"
EVALUE="${EVALUE:-0.0001}"
OUT_DIR="${OUT_DIR:-outputs/tblastn}"
DB_ROOT="${DB_ROOT:-data/ncbi_blastdbs}"
DB_NAMES="${DB_NAMES:-ref_prok_rep_genomes ref_euk_rep_genomes env_nt patnt tsa_nt nt}"

mkdir -p "$OUT_DIR"

SECONDS=0
echo "Starting TBLASTN pipeline at $(date)"

process_db_sharded() {
  DBNAME=$1
  OUT="$OUT_DIR/tblastn_${DBNAME}_${EVALUE}.out"

  echo "Processing ${DBNAME}..."

  tblastn \
    -query "$QUERY" \
    -db "$DB_ROOT/${DBNAME}/${DBNAME}" \
    -out "$OUT" \
    -evalue "$EVALUE" \
    -outfmt 6 \
    -max_target_seqs 50000 \
    -num_threads "$THREADS"

  echo "Finished ${DBNAME}."
  echo ""
}

for DBNAME in $DB_NAMES; do
  process_db_sharded "$DBNAME"
done

echo "Finished at: $(date)"
echo "Total elapsed time: $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m $((SECONDS % 60))s"
