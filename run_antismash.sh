#!/usr/bin/env bash
set -euo pipefail

# Optional local/HPC environment setup. Keep installs outside this script.
if [ -n "${MODULE_SETUP:-}" ]; then
  eval "$MODULE_SETUP"
fi

if [ -n "${CONDA_ENV:-}" ]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
  elif [ -n "${CONDA_EXE:-}" ] && [ -x "$CONDA_EXE" ]; then
    eval "$("$CONDA_EXE" shell.bash hook)"
  else
    echo "CONDA_ENV is set, but conda was not found. Put conda on PATH or set CONDA_EXE." >&2
    exit 1
  fi
  conda activate "$CONDA_ENV"
fi

command -v antismash >/dev/null || { echo "antismash is not available on PATH" >&2; exit 1; }
command -v parallel >/dev/null || { echo "GNU parallel is not available on PATH" >&2; exit 1; }

if [ "$#" -lt 3 ]; then
  echo "Usage: bash run_antismash.sh GB_DIR_OR_FILE FASTA_DIR OUT_DIR [RUN_NAME]" >&2
  exit 2
fi

GB_DIR="${1%/}"
FASTA_DIR="${2%/}"
OUT_DIR="${3%/}"
RUN_NAME="${4:-antismash}"

mkdir -p log

JOB_ID="${SLURM_JOB_ID:-local}"
exec > "log/${JOB_ID}_${RUN_NAME}.out"
exec 2> "log/${JOB_ID}_${RUN_NAME}.err"

THREADS_PER_JOB="${THREADS_PER_JOB:-8}"
JOBS_AT_ONCE="${JOBS_AT_ONCE:-4}"
FORCE_RERUN="${FORCE_RERUN:-0}"
FORCE_IDS="${FORCE_IDS:-}"
PY_FETCH_SCRIPT="${PY_FETCH_SCRIPT:-get_genbank.py}"
LOCK_STALE_SECONDS="${LOCK_STALE_SECONDS:-86400}"
FUNGAL_ANTISMASH_SUPPORTED="${FUNGAL_ANTISMASH_SUPPORTED:-auto}"

if [ "$FUNGAL_ANTISMASH_SUPPORTED" = "auto" ]; then
  ANTISMASH_HELP="$(antismash --help-showall 2>&1 || antismash --help 2>&1 || true)"
  if grep -qi 'glimmerhmm' <<< "$ANTISMASH_HELP"; then
    FUNGAL_ANTISMASH_SUPPORTED=1
  else
    FUNGAL_ANTISMASH_SUPPORTED=0
  fi
  unset ANTISMASH_HELP
fi
echo "Fungal antiSMASH support: $FUNGAL_ANTISMASH_SUPPORTED"

# Distinguishable skip files (job-scoped)
SKIP_MISSING_SEQUENCE="log/${JOB_ID}_skipped_missing_or_short_sequence.txt"
SKIP_UNSUPPORTED="log/${JOB_ID}_skipped_unsupported_taxon.txt"
SKIP_ALREADY_DONE="log/${JOB_ID}_skipped_already_complete.txt"
SKIP_STALE_OUTPUT="log/${JOB_ID}_rerun_stale_or_forced_output.txt"
FAILED_RUNS="log/${JOB_ID}_failed_antismash.txt"
RUN_MANIFEST="log/${JOB_ID}_${RUN_NAME}_run_manifest.tsv"

: > "$SKIP_MISSING_SEQUENCE"
: > "$SKIP_UNSUPPORTED"
: > "$SKIP_ALREADY_DONE"
: > "$SKIP_STALE_OUTPUT"
: > "$FAILED_RUNS"
printf "accession\tinput_file\toutput_dir\ttaxon\tgene_finding\n" > "$RUN_MANIFEST"

mkdir -p "$OUT_DIR"

export THREADS_PER_JOB OUT_DIR FASTA_DIR FORCE_RERUN FORCE_IDS PY_FETCH_SCRIPT LOCK_STALE_SECONDS FUNGAL_ANTISMASH_SUPPORTED RUN_MANIFEST FAILED_RUNS \
       SKIP_MISSING_SEQUENCE SKIP_UNSUPPORTED SKIP_ALREADY_DONE SKIP_STALE_OUTPUT

start_time=$(date +%s)
INPUT_LIST="log/${JOB_ID}_${RUN_NAME}_gb_inputs.txt"

run_one() {
  set -euo pipefail

  if [ "$#" -ne 1 ] || [ -z "${1:-}" ]; then
    echo "run_one expected exactly one .gb path, got $#" >&2
    exit 2
  fi

  gbk="$1"
  base=$(basename "$gbk" .gb)
  outpath="${OUT_DIR%/}/$base"
  workpath="${outpath}.__running__.$$"
  fasta_file="$FASTA_DIR/$base.fasta"

  has_complete_antismash() {
    d="$1"
    # Conservative completion test: require index.html and at least one additional expected artifact
    [ -s "$d/index.html" ] && ( [ -s "$d/regions.js" ] || ls "$d"/region*.gbk >/dev/null 2>&1 )
  }

  complete_from_expected_input() {
    d="$1"
    expected_input="$2"
    json="$d/$base.json"
    [ -s "$json" ] || return 1
    python -c 'import json, os, sys; data=json.load(open(sys.argv[1])); input_file=data.get("input_file", ""); expected=sys.argv[2]; sys.exit(0 if input_file == expected or os.path.basename(input_file) == os.path.basename(expected) else 1)' "$json" "$expected_input" 2>/dev/null
  }

  force_this_id() {
    [ "$FORCE_RERUN" = "1" ] && return 0
    [ -n "$FORCE_IDS" ] || return 1
    printf '%s\n' "$FORCE_IDS" | tr ',[:space:]' '\n' | grep -Fxq "$base"
  }

  accession_lockdir="${OUT_DIR%/}/.$base.antismash.lock"
  while ! mkdir "$accession_lockdir" 2>/dev/null; do
    lock_age=$(python -c 'import os, sys, time; p=sys.argv[1]; print(int(time.time() - os.path.getmtime(p)) if os.path.exists(p) else 0)' "$accession_lockdir")
    if [ "$lock_age" -gt "$LOCK_STALE_SECONDS" ]; then
      echo ">>> Breaking stale lock for $base (age ${lock_age}s)"
      rm -rf "$accession_lockdir"
      continue
    fi
    echo ">>> Waiting for active run of $base to finish"
    sleep 30
  done
  printf "%s\t%s\t%s\n" "${SLURM_JOB_ID:-local}" "$(hostname)" "$$" > "$accession_lockdir/owner.tsv"
  cleanup_run_one() {
    rm -f "$accession_lockdir/owner.tsv" 2>/dev/null || true
    rmdir "$accession_lockdir" 2>/dev/null || true
  }
  trap cleanup_run_one EXIT

  seq_len=$(awk '
    /^ORIGIN/ {found=1; next}
    /^\/\// {exit}
    found {gsub(/[0-9 \t]/, ""); total += length($0)}
    END {print total + 0}
  ' "$gbk")
  input_file="$gbk"
  input_is_fasta=0

  if [ -z "$seq_len" ] || [ "$seq_len" -lt 10 ]; then
    echo "*** GenBank has no full sequence, trying FASTA for: $base"
    mkdir -p "$FASTA_DIR"

    if [ ! -f "$fasta_file" ]; then
      fetch_lockdir="${FASTA_DIR%/}/.$base.fasta.fetch.lock"
      own_fetch_lock=0
      while [ "$own_fetch_lock" = "0" ]; do
        if mkdir "$fetch_lockdir" 2>/dev/null; then
          own_fetch_lock=1
          break
        fi
        fetch_lock_age=$(python -c 'import os, sys, time; p=sys.argv[1]; print(int(time.time() - os.path.getmtime(p)) if os.path.exists(p) else 0)' "$fetch_lockdir")
        if [ "$fetch_lock_age" -gt "$LOCK_STALE_SECONDS" ]; then
          echo ">>> Breaking stale FASTA fetch lock for $base (age ${fetch_lock_age}s)"
          rm -rf "$fetch_lockdir"
          continue
        fi
        sleep 5
        [ -f "$fasta_file" ] && break
      done

      if [ "$own_fetch_lock" = "1" ] && [ ! -f "$fasta_file" ]; then
        id_file="$fetch_lockdir/${base}.ids.txt"
        printf "%s\n" "$base" > "$id_file"
        echo ">>> FASTA missing for $base; fetching with $PY_FETCH_SCRIPT"
        if ! python "$PY_FETCH_SCRIPT" "$id_file" fasta "$FASTA_DIR"; then
          echo "<<< FASTA fetch failed for $base" >&2
        fi
      fi

      [ "$own_fetch_lock" = "1" ] && rm -rf "$fetch_lockdir"
    fi

    if [ -f "$fasta_file" ]; then
      seq_len=$(awk '!/^>/ {gsub(/[[:space:]]/, ""); total += length($0)} END {print total + 0}' "$fasta_file")
      if [ -n "$seq_len" ] && [ "$seq_len" -ge 10 ]; then
        input_file="$fasta_file"
        input_is_fasta=1
      fi
    fi
  fi

  if [ -z "$seq_len" ] || [ "$seq_len" -lt 10 ]; then
    echo "<<< Skipping: $base (no valid sequence in GenBank or FASTA)" >&2
    echo "$base" >> "$SKIP_MISSING_SEQUENCE"
    exit 0
  fi

  if has_complete_antismash "$outpath"; then
    if force_this_id; then
      echo "$base	forced" >> "$SKIP_STALE_OUTPUT"
    elif complete_from_expected_input "$outpath" "$input_file"; then
      echo "$base" >> "$SKIP_ALREADY_DONE"
      exit 0
    else
      echo "$base	stale_or_unknown_input" >> "$SKIP_STALE_OUTPUT"
    fi
  fi

  taxonomy=$({ grep -A3 '^  ORGANISM' "$gbk" || true; } | tr -d '\n')
  tool_opts=()
  taxon=""
  gene_finding=""

  if echo "$taxonomy" | grep -qiE 'Bacteria|Archaea'; then
    echo ">>> Processing prokaryote: $base"
    tool_opts=(--taxon bacteria --genefinding-tool prodigal)
    taxon="bacteria"
    gene_finding="prodigal"

  elif echo "$taxonomy" | grep -qi 'Fungi'; then
    if [ "$FUNGAL_ANTISMASH_SUPPORTED" != "1" ]; then
      echo "<<< Skipping fungus in FASTA/prodigal workflow; current antiSMASH does not support glimmerhmm: $base" >&2
      echo "$base" >> "$SKIP_UNSUPPORTED"
      exit 0
    fi

    echo ">>> Processing fungus: $base"
    taxon="fungi"
    if [ "$input_is_fasta" = "0" ] && grep -q '^     CDS ' "$gbk"; then
      echo ">>> CDS features found, skipping gene finding"
      tool_opts=(--taxon fungi --genefinding-tool none)
      gene_finding="none"
    else
      echo ">>> No CDS features found, using glimmerhmm for gene finding"
      tool_opts=(--taxon fungi --genefinding-tool glimmerhmm)
      gene_finding="glimmerhmm"
    fi

  else
    echo "<<< Skipping unsupported eukaryote (non-fungal): $base" >&2
    echo "$base" >> "$SKIP_UNSUPPORTED"
    exit 0
  fi

  # Remove leftovers from a prior timeout for this accession; final output is only replaced after completion.
  for stale_workpath in "${OUT_DIR%/}/$base.__running__."*; do
    [ -e "$stale_workpath" ] || continue
    [ "$stale_workpath" = "$workpath" ] && continue
    rm -rf "$stale_workpath" || true
  done
  rm -rf "$workpath"
  mkdir -p "$workpath"

  printf "%s\t%s\t%s\t%s\t%s\n" "$base" "$input_file" "$outpath" "$taxon" "$gene_finding" >> "$RUN_MANIFEST"
  echo ">>> Running antiSMASH on $input_file -> $outpath"
  if ! antismash "$input_file" --cpus "$THREADS_PER_JOB" "${tool_opts[@]}" --output-dir "$workpath"; then
    echo "<<< antiSMASH failed for $base; see this job's .err log" >&2
    echo "$base" >> "$FAILED_RUNS"
    exit 0
  fi

  if ! has_complete_antismash "$workpath"; then
    echo "<<< antiSMASH finished but expected outputs are incomplete for $base" >&2
    echo "$base" >> "$FAILED_RUNS"
    exit 0
  fi

  rm -rf "$outpath"
  mv "$workpath" "$outpath"
  echo ">>> Done: $base"
}
export -f run_one

if [ -f "$GB_DIR" ]; then
  printf '%s\n' "$GB_DIR" > "$INPUT_LIST"
else
  find "$GB_DIR" \( -type f -o -type l \) -name "*.gb" | sort > "$INPUT_LIST"
fi

if [ ! -s "$INPUT_LIST" ]; then
  echo "No .gb inputs found from: $GB_DIR" >&2
  echo "If you passed one accession, use a path to an existing .gb file." >&2
  exit 2
fi

parallel -j "$JOBS_AT_ONCE" --halt soon,fail=1 run_one :::: "$INPUT_LIST"

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total elapsed wall time: $((runtime/3600))h $(( (runtime/60)%60 ))m $((runtime%60))s"

echo "Skip files:"
echo "  $SKIP_MISSING_SEQUENCE"
echo "  $SKIP_UNSUPPORTED"
echo "  $SKIP_ALREADY_DONE"
echo "  $SKIP_STALE_OUTPUT"
echo "Failed antiSMASH runs:"
echo "  $FAILED_RUNS"
echo "Run manifest:"
echo "  $RUN_MANIFEST"
