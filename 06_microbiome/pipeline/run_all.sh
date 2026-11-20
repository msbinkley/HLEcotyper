#!/bin/bash
# =============================================================================
# run_all.sh  --  submit the FULL microbiome pipeline with SLURM dependencies
# =============================================================================
# Usage:   bash run_all.sh [maxconcurrent]
#   maxconcurrent : max array tasks running at once for steps 01 & 02 (default 8)
#
# Full workflow, each stage gated on the previous (afterok):
#   01 host removal (array)  ->  02 kraken2 (array)  ->  03 bracken+merge  ->  05 downstream R
#
# Resume vs clean rerun:
#   bash run_all.sh            -> submits the chain; step skip-guards reuse any
#                                 already-complete per-sample outputs (cheap resume).
#   FRESH=1 bash run_all.sh    -> clears ALL outputs first (01_host_removed + 02/03/05/06)
#                                 and reruns from scratch. REQUIRED when run-determining
#                                 parameters changed (host removal OR classification).
# Per-cohort overrides go in the ENVIRONMENT (see config.sh), e.g.:
#   FRESH=1 RUN_DIR=/oak/.../cHL_microbiome TRIM_DIR=/oak/.../cHL_fastp SAMPLE_MANIFEST= \
#       bash run_all.sh 8
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
cd "$SCRIPT_DIR"
mkdir -p "$LOG_DIR"   # per-cohort log dir (RUN_DIR/logs); SLURM .out/.err routed here via --chdir in submit()

MAXC="${1:-8}"

# Show the active cohort up front -- this pipeline is reusable across cohorts via
# environment overrides (see config.sh), so confirm you are running the RIGHT one.
echo "=============================================================="
echo " Microbiome pipeline (FULL: 01 -> 02 -> 03 -> 05)"
echo "   TRIM_DIR = $TRIM_DIR"
echo "   D_HOST   = $D_HOST"
echo "   RUN_DIR  = $RUN_DIR"
echo "   MANIFEST = ${SAMPLE_MANIFEST:-(none -> glob $TRIM_DIR)}"
echo "   REFS_DIR = $REFS_DIR"
echo "   KRAKEN2_CONFIDENCE=$KRAKEN2_CONFIDENCE  BRACKEN_THRESH=$BRACKEN_THRESH"
echo "=============================================================="

[ -d "$TRIM_DIR" ] || { echo "ERROR: TRIM_DIR does not exist: $TRIM_DIR" >&2; exit 1; }

# Run-determining parameters: ANY change invalidates existing outputs. Includes the
# host-removal knobs (step 01) AND the classification/abundance knobs (steps 02/03)
# so a change to either forces a clean rerun rather than silently reusing stale data
# (step skip-guards key only on file existence, not on parameters).
run_fingerprint=$(printf '%s\n' \
    "STAR_TWOPASS=$STAR_TWOPASS" \
    "STAR_FFPE_RELAXED=$STAR_FFPE_RELAXED" \
    "STAR_SEED_LMAX=$STAR_SEED_LMAX" \
    "STAR_MULTIMAP_NMAX=$STAR_MULTIMAP_NMAX" \
    "STAR_WINANCHOR_NMAX=$STAR_WINANCHOR_NMAX" \
    "LOWCOMPLEXITY_ENTROPY=$LOWCOMPLEXITY_ENTROPY" \
    "DEDUP=$DEDUP" \
    "T2T_BT2_INDEX=$T2T_BT2_INDEX" \
    "KRAKEN2_DB=$KRAKEN2_DB" \
    "KRAKEN2_CONFIDENCE=$KRAKEN2_CONFIDENCE" \
    "KRAKEN2_MIN_HIT_GROUPS=$KRAKEN2_MIN_HIT_GROUPS" \
    "REPORT_MINIMIZER_DATA=$REPORT_MINIMIZER_DATA" \
    "READ_LEN=$READ_LEN" \
    "BRACKEN_THRESH=$BRACKEN_THRESH" \
    "BRACKEN_LEVELS=$BRACKEN_LEVELS")
fingerprint_file="$RUN_DIR/.run_fingerprint"

ensure_run_subdir() {
    dir="$1"
    case "$dir/" in
        "$RUN_DIR"/*) ;;
        *) echo "ERROR: refusing to clean outside RUN_DIR: $dir (RUN_DIR=$RUN_DIR)" >&2; exit 1 ;;
    esac
}

outputs_present() {
    shopt -s nullglob
    local outputs=(
        "$D_HOST"/*_microbe_R1.fastq
        "$D_KRAKEN"/*.kreport*
        "$D_BRACKEN"/*
        "$D_RESULTS"/*
        "$D_FIGS"/*
    )
    shopt -u nullglob
    [ "${#outputs[@]}" -gt 0 ]
}

fingerprint_matches=0
if [ -s "$fingerprint_file" ] && [ "$(cat "$fingerprint_file")" = "$run_fingerprint" ]; then
    fingerprint_matches=1
fi

fingerprint_reason=""
if [ "${FRESH:-0}" = "1" ]; then
    fingerprint_reason="FRESH=1 requested"
elif [ "$fingerprint_matches" -eq 0 ]; then
    if [ -e "$fingerprint_file" ]; then
        fingerprint_reason="run-determining parameters changed since the previous run"
    elif outputs_present; then
        fingerprint_reason="no fingerprint exists but stage outputs are present"
    fi
fi

if [ -n "$fingerprint_reason" ]; then
    echo "WARNING: $fingerprint_reason." >&2
    echo "         Existing outputs under $RUN_DIR may be stale and would be reused by skip-guards." >&2
    if [ "${FRESH:-0}" != "1" ]; then
        echo "ERROR: refusing to delete existing outputs without explicit confirmation." >&2
        echo "       Re-run with: FRESH=1 bash run_all.sh $MAXC" >&2
        exit 1
    fi
    ensure_run_subdir "$D_HOST"
    ensure_run_subdir "$D_KRAKEN"
    ensure_run_subdir "$D_BRACKEN"
    ensure_run_subdir "$D_RESULTS"
    ensure_run_subdir "$D_FIGS"
    echo ">>> [fresh] clearing ALL stage outputs in $RUN_DIR (incl. host removal)" >&2
    rm -rf -- "$D_HOST"/*
    rm -f  -- "$D_KRAKEN"/*.kreport*
    rm -rf -- "$D_BRACKEN"/*
    rm -rf -- "$D_RESULTS"/*
    rm -rf -- "$D_FIGS"/*
fi
mkdir -p "$RUN_DIR"
printf '%s\n' "$run_fingerprint" > "$fingerprint_file"

# Build the authoritative sample sheet from the TRIMMED FASTQs (TRIM_DIR). Step 01
# reads its r1/r2 paths; step 02 reads only the sample NAME and derives the host-
# removed path itself, so the same sheet drives the whole chain.
echo ">>> [sample sheet] building from trimmed FASTQs in $TRIM_DIR"
bash "$SCRIPT_DIR/00_build_sample_sheet.sh"
[ -s "$SAMPLE_LIST" ] || { echo "ERROR: sample sheet is empty: $SAMPLE_LIST" >&2; exit 1; }
N=$(wc -l < "$SAMPLE_LIST")
[ "$N" -gt 0 ] || { echo "ERROR: no samples in $SAMPLE_LIST" >&2; exit 1; }
echo ">>> $N samples"

# Duplicate sample names would clobber each other's outputs (every stage keys on the
# bare sample name) and the step-03 completeness check could not detect it. Fail loudly.
dups=$(awk -F'\t' 'seen[$2]++{print $2}' "$SAMPLE_LIST" | sort -u)
if [ -n "$dups" ]; then
    echo "ERROR: duplicate sample names in $SAMPLE_LIST -- outputs would clobber:" >&2
    echo "$dups" | sed 's/^/       /' >&2
    exit 1
fi

# Pre-flight: verify modules/tools/DBs/inputs before submitting anything.
# Bypass with RUN_ALL_SKIP_PREFLIGHT=1 if you have already checked.
if [ "${RUN_ALL_SKIP_PREFLIGHT:-0}" != "1" ]; then
    echo ">>> [preflight] checking environment..."
    bash "$SCRIPT_DIR/preflight.sh" || {
        echo "ERROR: preflight failed. Fix the FAILs above, or bypass with" >&2
        echo "       RUN_ALL_SKIP_PREFLIGHT=1 bash run_all.sh" >&2
        exit 1; }
fi

# Export PIPELINE_DIR so each job can find config.sh regardless of SLURM's spool dir.
# --chdir="$RUN_DIR" routes each .sbatch's relative `--output=logs/...` into the
# COHORT's $RUN_DIR/logs. --parsable -> sbatch prints just the job id (cut strips any
# ";cluster" suffix) for the dependency chain.
submit() { sbatch --parsable --export=ALL,PIPELINE_DIR="$SCRIPT_DIR" --chdir="$RUN_DIR" "$@" | cut -d';' -f1; }

j01=$(submit --array=1-"$N"%"$MAXC" 01_host_removal.sbatch)
echo "    01 host removal     : $j01  (array 1-$N%$MAXC)"

# 02 starts only after EVERY host-removal task succeeds (afterok). A step-01 failure
# leaves 02/03/05 in DependencyNeverSatisfied and they cancel -- intentional, so a
# partial cohort is never classified. Fix the failed index and re-run.
j02=$(submit --dependency=afterok:"$j01" --array=1-"$N"%"$MAXC" 02_kraken2.sbatch)
echo "    02 kraken2          : $j02  (afterok $j01)"

j03=$(submit --dependency=afterok:"$j02" 03_bracken_merge.sbatch)
echo "    03 bracken+merge    : $j03  (afterok $j02)"

j05=$(submit --dependency=afterok:"$j03" 05_downstream.sbatch)
echo "    05 downstream R     : $j05  (afterok $j03)"

echo ">>> Submitted FULL pipeline. Track with:  squeue -u $USER"
echo ">>> Final tables -> $D_RESULTS ; figures -> $D_FIGS"
