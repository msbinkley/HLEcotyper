#!/bin/bash
# =============================================================================
# sweep_confidence_e2e.sh -- end-to-end Kraken2 --confidence sweep
# =============================================================================
# Runs the FULL pipeline (steps 02 -> 03 -> 05) once PER confidence value, each
# into an isolated output dir, all sharing the SINGLE step-01 host removal (host
# removal is confidence-independent, so it runs once).
#
#   step 01 (host removal, ONCE)  ->  for cs in CS_LIST:
#                                        02 kraken2 @ confidence=cs
#                                        03 bracken+merge
#                                        05 aggregate/filter
#                                     -> $RUN_DIR/confidence_sweep_e2e/cs_<cs>/05_results
#
# How isolation works: each cs reuses run_all.sh with RUN_DIR pointed at a per-cs
# subdir (so D_KRAKEN/D_BRACKEN/D_RESULTS auto-nest there) and KRAKEN2_CONFIDENCE
# set to cs. A symlink makes each cs dir's 01_host_removed point at the shared one.
#
# Usage (on Sherlock, after sourcing nothing -- it sources config.sh itself):
#   bash sweep_confidence_e2e.sh [maxconcurrent]       # default maxconcurrent=8
#   CS_LIST="0 0.1 0.2 0.3" bash sweep_confidence_e2e.sh 4
# For cHL set RUN_DIR (as for run_all.sh), e.g.:
#   RUN_DIR=/oak/.../Binkley_cHL/cHL_microbiome bash sweep_confidence_e2e.sh 8
#
# COST WARNING: this submits len(CS_LIST) independent Kraken2 arrays, each of which
# reloads the ~140 GB PlusPF DB per task. 11 confidences x N samples is a LOT of
# DB loads. Keep maxconcurrent modest, and/or trim CS_LIST.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

MAXC="${1:-8}"
# "0 to 1 by 0.1" -> 11 values. Override CS_LIST to subset.
CS_LIST="${CS_LIST:-0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0}"

BASE_RUN_DIR="$RUN_DIR"
HOST_SHARED="$D_HOST"                         # the single, shared step-01 output (absolute)
SWEEP_ROOT="$BASE_RUN_DIR/confidence_sweep_e2e"

echo "=============================================================="
echo " Confidence sweep (end-to-end 02->03->05)"
echo "   RUN_DIR (base) = $BASE_RUN_DIR"
echo "   D_HOST (shared)= $HOST_SHARED"
echo "   CS_LIST        = $CS_LIST"
echo "   maxconcurrent  = $MAXC   (per-cs Kraken2 array)"
echo "   sweep outputs  -> $SWEEP_ROOT/cs_<cs>/05_results"
echo "=============================================================="

# ---- Prerequisite: step 01 host removal must be done (it is cs-independent) ----
shopt -s nullglob
host_fastqs=("$HOST_SHARED"/*_microbe_R1.fastq)
shopt -u nullglob
if [ "${#host_fastqs[@]}" -eq 0 ]; then
    echo "ERROR: no host-removed FASTQs in $HOST_SHARED." >&2
    echo "       Step 01 is confidence-INDEPENDENT and must run ONCE before the sweep." >&2
    echo "       Run it first, then re-run this script:" >&2
    echo "         bash $SCRIPT_DIR/00_build_sample_sheet.sh" >&2
    echo "         N=\$(wc -l < \"$SAMPLE_LIST\")" >&2
    echo "         sbatch --chdir=\"$BASE_RUN_DIR\" --export=ALL,PIPELINE_DIR=\"$SCRIPT_DIR\" \\" >&2
    echo "                --array=1-\$N%$MAXC \"$SCRIPT_DIR/01_host_removal.sbatch\"" >&2
    exit 1
fi
echo ">>> step 01 host removal present: ${#host_fastqs[@]} samples in $HOST_SHARED"

# ---- Preflight once (env/DB/tools are identical across cs) ----
if [ "${SWEEP_SKIP_PREFLIGHT:-0}" != "1" ]; then
    echo ">>> [preflight] checking environment once..."
    bash "$SCRIPT_DIR/preflight.sh" || {
        echo "ERROR: preflight failed; fix the FAILs above or set SWEEP_SKIP_PREFLIGHT=1." >&2; exit 1; }
fi

# ---- Per-confidence chains ----
mkdir -p "$SWEEP_ROOT"
for cs in $CS_LIST; do
    csdir="$SWEEP_ROOT/cs_${cs}"
    mkdir -p "$csdir"
    # Share the single host removal: csdir/01_host_removed -> the real one.
    ln -sfn "$HOST_SHARED" "$csdir/01_host_removed"
    echo ""
    echo ">>> ===== confidence=$cs  ->  $csdir ====="
    # run_all.sh builds its own sample list from the (symlinked) host FASTQs, runs
    # preflight (skipped here), and submits 02->03->05 with SLURM dependencies.
    # --export=ALL inside run_all propagates RUN_DIR + KRAKEN2_CONFIDENCE to the jobs;
    # config.sh keeps both because each is `${VAR:-default}`.
    RUN_ALL_SKIP_PREFLIGHT=1 RUN_DIR="$csdir" KRAKEN2_CONFIDENCE="$cs" \
        bash "$SCRIPT_DIR/run_all.sh" "$MAXC"
done

echo ""
echo "=============================================================="
echo ">>> Submitted ${CS_LIST} confidence chains. Track: squeue -u $USER"
echo ">>> When all finish, tabulate with:"
echo "      bash $SCRIPT_DIR/summarize_confidence_sweep.sh"
echo ">>> Per-cs results: $SWEEP_ROOT/cs_<cs>/05_results/{high_confidence_taxa.csv,filter_funnel_summary.csv}"
echo "=============================================================="
