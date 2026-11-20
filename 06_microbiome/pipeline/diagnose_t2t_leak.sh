#!/bin/bash
# =============================================================================
# diagnose_t2t_leak.sh -- quantify T2T-alignable residual human in the microbial pool
# =============================================================================
# Step-01 B4 intentionally uses bowtie2 `--very-sensitive-local --un-conc`, which
# KEEPS non-concordant (discordant / one-mate-human) pairs by design -- that residual
# human is left for Kraken2 + the step-05 host filter to remove. This diagnostic
# MEASURES that retained-human fraction under a STRICTER, hypothetical both-mates-
# unmapped (`-f 12`) policy, WITHOUT rerunning the pipeline: it re-aligns the existing
# $D_HOST/<sample>_microbe_R*.fastq to T2T-CHM13 and counts pairs that are NOT
# both-mates-unmapped. It is a probe, not a description of current B4 behavior.
#
# Usage (interactive node or salloc; needs the T2T bowtie2 index + samtools):
#     bash diagnose_t2t_leak.sh                 # first 5 samples from SAMPLE_LIST
#     bash diagnose_t2t_leak.sh cHL01 cHL02     # explicit samples
# Override N default:  N_DIAG=8 bash diagnose_t2t_leak.sh
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "$SCRIPT_DIR/config.sh"

module load biology 2>/dev/null || true
module load samtools 2>/dev/null || true
module load biology bowtie2 2>/dev/null || true
[ -n "${BBMAP_BIN:-}" ] && export PATH="$BBMAP_BIN:$PATH"
have() { command -v "$1" >/dev/null 2>&1; }
have bowtie2  || { echo "ERROR: bowtie2 not found" >&2; exit 1; }
have samtools || { echo "ERROR: samtools not found" >&2; exit 1; }
[ -n "${T2T_BT2_INDEX:-}" ] || { echo "ERROR: T2T_BT2_INDEX unset in config.sh" >&2; exit 1; }
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---- pick samples ----
if [ "$#" -gt 0 ]; then
    samples=("$@")
else
    N_DIAG="${N_DIAG:-5}"
    mapfile -t samples < <(cut -f2 "$SAMPLE_LIST" | head -n "$N_DIAG")
fi

tmp=$(mktemp -d "${TMPDIR:-/tmp}/t2tleak.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

out="$RUN_DIR/t2t_leak_diagnostic.tsv"
printf "sample\tinput_pairs\tboth_unmapped_pairs\tt2t_aligned_pairs\tpct_t2t_aligned\n" > "$out"
echo "Writing -> $out"
echo "=================================================="

tot_in=0; tot_aln=0
for s in "${samples[@]}"; do
    r1="$D_HOST/${s}_microbe_R1.fastq"; r2="$D_HOST/${s}_microbe_R2.fastq"
    if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
        echo "[skip] $s : missing $r1 / $r2" >&2; continue
    fi
    in_pairs=$(($(wc -l < "$r1") / 4))

    # Re-align the current "microbial" survivors to T2T; count both-mates-unmapped.
    # Pairs that are NOT both-mates-unmapped = retained human a strict -f 12 policy
    # would drop (the current B4 --un-conc keeps them by design).
    bowtie2 -x "$T2T_BT2_INDEX" -1 "$r1" -2 "$r2" \
        --very-sensitive-local -p "$THREADS" 2>>"$tmp/${s}.bt2.log" \
      | samtools view -@ "$THREADS" -b -f 12 - \
      | samtools collate -@ "$THREADS" -O -u - "$tmp/coll_${s}" \
      | samtools fastq -@ "$THREADS" -n -1 "$tmp/u1.fq" -2 "$tmp/u2.fq" \
            -0 /dev/null -s /dev/null - 2>>"$tmp/${s}.fastq.log"
    both_unmapped=$(($(wc -l < "$tmp/u1.fq") / 4))
    t2t_aligned=$(( in_pairs - both_unmapped ))
    pct=$(awk -v a="$t2t_aligned" -v b="$in_pairs" 'BEGIN{ printf (b>0)? "%.2f":"NA", (b>0)?100*a/b:0 }')

    printf "%s\t%d\t%d\t%d\t%s\n" "$s" "$in_pairs" "$both_unmapped" "$t2t_aligned" "$pct" | tee -a "$out"
    tot_in=$(( tot_in + in_pairs )); tot_aln=$(( tot_aln + t2t_aligned ))
    rm -f "$tmp/u1.fq" "$tmp/u2.fq"
done

echo "=================================================="
awk -v a="$tot_aln" -v b="$tot_in" 'BEGIN{
    printf "TOTAL: input=%d  t2t_aligned=%d  (%.2f%% of the current microbial pool is T2T-alignable residual human that a strict both-mates-unmapped policy would drop)\n", b, a, (b>0)?100*a/b:0 }'
echo "NOTE: a probe of the CURRENT B4 --un-conc output, not current pipeline behavior."
echo "      The STAR --outFilterMultimapNmax effect needs a small step-01 rerun to measure."
