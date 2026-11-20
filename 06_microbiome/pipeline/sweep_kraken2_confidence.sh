#!/bin/bash
#SBATCH --job-name=mb_cs_sweep
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --partition=emoding
#SBATCH --output=logs/cs_sweep_%A.out
#SBATCH --error=logs/cs_sweep_%A.err
# =============================================================================
# sweep_kraken2_confidence.sh -- pick a Kraken2 --confidence cutoff empirically
# =============================================================================
# Runs Kraken2 on a handful of (ideally KNOWN EBV+) samples at several confidence
# values and tabulates the four numbers you need to choose a cutoff:
#   - EBV signal (positive control)  -- the cutoff MUST retain this
#   - unclassified %                 -- how fast you are losing reads
#   - Homo sapiens reads / %         -- host that confidence pushes to "unclassified"
#   - species richness (>=1, >=10 reads) -- where false-positive richness collapses
#
# Confidence only re-thresholds the SAME k-mer classification (fraction of a read's
# k-mers that must support the assigned lineage); higher CS -> fewer, more-trusted
# calls + more unclassified. It does NOT help host removal (borderline human goes to
# "unclassified", not "Homo sapiens").  See Kraken2 issue #265 and Sun et al. 2024
# (aBIOTECH, doi:10.1007/s42994-024-00178-0).
#
# Usage (interactive salloc or sbatch; needs the PlusPF DB in RAM):
#   bash sweep_kraken2_confidence.sh cHL01 cHL07 cHL11      # explicit (EBV+) samples
#   N_DIAG=4 bash sweep_kraken2_confidence.sh               # first 4 from SAMPLE_LIST
#   CS_LIST="0 0.05 0.1 0.15 0.2" bash sweep_kraken2_confidence.sh cHL01
# Reuse the SAME DB node: --memory-mapping keeps the DB in OS page cache so only the
# FIRST kraken2 invocation pays the ~140 GB load; the rest are fast.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
mkdir -p "$LOG_DIR" logs

THREADS="${SLURM_CPUS_PER_TASK:-8}"
CS_LIST="${CS_LIST:-0 0.1 0.2 0.3 0.4 0.5}"
read -r -a CS_VALUES <<< "$CS_LIST"     # array (avoids glob/word-split surprises)
FORCE="${FORCE:-false}"                  # FORCE=true ignores cached reports
# EBV positive control. Clade reads at the Lymphocryptovirus GENUS (taxid 10375)
# capture ALL EBV-lineage reads, including those an aggressive CS resolves only to
# genus; the species column (10376) is the stricter EBV-only count. Override if your
# DB taxonomy differs.
EBV_TAXID="${EBV_TAXID:-10375}"          # genus Lymphocryptovirus (robust retention signal)
EBV_SP_TAXID="${EBV_SP_TAXID:-10376}"    # species (EBV-specific; understates if reads lift to genus)
HUMAN_TAXID="${HUMAN_TAXID:-9606}"

# ---- pick samples ----
if [ "$#" -gt 0 ]; then
    samples=("$@")
else
    N_DIAG="${N_DIAG:-4}"
    mapfile -t samples < <(cut -f2 "$SAMPLE_LIST" | head -n "$N_DIAG")
fi

outdir="$RUN_DIR/confidence_sweep"
mkdir -p "$outdir/reports"
out="$outdir/confidence_sweep_summary.tsv"
printf "sample\tconfidence\ttotal_pairs\tunclassified\tpct_unclassified\thuman_reads\tpct_human_of_classified\tebv_genus_reads\tebv_species_reads\tn_species_ge1\tn_species_ge10\n" > "$out"
echo ">>> samples: ${samples[*]}"
echo ">>> confidence values: $CS_LIST"
echo ">>> writing -> $out"

for s in "${samples[@]}"; do
    r1="$D_HOST/${s}_microbe_R1.fastq"; r2="$D_HOST/${s}_microbe_R2.fastq"
    if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
        echo "[skip] $s : missing $r1 / $r2" >&2; continue
    fi
    for cs in "${CS_VALUES[@]}"; do
        report="$outdir/reports/${s}.cs${cs}.kreport"
        stamp="$report.params"
        # Cache key: rerun if anything that affects the report changed. Includes the
        # FASTQ size+mtime (not just the path) so regenerating $D_HOST/*_microbe_R*.fastq
        # in place (e.g. a step-01 rerun) invalidates the cache instead of silently
        # reusing a stale report. FORCE=true ignores the cache entirely.
        # finfo: portable "size:mtime" (GNU stat -c, then BSD stat -f, then NA).
        finfo() { stat -c '%s:%Y' "$1" 2>/dev/null || stat -f '%z:%m' "$1" 2>/dev/null || echo NA; }
        sig="$(printf 'db=%s\nsif=%s\nmin_hit_groups=%s\nr1=%s\nr1meta=%s\nr2=%s\nr2meta=%s\nconfidence=%s\n' \
            "$KRAKEN2_DB" "$KRAKEN2_SIF" "$KRAKEN2_MIN_HIT_GROUPS" \
            "$r1" "$(finfo "$r1")" "$r2" "$(finfo "$r2")" "$cs")"
        if [ "$FORCE" = "true" ] || [ ! -s "$report" ] || [ ! -s "$stamp" ] || ! cmp -s <(printf "%s" "$sig") "$stamp"; then
            echo ">>> kraken2 $s @ confidence=$cs"
            # Atomic write: a killed run must not leave a half-written .kreport that the
            # cache then trusts. Write to tmp, mv on success, stamp the params.
            tmp="$report.tmp.$$"
            "$CONTAINER_RUN" exec "$KRAKEN2_SIF" kraken2 \
                --db "$KRAKEN2_DB" \
                --threads "$THREADS" \
                --memory-mapping \
                --paired "$r1" "$r2" \
                --confidence "$cs" \
                --minimum-hit-groups "$KRAKEN2_MIN_HIT_GROUPS" \
                --report "$tmp" \
                --output /dev/null
            mv "$tmp" "$report"
            printf "%s" "$sig" > "$stamp"
        else
            echo ">>> [cache] $s @ confidence=$cs (report up to date)"
        fi

        # Parse the 6-column kreport. Match host/EBV/root by TAXID (robust to the
        # leading-space indentation in the name column). root (taxid 1) clade reads
        # = total classified; total pairs = classified + unclassified.
        awk -F'\t' -v s="$s" -v cs="$cs" -v ebv="$EBV_TAXID" -v ebvsp="$EBV_SP_TAXID" -v hum="$HUMAN_TAXID" '
            BEGIN{ OFS="\t"; uncl=0; root=0; human=0; ebvg=0; ebvs=0; sp1=0; sp10=0 }
            {
                # Guard: this script makes a plain 6-col report. A minimizer report
                # (8 cols) would shift rank/taxid to $6/$7 and silently corrupt metrics.
                if (NF != 6) { printf "ERROR: expected 6-column kreport, got %d fields in %s\n", NF, FILENAME > "/dev/stderr"; exit 2 }
                clade=$2; rank=$4; taxid=$5;
                if (rank=="U")        uncl=clade;
                if (taxid==1)         root=clade;
                if (taxid==hum)       human=clade;
                if (taxid==ebv)       ebvg=clade;
                if (taxid==ebvsp)     ebvs=clade;
                if (rank=="S") { if (clade>=1) sp1++; if (clade>=10) sp10++ }
            }
            END{
                total=root+uncl;
                pu=(total>0)?100*uncl/total:0;
                ph=(root>0)?100*human/root:0;
                printf "%s\t%s\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%d\t%d\t%d\n", s, cs, total, uncl, pu, human, ph, ebvg, ebvs, sp1, sp10
            }' "$report" | tee -a "$out"
    done
done

echo "=================================================="
echo ">>> DONE. Summary -> $out"
echo ">>> Choose the HIGHEST confidence that still robustly retains ebv_reads in"
echo "    known-EBV+ samples; watch pct_unclassified and n_species_ge10 for the"
echo "    sensitivity cost. (host removal is unaffected by --confidence.)"
