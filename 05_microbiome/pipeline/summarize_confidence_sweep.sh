#!/bin/bash
# =============================================================================
# summarize_confidence_sweep.sh -- tabulate the end-to-end confidence sweep
# =============================================================================
# Run AFTER sweep_confidence_e2e.sh chains have all finished. Reads each
# cs_<cs>/05_results/ and emits one row per confidence so you can pick a value:
#   tested_species         funnel stage-0 distinct species (input breadth)
#   hc_species / hc_calls  final HIGH-CONFIDENCE distinct species / sample-taxon calls
#   ebv_hc_rows            high-confidence rows matching EBV (Lymphocryptovirus) = positive control
#   median_pct_unclassified  median Kraken2 unclassified % (sensitivity cost as cs rises)
#
# Usage:  bash summarize_confidence_sweep.sh [SWEEP_ROOT]
#   (SWEEP_ROOT defaults to $RUN_DIR/confidence_sweep_e2e; set RUN_DIR for cHL.)
# -----------------------------------------------------------------------------
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

SWEEP_ROOT="${1:-$RUN_DIR/confidence_sweep_e2e}"
out="$SWEEP_ROOT/confidence_sweep_summary.tsv"
[ -d "$SWEEP_ROOT" ] || { echo "ERROR: no sweep dir $SWEEP_ROOT" >&2; exit 1; }

printf "confidence\ttested_species\thc_species\thc_calls\tebv_hc_rows\tmedian_pct_unclassified\tstatus\n" > "$out"
echo "Writing -> $out"
echo "=================================================="

shopt -s nullglob
csdirs=("$SWEEP_ROOT"/cs_*/)
[ "${#csdirs[@]}" -gt 0 ] || { echo "ERROR: no cs_* subdirs in $SWEEP_ROOT" >&2; exit 1; }

for csdir in "${csdirs[@]}"; do
    cs=$(basename "$csdir"); cs="${cs#cs_}"
    res="${csdir%/}/05_results"
    funnel="$res/filter_funnel_summary.csv"
    hc="$res/high_confidence_taxa.csv"
    qc="$res/QC_read_accounting.csv"

    if [ ! -s "$funnel" ]; then
        printf "%s\tNA\tNA\tNA\tNA\tNA\tincomplete\n" "$cs" | tee -a "$out"
        continue
    fi
    # filter_funnel_summary.csv: stage,calls,distinct_species,calls_dropped_here
    # (no stage name contains a comma). NR==2 = stage 0 (tested); last row = HIGH-CONFIDENCE.
    read -r tested hc_species hc_calls < <(
        awk -F, 'NR==2{t=$3} END{print t, $3, $2}' "$funnel")
    [ -n "${tested:-}" ] || { printf "%s\tNA\tNA\tNA\tNA\tNA\tincomplete\n" "$cs" | tee -a "$out"; continue; }

    ebv=0
    [ -s "$hc" ] && ebv=$(grep -Eic 'lymphocryptovirus|gammaherpesvirus 4' "$hc" || true)

    med="NA"
    if [ -s "$qc" ]; then
        med=$(awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="pct_unclassified") c=i; next}
                       c && $c!="" && $c!="NA" {print $c}' "$qc" \
              | sort -n \
              | awk '{a[n++]=$1} END{ if(n==0){print "NA"; exit}
                       if(n%2) printf "%.2f", a[int(n/2)];
                       else    printf "%.2f", (a[n/2-1]+a[n/2])/2 }')
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\tok\n" "$cs" "$tested" "$hc_species" "$hc_calls" "$ebv" "$med" | tee -a "$out"
done

echo "=================================================="
echo ">>> Done. Pick the HIGHEST confidence that keeps ebv_hc_rows > 0 (EBV is the"
echo "    positive control) while hc_species/hc_calls stay sensible and"
echo "    median_pct_unclassified hasn't exploded. Summary -> $out"
