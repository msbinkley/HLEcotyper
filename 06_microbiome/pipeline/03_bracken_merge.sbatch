#!/bin/bash
#SBATCH --job-name=mb_bracken_merge
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=emoding
#SBATCH --output=logs/03_bracken_%j.out
#SBATCH --error=logs/03_bracken_%j.err
# -----------------------------------------------------------------------------
# STEP 03 -- Bracken abundance re-estimation + cohort merge  (single job)
# -----------------------------------------------------------------------------
#     sbatch 03_bracken_merge.sbatch        # run after ALL step-02 array tasks
#
# Bracken redistributes Kraken2's reads to give proper species abundances.
# We then (a) merge per-sample reports into a BIOM matrix (for phyloseq/decontam)
# and (b) build a distinct-minimizer table from the Kraken2 reports so step 05
# can drop taxa supported by too few unique k-mers (KrakenUniq-style FP control).
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
mkdir -p "$D_BRACKEN" logs

# ---- A. Bracken per sample, per level ----
shopt -s nullglob
reports=("$D_KRAKEN"/*.kreport)
[ "${#reports[@]}" -gt 0 ] || { echo "ERROR: no .kreport files in $D_KRAKEN" >&2; exit 1; }

declare -A seen_reports=()
for rep in "${reports[@]}"; do
    seen_reports["$(basename "$rep" .kreport)"]=1
done
expected=0
missing_samples=()
while IFS=$'\t' read -r _batch sample _r1 _r2 || [ -n "${sample:-}" ]; do
    [ -z "${sample:-}" ] && continue
    expected=$((expected+1))
    if [ -z "${seen_reports[$sample]:-}" ]; then
        missing_samples+=("$sample")
    fi
done < "$SAMPLE_LIST"
if [ "${#missing_samples[@]}" -gt 0 ]; then
    found_expected=$((expected - ${#missing_samples[@]}))
    echo "ERROR: only $found_expected of $expected expected Kraken2 .kreport files found in $D_KRAKEN for $SAMPLE_LIST." >&2
    echo "       Missing samples: ${missing_samples[*]}" >&2
    exit 1
fi

bfail=0
for level in $BRACKEN_LEVELS; do
    suffix=$([ "$level" = "S" ] && echo "species" || echo "genus")
    echo ">>> Bracken level=$level ($suffix)"
    for rep in "${reports[@]}"; do
        base=$(basename "$rep" .kreport)
        out_b="$D_BRACKEN/${base}.${suffix}.bracken"
        out_k="$D_BRACKEN/${base}.${suffix}_bracken.kreport"
        [ -s "$out_k" ] && { echo "    [skip] $base ($suffix)"; continue; }
        # Write to temp, rename on success only. A failed Bracken must NOT leave a
        # partial *_bracken.kreport behind -- otherwise it would pass the completeness
        # gate below (which tests for the file's presence) and get merged.
        tb="$out_b.partial.$$"; tk="$out_k.partial.$$"
        if "$CONTAINER_RUN" exec "$BRACKEN_SIF" bracken \
                -d "$KRAKEN2_DB" -i "$rep" \
                -o "$tb" -w "$tk" \
                -r "$READ_LEN" -l "$level" -t "$BRACKEN_THRESH"; then
            mv -f "$tb" "$out_b"; mv -f "$tk" "$out_k"
        else
            echo "    WARNING: bracken failed for $base ($suffix)"; bfail=$((bfail+1))
            rm -f "$tb" "$tk"
        fi
    done
done
[ "$bfail" -gt 0 ] && echo "    NOTE: $bfail bracken sample-levels failed this run (see warnings above)."

# Fail-closed completeness gate: every expected sample must have a Bracken output for
# every requested level before we merge. A silently missing/failed sample would be
# merged as ABSENT and then mis-counted in step 05 as a true ZERO-CALL sample,
# deflating the prevalence denominator (n_total includes it via readcounts, but it
# contributes no calls). Set ALLOW_PARTIAL_BRACKEN=true to accept the gaps as genuine
# zero-call samples (e.g. a real low-biomass sample Bracken could not estimate).
ALLOW_PARTIAL_BRACKEN="${ALLOW_PARTIAL_BRACKEN:-false}"
bmissing=()
for level in $BRACKEN_LEVELS; do
    suffix=$([ "$level" = "S" ] && echo "species" || echo "genus")
    for rep in "${reports[@]}"; do
        base=$(basename "$rep" .kreport)
        [ -s "$D_BRACKEN/${base}.${suffix}_bracken.kreport" ] || bmissing+=("${base}.${suffix}")
    done
done
if [ "${#bmissing[@]}" -gt 0 ]; then
    echo "ERROR: ${#bmissing[@]} expected Bracken output(s) missing/empty; a failed sample must NOT be" >&2
    echo "       merged as a zero-call sample (it would bias the step-05 prevalence denominator)." >&2
    echo "       Missing: ${bmissing[*]}" >&2
    if [ "$ALLOW_PARTIAL_BRACKEN" = "true" ]; then
        echo "       ALLOW_PARTIAL_BRACKEN=true -> proceeding; treating those as genuine zero-call samples." >&2
    else
        echo "       Re-run the failed Bracken sample-levels, remove the sample(s) from $SAMPLE_LIST," >&2
        echo "       or set ALLOW_PARTIAL_BRACKEN=true to accept them as real zero-call. Aborting." >&2
        exit 1
    fi
fi

# ---- B. Merge into BIOM matrices ----
echo ">>> Merging Bracken reports into BIOM matrices (kraken-biom)"
module load python/3.9.0 2>/dev/null || true
pip install --user --quiet kraken-biom 2>/dev/null || true
KB="$HOME/.local/bin/kraken-biom"
command -v kraken-biom >/dev/null 2>&1 && KB="kraken-biom"

for level in $BRACKEN_LEVELS; do
    suffix=$([ "$level" = "S" ] && echo "species" || echo "genus")
    kfiles=("$D_BRACKEN"/*.${suffix}_bracken.kreport)
    if [ "${#kfiles[@]}" -gt 0 ] && [ -x "$KB" -o "$KB" = "kraken-biom" ]; then
        "$KB" "${kfiles[@]}" --fmt json -o "$D_BRACKEN/microbiome_${suffix}.biom" \
            && echo "    wrote microbiome_${suffix}.biom" \
            || echo "    WARNING: kraken-biom failed for $suffix"
    else
        echo "    WARNING: no ${suffix}_bracken.kreport files, or kraken-biom unavailable -- skipped $suffix BIOM"
    fi
done

# Step 05 hard-requires the species BIOM. Fail loudly HERE rather than letting
# the dependency chain "succeed" and then crashing inside R on a missing file.
if [[ " $BRACKEN_LEVELS " == *" S "* ]] && [ ! -s "$D_BRACKEN/microbiome_species.biom" ]; then
    echo "ERROR: species BIOM was not produced ($D_BRACKEN/microbiome_species.biom)." >&2
    echo "       Check the kraken-biom install (pip install --user kraken-biom) and the bracken output above." >&2
    exit 1
fi

# ---- C. Distinct-minimizer table (KrakenUniq-style FP control) ----
# From each 8-column minimizer report, keep species rows:
#   sample, taxid, name, taxon_reads, total_minimizers, distinct_minimizers
echo ">>> Building distinct-minimizer table"
mmz_out="$D_BRACKEN/minimizer_stats.tsv"
printf "sample\ttaxid\tname\ttaxon_reads\ttotal_minimizers\tdistinct_minimizers\n" > "$mmz_out"
for mr in "$D_KRAKEN"/*.kreport.minimizer; do
    [ -e "$mr" ] || { echo "    (no minimizer reports; skipping)"; break; }
    s=$(basename "$mr" .kreport.minimizer)
    # columns: 1%  2clade_reads  3taxon_reads  4tot_minz  5distinct_minz  6rank  7taxid  8name
    awk -F'\t' -v S="$s" 'BEGIN{OFS="\t"} $6=="S" {
        name=$8; gsub(/^ +/,"",name);
        print S,$7,name,$3,$4,$5 }' "$mr" >> "$mmz_out"
done

# ---- D. STAR total-read lookup (for RPM), mirroring the scratch logic ----
echo ">>> Building total_reads_lookup.tsv (from readcount.tsv input_pairs column)"
lookup="$D_BRACKEN/total_reads_lookup.tsv"
: > "$lookup"
for rc in "$D_HOST"/readcounts/*.readcount.tsv; do
    [ -e "$rc" ] || break
    awk -F'\t' 'NR>1 {print $1"\t"$2}' "$rc" >> "$lookup"   # sample, input_pairs (skip header)
done

# ---- E. Kraken2 classification QC: classified vs unclassified per sample ----
# Reads that fail the configured Kraken2 thresholds (or have no DB hit) are
# reported as "unclassified" (taxid 0). Root (taxid 1) clade-reads = classified.
echo ">>> Building kraken2_classification_QC.tsv"
kqc="$D_BRACKEN/kraken2_classification_QC.tsv"
printf "sample\tclassified_pairs\tunclassified_pairs\ttotal_pairs\tpct_unclassified\n" > "$kqc"
for rep in "${reports[@]}"; do
    base=$(basename "$rep" .kreport)
    u=$(awk -F'\t' '$5=="0"{print $2; exit}' "$rep"); u=${u:-0}   # unclassified (taxid 0)
    r=$(awk -F'\t' '$5=="1"{print $2; exit}' "$rep"); r=${r:-0}   # classified  (root, taxid 1)
    tot=$((u + r))
    pct=$(awk -v u="$u" -v t="$tot" 'BEGIN{ if(t>0) printf "%.2f",100*u/t; else printf "NA" }')
    printf "%s\t%s\t%s\t%s\t%s\n" "$base" "$r" "$u" "$tot" "$pct" >> "$kqc"
done

echo ">>> STEP 03 complete."
