#!/bin/bash
#
#SBATCH --job-name=cHL_trust4
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${MAIL_USER}
#SBATCH --partition=emoding,normal
#SBATCH --output=${LAB_ROOT}/sequencing/pipeline/runs/rnaseq/Binkley_cHL/TRUST4/logs/trust4_%j.out
#SBATCH --error=${LAB_ROOT}/sequencing/pipeline/runs/rnaseq/Binkley_cHL/TRUST4/logs/trust4_%j.err
# =============================================================================
# v3_K1_trust4_cHL_fastp.sh — TRUST4 on the fastp-trimmed cHL FASTQs (author 2026-11-03: "rerun trust4 on cHL_fastp under
# Binkley_cHL, refer to .../Binkley_NLPHL/TRUST4").
#
# Why: Fig S4a (TCR/BCR entropy) still comes from an older run on raw FASTQs (analysis/ajayss/tcellreceptor_analysis/cHL,
# TRUST4 v1.0.8): no fastp, the 3 Pico v2 bases still on R2, adapters not cut. cHL_fastp has all of these fixed
# (fastp --trim_front2 3 --detect_adapter_for_pe --correction).
#
# Same as Binkley_NLPHL/TRUST4/TRUST4.sh: TRUST4 install, references, -t 4, 4 samples at a time on 16 CPUs, 64 GB.
# Different from it:
#   - samples come from cHL_fastp/cHL_samples.txt (206); no batches and no merged technical replicates in this folder
#   - time limit 2 days, the most the normal partition allows (the NLPHL run took 23.4 h for about 186 samples,
#     close to its 1-day limit; at that rate 206 samples need about 26 h)
#   - partition emoding,normal
# Restart: submit the script again; samples that already have a report are skipped.
# Run on Sherlock (make the logs folder first; Slurm cannot create it):
#   mkdir -p ${LAB_ROOT}/sequencing/pipeline/runs/rnaseq/Binkley_cHL/TRUST4/logs
#   sbatch ${LAB_ROOT}/sequencing/pipeline/runs/rnaseq/Binkley_cHL/TRUST4/TRUST4_cHL_2026-11-03.sh
# =============================================================================
set -euo pipefail

# Step 1. Folders and files
export cohort_dir="${LAB_ROOT}/sequencing/pipeline/runs/rnaseq/Binkley_cHL"
export fastq_dir="$cohort_dir/cHL_fastp"                         # R1_<sample>.trimmed.fastq.gz, R2_<sample>.trimmed.fastq.gz
export sample_manifest="$fastq_dir/cHL_samples.txt"              # one sample per line (206)
export run_dir="$cohort_dir/TRUST4"
export output_dir="$run_dir/trust4_output"                       # one folder per sample
export logs_dir="$run_dir/trust4_logs"                           # one TRUST4 log per sample
export trust4_dir="${PROJECT_ROOT}/tools/TRUST4"    # v1.1.10-r630, same as the NLPHL run
export bcrtcr_ref="$trust4_dir/hg38_bcrtcr.fa"
export imgt_ref="$trust4_dir/human_IMGT+C.fa"
mkdir -p "$output_dir" "$logs_dir"

# Step 2. Tools (same modules as the NLPHL run)
module --force purge
module load biology
module load samtools
export PATH="$trust4_dir:$PATH"
for tool in run-trust4 samtools parallel; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: '$tool' not found on PATH" >&2
        exit 1
    fi
done
for reference in "$bcrtcr_ref" "$imgt_ref"; do
    if [ ! -f "$reference" ]; then
        echo "ERROR: TRUST4 reference $reference not found" >&2
        exit 1
    fi
done
echo ">>> TRUST4: $(which run-trust4)"
echo ">>> Output: $output_dir"

# Step 3. Build the sample table (sample, R1, R2). A missing FASTQ stops the run before anything starts.
sample_list="$logs_dir/sample_list.tsv"
> "$sample_list"
while read -r sample_name; do
    if [ -z "$sample_name" ]; then
        continue
    fi
    r1="$fastq_dir/R1_${sample_name}.trimmed.fastq.gz"
    r2="$fastq_dir/R2_${sample_name}.trimmed.fastq.gz"
    if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
        echo "ERROR: FASTQ missing for $sample_name" >&2
        exit 1
    fi
    printf "%s\t%s\t%s\n" "$sample_name" "$r1" "$r2" >> "$sample_list"
done < "$sample_manifest"
n_samples=$(wc -l < "$sample_list")
n_unique=$(cut -f1 "$sample_list" | sort -u | wc -l)
echo ">>> Samples to process: $n_samples (unique names: $n_unique)"
if [ "$n_samples" -ne 206 ] || [ "$n_unique" -ne 206 ]; then
    echo "ERROR: expected 206 unique samples" >&2
    exit 1
fi

# Step 4. TRUST4 for one sample. The report is the last file written, so it marks a finished sample.
run_one_trust4() {
    local sample_name=$1
    local r1=$2
    local r2=$3
    local sample_dir="$output_dir/$sample_name"
    local prefix="TRUST_${sample_name}"
    mkdir -p "$sample_dir"
    if [ -f "$sample_dir/${prefix}_report.tsv" ]; then
        echo "    [SKIP] $sample_name - already complete"
        return 0
    fi
    run-trust4 \
        -1 "$r1" -2 "$r2" \
        -f "$bcrtcr_ref" \
        --ref "$imgt_ref" \
        -o "$sample_dir/$prefix" \
        -t 4 \
        > "$logs_dir/${prefix}_trust4.log" 2>&1
}
export -f run_one_trust4

# Step 5. Run 4 samples at a time with 4 threads each (16 CPUs, as NLPHL).
#         A failed sample does not stop the others; parallel's exit status is the number of failed samples.
echo ">>> Running TRUST4 (4 jobs x 4 threads) from $(date)"
failed_samples=0
parallel --jobs 4 --colsep '\t' run_one_trust4 {1} {2} {3} :::: "$sample_list" || failed_samples=$?

# Step 6. Summary
n_done=$(find "$output_dir" -name "*_report.tsv" | wc -l)
echo "=========================================="
echo "  TRUST4 finished: $(date)"
echo "  Samples completed: $n_done / $n_samples"
echo "  parallel exit status (failed samples): $failed_samples"
echo "  Output: $output_dir/<sample>/TRUST_<sample>_*"
echo "  Logs:   $logs_dir/"
echo "=========================================="
