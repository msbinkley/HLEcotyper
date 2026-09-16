#!/bin/bash
#SBATCH --job-name=mb_host_removal
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --partition=emoding
#SBATCH --output=logs/01_host_%A_%a.out
#SBATCH --error=logs/01_host_%A_%a.err
# -----------------------------------------------------------------------------
# STEP 01 -- Host (human) read removal  ==> candidate-microbial FASTQs
# -----------------------------------------------------------------------------
# Submit AFTER building the sample sheet:
#     sbatch --array=1-<N>%8 01_host_removal.sbatch
#
# What it does (and WHY), for each sample:
#   A. STAR -> GRCh38 (splice-aware). Reads that map are human; the UNMAPPED
#      reads are the candidate microbial pool. STAR's Log.final.out gives the
#      total input-read count we use later for RPM normalization.
#   B. (recommended, auto-skipped if tools/refs absent) clean the unmapped pool:
#         - bbduk rRNA k-mer removal (host rRNA survives step A and otherwise
#           dominates Pico-kit libraries),
#         - bbduk low-complexity (entropy) filter -- the #1 source of bogus
#           taxonomic hits on FFPE data,
#         - clumpify exact-duplicate (PCR) removal (low-input FFPE -> many dups;
#           optical dedup is not enabled -- see DEDUP note in config.sh),
#         - optional bowtie2 vs T2T-CHM13 to trap residual human GRCh38 misses.
#   C. Emit <sample>_microbe_R1/R2.fastq and a per-sample read-accounting line.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

# ---- modules / tools (faithful to lab convention; optional tools guarded) ----
# NOTE: do NOT `module --force purge` here -- it breaks the Sherlock module env.
module load biology 2>/dev/null || true
module load samtools 2>/dev/null || true
module load biology bbmap 2>/dev/null || true
module load biology bowtie2 2>/dev/null || true

# Use the locally-installed STAR 2.7.11b (matches bs_STAR_align.sh and the index
# build), NOT the older `star` module which would shadow it.
export PATH="$STAR_BIN_DIR:$PATH"
# Optional BBMap install (bbduk.sh, clumpify.sh) if the module didn't provide them.
[ -n "${BBMAP_BIN:-}" ] && export PATH="$BBMAP_BIN:$PATH"
have() { command -v "$1" >/dev/null 2>&1; }
have STAR || { echo "ERROR: STAR not found (check STAR_BIN_DIR=$STAR_BIN_DIR)" >&2; exit 1; }
echo ">>> STAR: $(STAR --version 2>&1) from $(command -v STAR)"

mkdir -p "$D_HOST" "$D_HOST/star" "$D_HOST/readcounts" "$LOG_DIR" logs

# ---- pull this array task's sample from the sheet ----
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
batch=$(echo "$LINE" | cut -f1)
sample=$(echo "$LINE" | cut -f2)
r1=$(echo "$LINE" | cut -f3)
r2=$(echo "$LINE" | cut -f4)
THREADS=${SLURM_CPUS_PER_TASK:-12}

echo "=================================================="
echo " Sample : $sample   (batch $batch, array ${SLURM_ARRAY_TASK_ID:-NA})"
echo "=================================================="

out_r1="$D_HOST/${sample}_microbe_R1.fastq"
out_r2="$D_HOST/${sample}_microbe_R2.fastq"
rc_tsv="$D_HOST/readcounts/${sample}.readcount.tsv"
done_flag="$D_HOST/readcounts/${sample}.done"
# Resume ONLY when the sample is fully complete: both FASTQs, the readcount row, AND
# the .done sentinel (written last). A job killed mid-finalize leaves no sentinel, so
# it is redone rather than skipped with truncated/accounting-missing data.
if [ -s "$out_r1" ] && [ -s "$out_r2" ] && [ -s "$rc_tsv" ] && [ -f "$done_flag" ]; then
    echo "[SKIP] $sample already host-removed (complete)."; exit 0
fi

work="$D_HOST/star/${sample}_"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/${sample}.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

# ---- ONE consolidated per-sample log: every step-01 tool writes here, each under a
# labelled banner, so you can tell which block came from which sub-step. ----
slog="$D_HOST/star/${sample}.step01.log"
banner() { printf '\n========================= %s =========================\n' "$*" >> "$slog"; }
{
  echo "############################################################"
  echo "# STEP 01 host removal  --  sample: $sample  (batch: $batch, array task: ${SLURM_ARRAY_TASK_ID:-NA})"
  echo "# date : $(date '+%Y-%m-%d %H:%M:%S %Z')   node: $(hostname)"
  echo "# params: STAR_TWOPASS=$STAR_TWOPASS  FFPE_RELAXED=$STAR_FFPE_RELAXED  MULTIMAP=${STAR_MULTIMAP_NMAX:-default(10)}  WINANCHOR=${STAR_WINANCHOR_NMAX:-default(50)}"
  echo "#         ENTROPY=$LOWCOMPLEXITY_ENTROPY  DEDUP=$DEDUP"
  echo "# sections: [A] STAR | [B1] bbduk rRNA | [B2] bbduk entropy | [B3] clumpify dedup | [B4] bowtie2 T2T | [C] read accounting"
  echo "############################################################"
} > "$slog"

# ---- STEP A: STAR de-hosting ----
# Mirrors bs_STAR_align.sh (2-pass, FFPE-relaxed) but adds --outReadsUnmapped
# Fastx (the unmapped mates are the microbial pool) and omits TranscriptomeSAM.
# The prebuilt index already contains the annotation, so a GTF is only added if
# STAR_GTF is set. For host removal, relaxed filters are ON by default so that
# more short FFPE host reads are mapped (= fewer human reads leak downstream).
echo ">>> [A] STAR -> GRCh38 (GENCODE v47)"
# --outSAMtype None: we never need the host alignment BAM, only the UNMAPPED
# reads (the microbial pool) and Log.final.out (read-count denominator for RPM).
# STAR still aligns and reports stats; it just skips writing/sorting the BAM,
# which is faster and avoids gigabytes of host reads we would only delete.
STAR_OPTS=()
[ -n "${STAR_GTF:-}" ] && [ -s "${STAR_GTF:-/nonexistent}" ] && STAR_OPTS+=(--sjdbGTFfile "$STAR_GTF")
[ "${STAR_TWOPASS:-false}" = "true" ] && STAR_OPTS+=(--twopassMode Basic)
if [ "${STAR_FFPE_RELAXED:-false}" = "true" ]; then
    STAR_OPTS+=(--outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 --outFilterMatchNmin 20)
fi
# Dense seeding (seedSearchStartLmax 25) was the main slowdown; default 50 is ~10x
# faster with negligible host-capture loss. Only added if STAR_SEED_LMAX is set.
[ -n "${STAR_SEED_LMAX:-}" ] && STAR_OPTS+=(--seedSearchStartLmax "$STAR_SEED_LMAX")
# Hyper-repetitive human reads that exceed STAR's multimap cap are emitted as
# UNMAPPED ("too many loci") and would leak into the microbial pool, where Kraken2
# then traps them as Homo sapiens. Raising the cap lets them map (to human repeats)
# so STAR removes them as host. Empty = STAR default (10).
[ -n "${STAR_MULTIMAP_NMAX:-}" ] && STAR_OPTS+=(--outFilterMultimapNmax "$STAR_MULTIMAP_NMAX")
# Anchor detection cap. MUST be >= outFilterMultimapNmax (STAR manual 4.1): otherwise
# seeds for hyper-repetitive human map to too many loci, get dropped as anchors, and
# the read leaks as "unmapped" before the output cap above ever applies.
[ -n "${STAR_WINANCHOR_NMAX:-}" ] && STAR_OPTS+=(--winAnchorMultimapNmax "$STAR_WINANCHOR_NMAX")
# Enforce that invariant at runtime so an env override (e.g. WINANCHOR=50 with
# MULTIMAP=1000) cannot silently make the multimap-cap raise a no-op.
if [ -n "${STAR_MULTIMAP_NMAX:-}" ] && [ -n "${STAR_WINANCHOR_NMAX:-}" ] \
   && [ "$STAR_WINANCHOR_NMAX" -lt "$STAR_MULTIMAP_NMAX" ]; then
    echo "ERROR: STAR_WINANCHOR_NMAX ($STAR_WINANCHOR_NMAX) < STAR_MULTIMAP_NMAX ($STAR_MULTIMAP_NMAX);" \
         "winAnchorMultimapNmax must be >= outFilterMultimapNmax (STAR manual 4.1)." >&2
    exit 1
fi
STAR --runThreadN "$THREADS" \
     --runMode alignReads \
     --genomeDir "$STAR_GENOME_DIR" \
     --readFilesIn "$r1" "$r2" \
     --readFilesCommand zcat \
     --outSAMtype None \
     --outReadsUnmapped Fastx \
     --outFileNamePrefix "$work" \
     ${STAR_OPTS[@]+"${STAR_OPTS[@]}"}

mate1="${work}Unmapped.out.mate1"
mate2="${work}Unmapped.out.mate2"
[ -s "$mate1" ] || { echo "ERROR: STAR produced no unmapped mate1 for $sample" >&2; exit 1; }

# total input reads (read PAIRS) for RPM normalization later
total_reads=$(grep "Number of input reads" "${work}Log.final.out" | awk '{print $NF}')
star_unmapped=$(($(wc -l < "$mate1") / 4))
echo "    input read pairs=$total_reads ; unmapped pairs=$star_unmapped"

# Fold STAR's Log.final.out into the consolidated log (the standalone file is removed
# in cleanup below, so this log is the single place to read the STAR summary).
banner "[A] STAR -> GRCh38  (Log.final.out)"
cat "${work}Log.final.out" >> "$slog" 2>/dev/null || echo "(Log.final.out missing)" >> "$slog"

cur1="$mate1"; cur2="$mate2"

# ---- STEP B: optional cleanup of the unmapped pool ----
after_rrna="$star_unmapped"; after_entropy="$star_unmapped"; after_dedup="$star_unmapped"; after_t2t="$star_unmapped"

# B1. rRNA removal (bbduk) -- removes pairs whose reads k-mer-match host rRNA
if have bbduk.sh && [ -s "$HUMAN_RRNA_FASTA" ]; then
    echo ">>> [B1] bbduk rRNA removal"
    banner "[B1] bbduk rRNA removal (k=31, ref=$HUMAN_RRNA_FASTA)"
    bbduk.sh -Xmx32g threads="$THREADS" \
        in1="$cur1" in2="$cur2" \
        out1="$tmp/no_rrna_1.fq" out2="$tmp/no_rrna_2.fq" \
        ref="$HUMAN_RRNA_FASTA" k=31 1>>"$slog" 2>&1
    cur1="$tmp/no_rrna_1.fq"; cur2="$tmp/no_rrna_2.fq"
    after_rrna=$(($(wc -l < "$cur1") / 4))
else
    echo ">>> [B1] SKIP rRNA removal (bbduk or HUMAN_RRNA_FASTA missing)"
    banner "[B1] SKIPPED rRNA removal (bbduk or HUMAN_RRNA_FASTA missing)"
fi

# B2. low-complexity filter (bbduk entropy)
if have bbduk.sh && [ "$LOWCOMPLEXITY_ENTROPY" != "0" ]; then
    echo ">>> [B2] bbduk entropy filter (>$LOWCOMPLEXITY_ENTROPY)"
    banner "[B2] bbduk low-complexity / entropy filter (entropy=$LOWCOMPLEXITY_ENTROPY, window=50, k=5)"
    bbduk.sh -Xmx32g threads="$THREADS" \
        in1="$cur1" in2="$cur2" \
        out1="$tmp/hc_1.fq" out2="$tmp/hc_2.fq" \
        entropy="$LOWCOMPLEXITY_ENTROPY" entropywindow=50 entropyk=5 1>>"$slog" 2>&1
    cur1="$tmp/hc_1.fq"; cur2="$tmp/hc_2.fq"
    after_entropy=$(($(wc -l < "$cur1") / 4))
else
    echo ">>> [B2] SKIP entropy filter"; after_entropy="$after_rrna"
    banner "[B2] SKIPPED entropy filter (LOWCOMPLEXITY_ENTROPY=0 or bbduk missing)"
fi

# B3. duplicate removal (clumpify)
if have clumpify.sh && [ "$DEDUP" = "true" ]; then
    echo ">>> [B3] clumpify dedup"
    banner "[B3] clumpify exact-duplicate (PCR) removal (dedupe=t)"
    clumpify.sh -Xmx32g threads="$THREADS" \
        in1="$cur1" in2="$cur2" \
        out1="$tmp/dd_1.fq" out2="$tmp/dd_2.fq" \
        dedupe=t 1>>"$slog" 2>&1
    cur1="$tmp/dd_1.fq"; cur2="$tmp/dd_2.fq"
    after_dedup=$(($(wc -l < "$cur1") / 4))
else
    echo ">>> [B3] SKIP dedup"; after_dedup="$after_entropy"
    banner "[B3] SKIPPED clumpify dedup (DEDUP=$DEDUP or clumpify missing)"
fi

# B4. optional second-host removal vs T2T-CHM13.
# bowtie2 --very-sensitive-local (local mode soft-clips partial-human reads, so
# adapter/FFPE-damaged ends still align). We then keep ONLY pairs where BOTH mates
# are unmapped (SAM flag 12 = read unmapped + mate unmapped). This DISCARDS:
#   - concordant human (clean pairs),
#   - DISCORDANT human (both ends hit T2T, odd geometry -- "aligned discordantly 1 time"),
#   - one-mate-human pairs.
# i.e. a pair survives only if NEITHER end matched the complete human (T2T) assembly.
# Cost: a genuine microbial read whose mate happens to match human is also dropped --
# rare for paired reads from one fragment; the bias here is toward maximal host removal.
# `samtools collate` regroups mates so `samtools fastq` emits intact pairs. Requires
# samtools (loaded at top); if absent, SKIP rather than silently keep human.
if [ -n "$T2T_BT2_INDEX" ] && have bowtie2 && have samtools \
   && { [ -s "${T2T_BT2_INDEX}.1.bt2" ] || [ -s "${T2T_BT2_INDEX}.1.bt2l" ]; }; then
    echo ">>> [B4] bowtie2 vs T2T-CHM13 (keep ONLY pairs where BOTH mates are unmapped)"
    banner "[B4] bowtie2 vs T2T-CHM13 (--very-sensitive-local; keep both-mates-unmapped via samtools -f 12)"
    bowtie2 -x "$T2T_BT2_INDEX" -1 "$cur1" -2 "$cur2" \
        --very-sensitive-local -p "$THREADS" 2>>"$slog" \
      | samtools view -@ "$THREADS" -b -f 12 - \
      | samtools collate -@ "$THREADS" -O -u - "$tmp/collate_t2t" \
      | samtools fastq -@ "$THREADS" -n -1 "$tmp/noT2T_1.fq" -2 "$tmp/noT2T_2.fq" \
            -0 /dev/null -s /dev/null - 2>>"$slog"
    cur1="$tmp/noT2T_1.fq"; cur2="$tmp/noT2T_2.fq"
    after_t2t=$(($(wc -l < "$cur1") / 4))
elif [ -n "$T2T_BT2_INDEX" ] && have bowtie2 \
     && { [ -s "${T2T_BT2_INDEX}.1.bt2" ] || [ -s "${T2T_BT2_INDEX}.1.bt2l" ]; }; then
    echo ">>> [B4] SKIP T2T removal (samtools missing -- needed for both-mates-unmapped extraction)"; after_t2t="$after_dedup"
    banner "[B4] SKIPPED T2T removal (samtools missing)"
else
    echo ">>> [B4] SKIP T2T removal (T2T_BT2_INDEX unset, index files absent, or bowtie2 missing)"; after_t2t="$after_dedup"
    banner "[B4] SKIPPED T2T removal (T2T_BT2_INDEX unset, index files absent, or bowtie2 missing)"
fi

# ---- STEP C: finalize (atomic) ----
# Stage to .partial, then rename within $D_HOST (same filesystem -> atomic) only after
# the readcount is staged too; drop the .done sentinel LAST. The skip-guard requires
# all of {R1, R2, readcount, .done}, so an interrupted finalize never looks complete.
mv -f "$cur1" "$out_r1.partial"
mv -f "$cur2" "$out_r2.partial"
final=$(($(wc -l < "$out_r1.partial") / 4))

{ printf "sample\tinput_pairs\tstar_unmapped\tafter_rrna\tafter_entropy\tafter_dedup\tafter_t2t\tfinal_microbial\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sample" "$total_reads" "$star_unmapped" "$after_rrna" "$after_entropy" "$after_dedup" "$after_t2t" "$final"
} > "$rc_tsv.partial"

rm -f "$done_flag"                       # clear any stale sentinel before committing
mv -f "$out_r1.partial" "$out_r1"
mv -f "$out_r2.partial" "$out_r2"
mv -f "$rc_tsv.partial" "$rc_tsv"
: > "$done_flag"                         # sentinel LAST = sample is fully complete

# Mirror the accounting into the consolidated log so one file tells the whole story.
banner "[C] read accounting (pairs surviving each stage)"
{
  printf "  input_pairs      = %s\n" "$total_reads"
  printf "  star_unmapped    = %s   (STAR removed %s as host)\n" "$star_unmapped" "$((total_reads - star_unmapped))"
  printf "  after_rrna       = %s   (bbduk rRNA removed %s)\n"   "$after_rrna"    "$((star_unmapped - after_rrna))"
  printf "  after_entropy    = %s   (entropy removed %s)\n"      "$after_entropy" "$((after_rrna - after_entropy))"
  printf "  after_dedup      = %s   (dedup removed %s)\n"        "$after_dedup"   "$((after_entropy - after_dedup))"
  printf "  after_t2t        = %s   (T2T removed %s)\n"          "$after_t2t"     "$((after_dedup - after_t2t))"
  printf "  final_microbial  = %s\n" "$final"
} >> "$slog"

# Discard STAR scratch. Log.final.out is now folded into the consolidated log
# (${sample}.step01.log), so the standalone copy is removed too -- that consolidated
# log is the single place to read every step-01 tool's output.
rm -rf "${work}Unmapped.out.mate"* "${work}_STARpass1" "${work}_STARtmp" \
       "${work}_STARgenome" "${work}SJ.out.tab" "${work}Log.out" \
       "${work}Log.progress.out" "${work}Aligned"* "${work}Log.final.out" \
       "${work}clean.log" 2>/dev/null || true

echo ">>> DONE $sample : final candidate-microbial pairs = $final"
echo ">>> full per-step log: $slog"
