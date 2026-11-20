#!/bin/bash
# =============================================================================
# config.sh  --  Central configuration for the tumor-microbiome pipeline
# =============================================================================
# Every step script sources this file, so all paths/parameters live in ONE place.
# Edit the values below to match your environment, then submit the steps.
#
# Convention note: this mirrors the lab's existing RNA-seq jobs
# (run_sf_array.sh, TRUST4.sh): SLURM array jobs, partition=emoding,
# mail-user=$USER, conda envs under brian/conda_envs, tools under brian/tools,
# and a tab-separated sample_list.tsv ( batch <TAB> sample <TAB> r1 <TAB> r2 ).
# =============================================================================

# Resolve the directory this pipeline lives in (so scripts can find siblings).
export PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# ---------------------------------------------------------------------------
# 1. COHORT INPUT / OUTPUT
# ---------------------------------------------------------------------------
# --- COHORT SELECTION --------------------------------------------------------
# Override these in the ENVIRONMENT to run a new cohort from a shared copy of this
# pipeline WITHOUT editing this file, e.g.:
#   TRIM_DIR=/oak/.../Binkley_cHL/cHL_fastp \
#   RUN_DIR=/oak/stanford/groups/emoding/sequencing/pipeline/runs/rnaseq/Binkley_cHL/cHL_microbiome \
#   SAMPLE_MANIFEST= \
#   bash run_all.sh 8
# cHL is a special case: its real run directory is cHL_microbiome, not the
# default $COHORT_DIR/microbiome, so set RUN_DIR explicitly as shown above.
# Each value falls back to the NLPHL default when unset, so the original NLPHL run
# still works with no environment variables. run_all.sh submits with
# `sbatch --export=ALL`, so exported overrides reach every step at runtime.
export COHORT_DIR="${COHORT_DIR:-/oak/stanford/groups/emoding/sequencing/pipeline/runs/rnaseq/Binkley_NLPHL}"

# Demuxed + fastp-trimmed paired FASTQs (the same inputs used for STAR-RSEM / salmon).
export TRIM_DIR="${TRIM_DIR:-$COHORT_DIR/NovaSeq_trimmed}"

# All pipeline outputs are written under here.
export RUN_DIR="${RUN_DIR:-$COHORT_DIR/microbiome}"
# LOG_DIR / SAMPLE_LIST are derived STRICTLY from RUN_DIR (unconditional `=`, NOT
# `:-`), exactly like the D_* dirs below. With `:-`, a value left exported in your
# shell by an earlier `source config.sh` for another cohort would leak into this
# run -- which is how a cHL run picked up the NLPHL sample_list.tsv. SAMPLE_LIST is
# an internal build artifact; to run a custom sample SOURCE, set SAMPLE_MANIFEST
# (or hand-edit the generated sample_list.tsv), never point SAMPLE_LIST elsewhere.
export LOG_DIR="$RUN_DIR/logs"
export SAMPLE_LIST="$LOG_DIR/sample_list.tsv"

# Authoritative sample list. If set AND present, the sample sheet is built from
# this manifest (one sample name per line; column 1) instead of globbing the
# directory -- this matches the expression pipeline (bs_STAR_align.sh reads the
# same file) and avoids stray/duplicate FASTQs in subfolders. Empty => glob.
# The `-` (not `:-`) lets an explicit `SAMPLE_MANIFEST=` force globbing; an unset
# value falls back to the NLPHL manifest. A manifest path that does not exist also
# falls back to globbing (see 00_build_sample_sheet.sh), so a new cohort without a
# manifest works either way.
export SAMPLE_MANIFEST="${SAMPLE_MANIFEST-$TRIM_DIR/LP_samples.txt}"

# Per-step output sub-directories.
export D_HOST="$RUN_DIR/01_host_removed"      # candidate-microbial FASTQs + STAR logs
export D_KRAKEN="$RUN_DIR/02_kraken2"         # .kraken / .kreport
export D_BRACKEN="$RUN_DIR/03_bracken"        # .bracken + merged matrices/biom
export D_RESULTS="$RUN_DIR/05_results"        # filtered Kraken2/Bracken tables
export D_FIGS="$RUN_DIR/06_figures"           # plots / PDFs

# FASTQ naming. Trimmed files are R1_<sample>.trimmed.fastq.gz (lab convention).
# The sample-sheet builder also tries a few common alternates automatically.
export R1_TAG="${R1_TAG:-R1_}"
export R1_SUFFIX="${R1_SUFFIX:-.trimmed.fastq.gz}"

# ---------------------------------------------------------------------------
# 2. COMPUTE (SLURM)
# ---------------------------------------------------------------------------
export PARTITION="emoding"
export MAIL_USER="${MAIL_USER:-$USER}"

# ---------------------------------------------------------------------------
# 3. HOST / HUMAN REMOVAL REFERENCES   (step 01)
# ---------------------------------------------------------------------------
# STAR index used for host depletion (GRCh38 / GENCODE v47 — the SAME index your
# bs_STAR_align.sh expression pipeline uses). Built WITH its annotation (sjdb
# baked in), so a GTF is NOT required at mapping time; leave STAR_GTF empty.
export STAR_GENOME_DIR="/oak/stanford/groups/emoding/scripts/rnaseq_scripts/index/STARIndex/gencode47_GRCh38"
export STAR_GTF=""

# STAR binary: use the locally-installed STAR 2.7.11b that your expression pipeline
# uses (and that built the index) — NOT the `star` module (2.7.10b), which would
# shadow it and refuse to load the genome.
export STAR_BIN_DIR="/oak/stanford/groups/emoding/analysis/brian/tools/bin"
export STAR_TWOPASS="false"      # 1-pass. (2-pass added ~7x runtime with NO host-capture gain on the cHL01 test.)
export STAR_FFPE_RELAXED="true"  # relaxed match/score thresholds -> capture more short FFPE host reads
export STAR_SEED_LMAX=""         # empty = STAR default (50, fast). Set 25 for max sensitivity (much SLOWER).
export STAR_MULTIMAP_NMAX="50"   # output multimap cap (STAR default 10). 50 captures reads mapping to 11-50
                                 # loci as host (instead of leaking as "too many loci"), WITHOUT touching
                                 # winAnchor: 50 == the default winAnchor, so the SEARCH is unchanged and fast
                                 # (no exhaustive anchor scan). Do NOT raise above 50 without also raising
                                 # STAR_WINANCHOR_NMAX -- that combo (1000/1000) backfired on cHL01: it killed
                                 # "too many loci" but exploded "unmapped: other" (~all human rRNA), GREW the
                                 # leaked pool +27%, and ran ~7x slower. rDNA/rRNA reads don't map to GRCh38 at
                                 # all (assembly gap), so no cap recovers them -- bbduk-rRNA/T2T/Kraken2 do.
export STAR_WINANCHOR_NMAX=""    # empty = STAR default (50). Leave default: it equals STAR_MULTIMAP_NMAX (50),
                                 # so the manual's winAnchor>=cap invariant holds with no search-speed cost.

# BBMap tools (bbduk.sh, clumpify.sh) for step-01 rRNA / low-complexity / dedup
# cleanup. `module load biology bbmap` does NOT expose them on this cluster, so
# install once and point BBMAP_BIN at the env's bin/ (empty = skip those steps):
#   conda create -p /oak/stanford/groups/emoding/analysis/brian/conda_envs/bbmap -c bioconda bbmap
#   then point BBMAP_BIN at the env's bin/ (it also contains java, so PATH alone works).
export BBMAP_BIN="/oak/stanford/groups/emoding/analysis/brian/conda_envs/bbmap/bin"

# Second-pass cleanup (recommended). Each is OPTIONAL and auto-skipped if the
# tool or reference is missing (the step logs what it actually ran).
#   - Human rRNA FASTA (e.g. SILVA human entries / NR_146144 etc.) for bbduk.
#   - Optional second host reference: T2T-CHM13 bowtie2 index prefix, which
#     traps residual human reads that GRCh38 misses.
# Human-only rRNA reference for the step-01 rRNA scrub. Build it once with
# build_human_rrna.sh (writes to this path). If absent, the rRNA sub-step skips.
# Refs are cohort-INDEPENDENT: keep them in a shared REFS_DIR and reuse for every
# cohort (build once, point every run here). Override REFS_DIR or either path in
# the environment if your install lives elsewhere.
export REFS_DIR="${REFS_DIR:-/oak/stanford/groups/emoding/analysis/brian/tools/microbiome/refs}"
export HUMAN_RRNA_FASTA="${HUMAN_RRNA_FASTA:-$REFS_DIR/human_rRNA.fa}"
export T2T_BT2_INDEX="${T2T_BT2_INDEX:-$REFS_DIR/chm13v2.0}"   # build once with build_t2t_index.sh (B4 skips if absent)
export LOWCOMPLEXITY_ENTROPY="0.2"    # bbduk entropy filter (0 disables). 0.2 = lenient: keeps more reads
                                      # (incl. AT-rich microbes); step-05 distinct-minimizer filter backstops FPs.
export DEDUP="false"                   # clumpify exact-duplicate (PCR) removal. OFF: the step-05
                                       # distinct-minimizer / minimizer-coverage filters are PCR-dup-robust
                                       # (duplicate reads share minimizers), so they cover the false-positive
                                       # role dedup played. NOTE: with dedup off, read counts / RPM / microbial
                                       # load are duplicate-INFLATED -- abundance is now "raw", not dedup'd.
                                       # (optical dedup also NOT enabled; add optical=t dupedist=N in step 01.)

# ---------------------------------------------------------------------------
# 4. KRAKEN2 + BRACKEN   (steps 02, 03)
# ---------------------------------------------------------------------------
export KRAKEN2_SIF="/oak/stanford/groups/emoding/analysis/brian/tools/microbiome/containers/kraken2_bracken_2.17.1_3.1.sif"
export BRACKEN_SIF="$KRAKEN2_SIF"
# DB MUST contain the human genome (so host reads are trapped as Homo sapiens,
# not misassigned to microbes) and must be a complete build (taxo.k2d present).
export KRAKEN2_DB="/oak/stanford/groups/emoding/analysis/brian/tools/microbiome/db/k2_pluspf_20260226"

# Container runtime for the Kraken2/Bracken .sif images (steps 02 & 03). Detected
# at runtime -- apptainer preferred, singularity as a fallback -- so a node that
# only provides one of them still works. Steps must use "$CONTAINER_RUN", never a
# hardcoded `apptainer`, or preflight could pass on singularity while the job dies.
if command -v apptainer >/dev/null 2>&1;   then export CONTAINER_RUN="apptainer"
elif command -v singularity >/dev/null 2>&1; then export CONTAINER_RUN="singularity"
else export CONTAINER_RUN="apptainer"; fi   # fall back; preflight flags a genuine absence

# Kraken2 calling: use Kraken2 defaults, then control false positives downstream
# with Bracken species abundance, distinct-minimizer, prevalence, host, and
# contaminant-blocklist filters.
export KRAKEN2_CONFIDENCE="${KRAKEN2_CONFIDENCE:-0.1}"  # env-overridable (sweep sets it per run); default 0.1
                                                       # (0.1 = literature floor for short-read FFPE vs a human-
                                                       #  containing DB; sweep 0-0.3 later to tune. See review.)
export KRAKEN2_MIN_HIT_GROUPS="2"     # Kraken2 default minimum hit groups
export REPORT_MINIMIZER_DATA="true"   # adds distinct-minimizer cols (KrakenUniq-style)

export READ_LEN="100"                 # FFPE median read length ~110; nearest Bracken preset
export BRACKEN_THRESH="5"             # drop species below this many reads
export BRACKEN_LEVELS="S G"           # space-separated LIST, one Bracken run per level.
                                      # "S"=species (REQUIRED -- step 05 needs microbiome_species.biom),
                                      # "G"=genus summary. Valid ranks: D P C O F G S (NOT "SG").
                                      # NOTE: 03_bracken_merge labels every non-S level "genus", so do
                                      # not add F/O/etc. here without generalizing that suffix map.

# ---------------------------------------------------------------------------
# 5. DOWNSTREAM FILTERING / DECONTAMINATION   (step 05)
# ---------------------------------------------------------------------------
export CONTAMINANT_BLOCKLIST="$PIPELINE_DIR/contaminants_blocklist.txt"
export MIN_READS="5"                   # min Kraken2/Bracken reads for a taxon in a sample
export MIN_RPM="0.5"                   # min reads per million INPUT read pairs (depth-normalized:
                                       # denominator = STAR "Number of input reads", i.e. the full
                                       # pre-host library, NOT the post-host microbial pool)
export MIN_PREVALENCE_FRAC="0.01"      # drop singletons; ~>=2 samples in a ~190-sample cohort
export MIN_MINIMIZER_COVERAGE="0.01"   # Kraken2 distinct-minimizers / total k-mers floor
export MIN_DISTINCT_MINIMIZERS="5"     # min distinct minimizers to keep a taxon (KrakenUniq-style)

# OPTIONAL decontam (R package). Provide a 2-column TSV ( sample <TAB> type )
# where type is "sample" or "control"; if set, step 05 runs decontam against the
# negative controls. Leave empty to skip (a kitome blocklist is used instead).
export CONTROLS_LIST=""
