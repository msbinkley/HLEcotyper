# =============================================================================
# v2_A_gene_universe.R — step A of REANALYSIS_PLAN_v2: the shared gene universe.
#
# Runs on Sherlock (R/4.2.0).
#
# What it does:
#   1. Builds the GENCODE v47 gene table.
#   2. Maps every scRNA dataset's features to stable v47 gene IDs, through the annotation each dataset was
#      quantified with (table below).
#   3. Defines the shared gene universe.
#   (The bulk TPM is no longer built here: since 2026-10-21 it is built by v2_A4_bulk_TPM_coding.R;
#    author: the bulk keeps ALL protein-coding genes, not only the universe.)
#
#   dataset   features in object           annotation used for mapping
#   lab       symbols (sc 3' + sn Flex1.0)  Cell Ranger GRCh38-2020-A (GENCODE v32)  [lab sc reference confirmed by overlap]
#   Aoki      symbols + feature meta `id`  ENSG from the object (Ensembl 84 era)
#   Stewart   ENSG                          ENSG
#   10x       symbols                       GRCh38-2020-A (confirmed by overlap)
#   S1        symbols (Flex probe v1.1)     GRCh38-2024-A (GENCODE v44)
#
# Gene universe = v47 protein-coding genes + IG/TR constant-region genes (IG_C_gene, TR_C_gene)  (no PAR_Y copies)
#                   [author 2026-10-14: the IG/TR constant genes encode protein but GENCODE types them separately,
#                    and both Flex panels probe 16 of the 20]
#               ∩ Flex probe set v1.0 genes ∩ Flex probe set v1.1 genes   (measured by both Flex chemistries)
#               ∩ genes present in every dataset after mapping           (lab, Aoki, Stewart, 10x, S1; but see the zero-fill note in Step 4)
#               ∩ genes in the new salmon ENSG TPM
#
# Outputs (in <W>/ref):
#   gene_table_v47.tsv, feature_map_<dataset>.tsv, gene_universe.tsv, harmonisation_report.tsv
#   (bulk_TPM_v47_universe_187.tsv + bulk_sample_check.tsv were written here until 2026-10-21; now v2_A4)
#
# Run: Rscript v2_A_gene_universe.R   (optional env var TPM_ENSG = path of the bulk ENSG TPM table)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
})

# ---- Settings: paths -------------------------------------------------------------------------------
W <- "${PROJECT_ROOT}/cHL_v2"            # project folder on Sherlock

# gene annotations (GTF): v47 = target annotation; v32 / v44 = the Cell Ranger references the symbol-only datasets used
GTF <- c(v47 = "${LAB_ROOT}/scripts/rnaseq_scripts/index/SalmonIndex/gencodeIndex47_k21_v1.11.4/gencode.v47.annotation.gtf",
         v32 = "${LAB_ROOT}/scripts/human_indices/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
         v44 = "${LAB_ROOT}/scripts/human_indices/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz")

# 10x Flex probe sets (v1.0 = lab snRNA; v1.1 = ex vivo S1)
PROBE <- c(flex_v1.0 = "/home/groups/emoding/cellranger_install/cellranger-7.0.1/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv",
           flex_v1.1 = "${LAB_ROOT}/sequencing/singlecell/2026-10-01_cHL_exvivo/cHL_exvivo_Pool1/outs/probe_set.csv")

# bulk ENSG TPM (used only for the universe's "present in bulk" filter): the 2026-10-14 `-l ISR --allowDovetail` run
TPM_ENSG <- Sys.getenv("TPM_ENSG", "${LAB_ROOT}/sequencing/pipeline/runs/rnaseq/Binkley_cHL/cHL_salmon_ISR_dovetail/cHL_salmon_ISR_dovetail_gene_ensg_TPM.tsv")

# ex vivo Flex object (S1 = its Sample1)
S1_OBJ <- "${PROJECT_ROOT}/scRNA/2026-10-04_cHL_ex_vivo/output/2026-10-04_chl_exvivo_final_seurat_obj_f6.rds"

OUT <- file.path(W, "ref")
dir.create(OUT, showWarnings = FALSE)

# stable(): drop the version from a gene ID, e.g. ENSG00000123.4 -> ENSG00000123 (a "_PAR_Y" ending is kept)
stable <- function(x) sub("\\.[0-9]+(_PAR_Y)?$", "\\1", x)


# ---- Step 1. Gene tables from the GTF files ---------------------------------------------------------
# read_gtf_genes(): read the "gene" lines of a GTF (gzipped or not) and return gene_id, gene_name, gene_type
read_gtf_genes <- function(gtf_file) {
  cat_cmd <- if (grepl("[.]gz$", gtf_file)) "gzip -dc" else "cat"
  # awk keeps column 9 (the attribute text) of every line whose feature type (column 3) is "gene"
  attribute_lines <- fread(cmd = sprintf("%s '%s' | awk -F'\\t' '$3==\"gene\"{print $9}'", cat_cmd, gtf_file),
                           sep = "\n", header = FALSE)$V1
  # get_attribute(): pull the value of one attribute, e.g. gene_name "CD3E" -> CD3E
  get_attribute <- function(key) sub(sprintf('.*%s "([^"]+)".*', key), "\\1", attribute_lines)
  data.table(gene_id = get_attribute("gene_id"), gene_name = get_attribute("gene_name"),
             gene_type = get_attribute("gene_type"))
}

gene_tables <- list()                                   # one table per annotation: v47, v32, v44 (in this order)
for (annotation in names(GTF)) {
  gene_tables[[annotation]] <- read_gtf_genes(GTF[[annotation]])
}

# the v47 gene table: versioned ID, stable ID, symbol, gene type, and whether it is a PAR_Y copy
v47 <- gene_tables$v47[, .(ensg_versioned = gene_id, ensg = stable(gene_id), symbol = gene_name, type = gene_type)]
v47[, par_y := grepl("_PAR_Y$", ensg)]
fwrite(v47, file.path(OUT, "gene_table_v47.tsv"), sep = "\t")
cat(sprintf("v47: %d genes (%d protein_coding, %d PAR_Y)\n", nrow(v47), sum(v47$type == "protein_coding"), sum(v47$par_y)))

# to47(): stable ID if it is a v47 gene (not a PAR_Y copy), otherwise NA
to47 <- function(ensg) {
  stable_id <- stable(ensg)
  ifelse(stable_id %in% v47$ensg[!v47$par_y], stable_id, NA_character_)
}

# symbol -> ENSG lookups for the Cell Ranger references.
# Only symbols that occur exactly once are kept; a symbol shared by several genes is ambiguous and dropped.
sym_lookup <- function(gene_table) {
  repeated_names <- gene_table$gene_name[duplicated(gene_table$gene_name)]
  unique_name <- !duplicated(gene_table$gene_name) & !(gene_table$gene_name %in% repeated_names)
  gene_table <- gene_table[unique_name]
  setNames(stable(gene_table$gene_id), gene_table$gene_name)
}
lookup_v32 <- sym_lookup(gene_tables$v32)      # GRCh38-2020-A (GENCODE v32)
lookup_v44 <- sym_lookup(gene_tables$v44)      # GRCh38-2024-A (GENCODE v44)


# ---- Step 2. Each dataset's features -----------------------------------------------------------------
# map_symbols(): symbol -> ENSG through a lookup. A symbol not found is tried again without a Seurat
# make.unique suffix (e.g. "TBCE.1" -> "TBCE").
map_symbols <- function(symbols, lookup) {
  ensg <- unname(lookup[symbols])
  not_found <- is.na(ensg)
  symbols_no_suffix <- sub("\\.[0-9]+$", "", symbols)
  ensg[not_found] <- unname(lookup[symbols_no_suffix[not_found]])
  ensg
}

# read each object only to get its feature (gene) names; drop the object straight after
features <- list()
lab <- readRDS(file.path(W, "input/lab_scRNA_seurat.rds"))
features$lab <- rownames(lab)
rm(lab)
tenx <- readRDS(file.path(W, "input/tenx_filtered.rds"))
features$tenx <- rownames(tenx)
rm(tenx)
stewart <- readRDS(file.path(W, "input/stewart_filtered.rds"))
features$stewart <- rownames(stewart)
rm(stewart)
aoki <- readRDS(file.path(W, "input/aoki_filtered.rds"))
aoki_feature_meta <- aoki[["originalexp"]][[]]          # Aoki feature metadata: column `id` holds the ENSG
features$aoki <- rownames(aoki)
rm(aoki)
s1 <- readRDS(S1_OBJ)
features$s1 <- rownames(s1)
rm(s1)
invisible(gc())

# Which reference fits the symbol-only datasets best (lab sc, 10x, S1)? Print the % of symbols found in each.
for (dataset in c("lab", "tenx", "s1")) {
  dataset_symbols <- features[[dataset]]
  pct_v32 <- 100 * mean(dataset_symbols %in% names(lookup_v32) | sub("\\.[0-9]+$", "", dataset_symbols) %in% names(lookup_v32))
  pct_v44 <- 100 * mean(dataset_symbols %in% names(lookup_v44))
  pct_v47 <- 100 * mean(dataset_symbols %in% v47$symbol)
  cat(sprintf("%-5s symbols found: v32 %.1f%%  v44 %.1f%%  v47 %.1f%%\n", dataset, pct_v32, pct_v44, pct_v47))
}


# ---- Step 3. Per-dataset feature -> v47 maps -------------------------------------------------------------
# ENSG of every feature, through the annotation each dataset was quantified with (see table in the header)
src_ensg <- list(
  lab     = map_symbols(features$lab,  lookup_v32),
  tenx    = map_symbols(features$tenx, lookup_v32),
  s1      = map_symbols(features$s1,   lookup_v44),
  aoki    = stable(aoki_feature_meta[features$aoki, "id"]),
  stewart = stable(features$stewart))

# one map table per dataset: dataset, feature, source ENSG, v47 ID (NA = not a v47 gene)
maps <- list()
for (dataset in names(src_ensg)) {
  maps[[dataset]] <- data.table(dataset = dataset, feature = features[[dataset]],
                                src_ensg = src_ensg[[dataset]], v47 = to47(src_ensg[[dataset]]))
}
for (dataset in names(maps)) {
  fwrite(maps[[dataset]], file.path(OUT, sprintf("feature_map_%s.tsv", dataset)), sep = "\t")
}


# ---- Step 4. Gene sets that make up the universe -------------------------------------------------------
# genes in a Flex probe set (only probes marked included), as v47 IDs
probe_genes <- function(probe_file) {
  probes <- fread(probe_file, skip = "gene_id")
  included_genes <- unique(probes[included %in% c(TRUE, "TRUE"), gene_id])
  to47(included_genes)
}
probe_sets <- list()
for (panel in names(PROBE)) {
  probe_sets[[panel]] <- probe_genes(PROBE[[panel]])
}

tpm_genes <- stable(fread(TPM_ENSG, select = 1)[[1]])           # genes in the bulk TPM table

# protein-coding + IG/TR constant-region genes, no PAR_Y copies
protein_coding <- v47[type %in% c("protein_coding", "IG_C_gene", "TR_C_gene") & !par_y, ensg]

## A gene absent from a whole-transcriptome object (lab, 10x, Stewart, S1 after f6 filtering) had zero counts there, which is a
## measurement, so those datasets are zero-filled rather than intersected. Aoki is the exception: its object carries the full
## GRCh38-1.2.0 (Ensembl 84) gene list, so a gene missing there was never annotated and is excluded.
sets <- list(protein_coding_IG_TR_C_v47 = protein_coding)
for (panel in names(probe_sets)) {
  panel_genes <- probe_sets[[panel]]
  sets[[panel]] <- panel_genes[!is.na(panel_genes)]
}
sets$aoki_annotation <- unique(na.omit(maps$aoki$v47))
sets$bulk_tpm <- tpm_genes


# ---- Step 5. The universe = genes in every one of the sets -------------------------------------------------
universe <- sets[[1]]
for (i in 2:length(sets)) {
  universe <- intersect(universe, sets[[i]])
}

# how many universe genes each whole-transcriptome dataset lacks (these get zero counts in step B)
zero_filled <- integer(0)
for (dataset in c("lab", "tenx", "stewart", "s1")) {
  zero_filled[dataset] <- length(setdiff(universe, maps[[dataset]]$v47))
}
ig_tr_kept <- sort(v47[ensg %in% universe & type != "protein_coding", symbol])
cat(sprintf("universe: %d genes; IG/TR constant genes kept: %s\n", length(universe), paste(ig_tr_kept, collapse = " ")))
cat("universe genes absent from (zero-filled in):", paste(names(zero_filled), zero_filled, collapse = ", "), "\n")
fwrite(data.table(ensg = universe, symbol = v47$symbol[match(universe, v47$ensg)]),
       file.path(OUT, "gene_universe.tsv"), sep = "\t")


# ---- Step 6. Harmonisation report ---------------------------------------------------------------------------
# one row per dataset: features, how many map to v47 (number and %), how many are in the universe
report_rows <- list()
for (dataset in names(maps)) {
  dataset_v47 <- maps[[dataset]]$v47
  report_rows[[dataset]] <- data.table(dataset = dataset, features = nrow(maps[[dataset]]),
                                       mapped_to_v47 = sum(!is.na(dataset_v47)),
                                       pct_mapped = round(100 * mean(!is.na(dataset_v47)), 2),
                                       in_universe = sum(dataset_v47 %in% universe))
}
report <- rbindlist(report_rows)

# then one row per gene set and a final UNIVERSE row (gene counts only)
set_rows <- data.table(dataset = c(names(sets), "UNIVERSE"), features = c(lengths(sets), length(universe)),
                       mapped_to_v47 = NA, pct_mapped = NA, in_universe = NA)
report <- rbind(report, set_rows)
fwrite(report, file.path(OUT, "harmonisation_report.tsv"), sep = "\t")
print(report)

## ---- 4. bulk TPM: moved to v2_A4_bulk_TPM_coding.R (all protein-coding genes; author 2026-10-21) ----------
