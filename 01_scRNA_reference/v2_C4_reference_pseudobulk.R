# =============================================================================
# v2_C4_reference_pseudobulk.R
#
# PURPOSE
#   Step C of REANALYSIS_PLAN_v2 (runs on Sherlock). Builds the CIBERSORTx single-cell reference files and the
#   pseudo-bulk benchmark mixtures from the final cleaned merged cHL object.
#
# STANDARD SINCE 2026-10-23 (author: "use round 14 as standard"; "make the remaining Figure 2 without downsampling")
#   Defaults now: OBJ = cHL_merged_harmony_v2_r14.rds, LABELS = output/stepC/inputs/2026-10-23_round14_labels.tsv.gz,
#   CAP_TYPE = 0 (no downsampling: every labelled cell goes into the reference), BATCH_FILTER = 0 (off),
#   INCLUDE_ALL_SOURCES = 1 (no cell is left out because of its source: the ref_exclude flags, e.g. the 157 10x cells
#   published as generic "B cell", are ignored; author: "use the cell labels of the merged reference, don't worry about
#   which sources the cells come from"), and a sixth pseudo-bulk set "Merged" (all cells of every dataset; was v2_C9).
#   The round-13 run used: OBJ=cHL_merged_harmony_v2_r13.rds LABELS=<...>/2026-10-21_round13_labels.tsv.gz CAP_TYPE=10000
#   BATCH_FILTER=0, with the 10x generic-B cells excluded (outputs in output/superseded/round13/stepC).
#
# WHICH OBJECT / LABELS (env variables)
#   OBJ    = file name in output/ (default round 14), or "rolling": the r<N>_rolling object of the rolling cell-type
#            unstable loop (author 2026-10-14).
#   LABELS = file of final per-cell labels (round 12 and later; default round 14). The labels are used as they are.
#            Set LABELS="" to fall back to the old cluster-annotation rules below.
#   Round 13 (author 2026-10-21, FINAL object):
#            OBJ=cHL_merged_harmony_v2_r13.rds  LABELS=<...>/2026-10-21_round13_labels.tsv.gz
#            (columns final_label / cluster / cluster_type; T_MKI67 is now "Cycling" = cycling T + cycling GC B).
#   Round 12 (author 2026-10-20):
#            OBJ=cHL_merged_harmony_v2_r12.rds  LABELS=<...>/2026-10-20_round12_labels.tsv.gz
#            -> labels taken as they are (r12 cluster annotation + sub-UMAP T_gd; HRS / ILC / cycling B_GC follow their
#            clusters); the sub-UMAP removal and relabelling below are skipped because r12 already contains them.
#
# LABEL RULES (author 2026-10-15)
#   - CLUSTER ANNOTATION FOR ALL CELLS. Every cell (S1 and non-S1) takes the cell type of its final cluster.
#   - Published labels serve only to annotate the clusters: cluster type = majority published label of the cluster's
#     non-S1 cells. (The loop annotation is not used, because no r3 cluster had an ILC majority.)
#   - Exception: protected rare types keep their own published label (HRS, T_gd, ILC; T_gd never forms a cluster).
#   - S1's ex vivo labels are not used.
#   - T_Treg is its own type (never merged with T_CD4).
#   - CD8 / NK / ILC / T_gd: refined in the cluster 4/8/16 sub-UMAP (945 cells removed; T_gd from sub-UMAP clusters),
#     see Step 0.
#
# WHAT THE SCRIPT DOES
#   0. Counts on the step-A gene universe (17,807 protein-coding + 16 IG/TR constant genes). The 16 IG/TR genes are
#      added here from each dataset's input object for the r10 cells, because r10 was built on the protein-coding-only
#      universe.
#   1. (C1) Reference cells: cluster annotation (protected types: own label); the 10x generic B_cell cells are left out
#      (S2D convention); at most 10,000 cells per type, simple random sample (author 2026-10-21; the plan had said
#      5,000 + 25 % per sample).
#   2. (C2) Batch-confounded genes, removed before variable-gene selection (rule and count reported).
#      Author 2026-10-21: OFF (BATCH_FILTER=0), as v1 - the rule removed canonical markers (CD163, C1QA, LYZ, MS4A1,
#      FOXP3, GNLY, CCL17, ...).
#   3. 4,000 vst variable genes, vst per dataset (split layers, as v1 split scRNA / snRNA). All are universe genes, so
#      all are in the bulk TPM (v1: 605 of 4,000 were not).
#   4. (C3) Reference files: raw UMI counts (signature matrix + HiRes refsample) and log-normalised (shown once, as v1).
#      Merged arm only since 2026-10-23 (author); round 13 also had lab_only (C6 control) and merged_noS1 (sensitivity).
#   5. (C5) Pseudo-bulk benchmarks: cHL (lab sc + sn), S1, Aoki, 10x, Stewart (author 2026-10-20). v1 make_pseudobulk()
#      unchanged (0-40 % in 2 % steps, 10 bootstraps with set.seed(1..10), CPM); labelled cells only; all universe genes;
#      10x B_cell kept as "B_cell" (S2D).
#   Gene identifiers in every output: v47 symbols (unique within the universe), as v1.
#
# OTHER SWITCHES (env variables, author 2026-10-21)
#   REUSE_COUNTS    = path of a saved counts + labels rds -> skips Step 0.
#   BATCH_FILTER    = 0 turns the C2 batch-gene filter off (as v1; it removed canonical markers such as CD163 / MS4A1 /
#                     FOXP3 because ambient RNA differs between datasets). Default 0 (off) since 2026-10-23.
#   SKIP_PSEUDOBULK = 1 writes the references only.
#   CAP_TYPE        = max cells per cell type; 0 = no cap (default since 2026-10-23; the round-13 run used 10000).
#   INCLUDE_ALL_SOURCES = 1 (default) ignores the ref_exclude flags; 0 keeps the old exclusion of flagged cells.
#   CAP_SAMPLE      = per-sample cap (default 1 = no cap). Values < 1 restore the capped-proportional allocation.
#   SUBUMAP_REFINE  = folder with the cluster 4/8/16 sub-UMAP files (default output/stepC/inputs).
#   OUT_SUBDIR      = output folder relative to W (default output/stepC). Round 16 (2026-10-26): output/r16/stepC.
#
# INPUTS : output/<OBJ>, LABELS file, ref/gene_table_v47.tsv, ref/gene_universe.tsv, ref/feature_map_<dataset>.tsv,
#          input/{lab_scRNA_seurat,aoki_filtered,stewart_filtered,tenx_filtered}.rds, the S1 f6 object (S1_OBJ)
# OUTPUTS: W/output/stepC/ ... reference/, pseudobulk/, qc/, cHL_v2_r10_counts_universe.rds
# RUN    : sbatch v2_run.sh v2_C4_reference_pseudobulk.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

W      <- "${PROJECT_ROOT}/cHL_v2"     # project root on OAK
REF    <- file.path(W, "ref")                                       # gene tables
OUT    <- file.path(W, Sys.getenv("OUT_SUBDIR", "output/stepC"))    # all outputs of this step (OUT_SUBDIR: author
                                                                    #   2026-10-26, round 16 writes to output/r16/stepC
                                                                    #   so the round-14 files are not overwritten)
S1_OBJ <- "${PROJECT_ROOT}/scRNA/2026-10-04_cHL_ex_vivo/output/2026-10-04_chl_exvivo_final_seurat_obj_f6.rds"
source(file.path(W, "script/utils_scRNA.R"))                 # v1 helpers: make_pseudobulk(), write_pseudobulk(), write_dense_tsv()

# Print a time-stamped log line
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

for (sub_dir in c("reference", "pseudobulk", "qc")) {
  dir.create(file.path(OUT, sub_dir), recursive = TRUE, showWarnings = FALSE)
}

## Caps (author 2026-10-21): at most CAP_TYPE cells per type, simple random sample (seed per type).
## The per-sample cap of the plan (25 %) is REMOVED by default (CAP_SAMPLE = 1).
## Setting CAP_SAMPLE < 1 restores the capped-proportional allocation of Step 1.
CAP_TYPE   <- as.integer(Sys.getenv("CAP_TYPE", "0"))
if (CAP_TYPE <= 0) {
  CAP_TYPE <- .Machine$integer.max                           # no cap: every labelled cell is kept (no random draw happens)
}
CAP_SAMPLE <- as.numeric(Sys.getenv("CAP_SAMPLE", "1"))
N_HVG      <- 4000                                            # number of variable genes in the reference
if (CAP_SAMPLE >= 1) {
  sample_cap_text <- "none (simple random sample)"
} else {
  sample_cap_text <- sprintf("%.0f%%", 100 * CAP_SAMPLE)
}
log_("caps: %d cells per type; per-sample cap %s", CAP_TYPE, sample_cap_text)

# Gene universe (Ensembl IDs) and their v47 symbols; symbols must be complete and unique
v47      <- fread(file.path(REF, "gene_table_v47.tsv"))
universe <- fread(file.path(REF, "gene_universe.tsv"))$ensg
sym      <- v47$symbol[match(universe, v47$ensg)]
stopifnot(!anyNA(sym), !anyDuplicated(sym))


## ---- Step 0. Counts on the gene universe + the 16 IG/TR constant genes ----------------------------------------
REUSE_COUNTS <- Sys.getenv("REUSE_COUNTS", "")
BATCH_FILTER <- Sys.getenv("BATCH_FILTER", "0") == "1"
SKIP_PB      <- Sys.getenv("SKIP_PSEUDOBULK", "0") == "1"

if (nzchar(REUSE_COUNTS)) {
  # Step 0 (shortcut): read the counts + cell table saved by an earlier run
  saved      <- readRDS(REUSE_COUNTS)
  cnt        <- saved$counts
  cell_table <- data.table(saved$meta)
  rm(saved)
  invisible(gc())
  stopifnot(identical(colnames(cnt), cell_table$cell), identical(rownames(cnt), sym))
  log_("reused counts %s: %d genes x %d cells", basename(REUSE_COUNTS), nrow(cnt), ncol(cnt))

} else {
  # Step 0a. Read the final object ("rolling" = the single r<N>_rolling object in output/)
  OBJ <- Sys.getenv("OBJ", "cHL_merged_harmony_v2_r14.rds")       # round 14 = standard (author 2026-10-23)
  if (OBJ == "rolling") {
    OBJ <- list.files(file.path(W, "output"), pattern = "^cHL_merged_harmony_v2_r[0-9]+_rolling[.]rds$")
    stopifnot(length(OBJ) == 1)
  }
  log_("final object: %s", OBJ)
  obj <- readRDS(file.path(W, "output", OBJ))

  # One row per cell, with the columns we need (the loop's cell_type is kept as loop_type)
  all_meta   <- data.table(obj@meta.data, keep.rownames = "cell")
  cell_table <- all_meta[, .(cell, dataset, sample, chemistry, label, ref_exclude, s1_cell_type, loop_type = cell_type)]
  rm(all_meta)

  # round 12+ (author 2026-10-20): per-cell labels already final in a file -> use them as they are; default round 14
  LABELS <- Sys.getenv("LABELS", file.path(W, "output/stepC/inputs/2026-10-23_round14_labels.tsv.gz"))
  if (nzchar(LABELS)) {
    # Step 0b (round 12+). Final labels come from the LABELS file, unchanged
    label_file <- fread(cmd = sprintf("gzip -dc %s", LABELS))
    stopifnot(setequal(label_file$cell, cell_table$cell))
    if ("r12_label" %in% names(label_file)) {           # round-12 file uses r12_* column names
      setnames(label_file, c("r12_label", "r12_cluster", "r12_cluster_type"), c("final_label", "cluster", "cluster_type"))
    }
    stopifnot(all(c("final_label", "cluster", "cluster_type") %in% names(label_file)))
    cell_table[, label_published := label]               # keep the published label for the QC tables
    # Join on cell: label <- final_label, final_cluster <- cluster, cluster_type <- cluster_type
    cell_table[label_file, on = "cell", `:=`(label         = i.final_label,
                                             final_cluster = as.character(i.cluster),
                                             cluster_type  = i.cluster_type)]
    cell_table[dataset == "s1", ref_exclude := NA_character_]     # S1 f6 state flag no longer applies
    # Author 2026-10-23: use every labelled cell whatever its source -> clear all ref_exclude flags
    #   (the 10x cells published as generic "B cell" then enter the reference and the pseudo-bulk with their final labels)
    if (Sys.getenv("INCLUDE_ALL_SOURCES", "1") == "1") {
      log_("INCLUDE_ALL_SOURCES: %d ref_exclude flags ignored", sum(!is.na(cell_table$ref_exclude)))
      cell_table[, ref_exclude := NA_character_]
    }
    # Table of cluster -> cluster type, in numeric cluster order
    ctype <- unique(cell_table[, .(final_cluster, cluster_type)])
    ctype <- ctype[order(as.integer(final_cluster))]
    log_("labels from %s: %d cells; T_gd %d (%d S1)", basename(LABELS), nrow(cell_table),
         cell_table[label == "T_gd", .N], cell_table[label == "T_gd" & dataset == "s1", .N])

  } else {
    # Step 0b (before round 12). Cluster annotation from the published labels
    cell_table[, final_cluster := as.character(obj$RNA_snn_res.0.6)]
    # Count published labels per cluster, using non-S1 cells with a label only
    label_counts <- cell_table[dataset != "s1" & !is.na(label), .N, by = .(final_cluster, label)]
    label_counts <- label_counts[order(final_cluster, -N)]          # most common label first within each cluster
    # Cluster type = the most common label; also record its %, the runner-up and the number of non-S1 cells
    ctype <- label_counts[, .(cluster_type = label[1],
                              pct          = round(100 * N[1] / sum(N), 1),
                              second       = paste(label[2], round(100 * N[2] / sum(N), 1)),
                              nonS1        = sum(N)),
                          by = final_cluster]
    stopifnot(all(unique(cell_table$final_cluster) %in% ctype$final_cluster))     # every final cluster has non-S1 cells
    # Attach the cluster type to every cell, then put the cells back in object order
    cell_table <- merge(cell_table, ctype[, .(final_cluster, cluster_type)], by = "final_cluster", sort = FALSE)
    cell_table <- cell_table[match(colnames(obj), cell)]

    # Cluster annotation for all cells; protected rare types (non-S1) keep their own published label
    PROTECTED <- c("HRS", "T_gd", "ILC")
    cell_table[, label_published := label]
    cell_table[, label := fifelse(dataset != "s1" & label_published %in% PROTECTED, label_published, cluster_type)]
    cell_table[dataset == "s1", ref_exclude := NA_character_]                       # S1 f6 state flag no longer applies

    ## CD8 / NK / ILC / T_gd refinement (author 2026-10-15/16; v2_C5c/C5f/C5h):
    ##   (1) remove the 945 cells dropped in the cluster 4/8/16 sub-UMAP cleanup rounds
    ##       (r2: 253 heat-shock CD8; r3: 104 MT-dominated; r4: 588 = generic weak + CD4 T cells in the CD8 cluster);
    ##   (2) take those cells' final labels from the sub-UMAP: round-4 subcluster 14 -> T_gd, subcluster 3
    ##       sub-subclusters 2/4/6 -> T_gd, all other cells their cluster annotation (protected types keep their own label).
    SUBR <- Sys.getenv("SUBUMAP_REFINE", file.path(W, "output/stepC/inputs"))
    subumap_removed <- fread(cmd = sprintf("gzip -dc %s", file.path(SUBR, "2026-10-16_cl4_8_16_removed_cells_all_rounds.tsv.gz")))$cell
    subumap_labels  <- fread(cmd = sprintf("gzip -dc %s", file.path(SUBR, "2026-10-16_cl4_8_16_final_labels.tsv.gz")))
    stopifnot(length(subumap_removed) == 945,
              nrow(subumap_labels) == 13619,
              all(subumap_removed %in% cell_table$cell),
              all(subumap_labels$cell %in% cell_table$cell),
              !any(subumap_labels$cell %in% subumap_removed))
    cell_table <- cell_table[!cell %in% subumap_removed]
    cell_table[subumap_labels, on = "cell", label := i.final_label]       # sub-UMAP label replaces the cluster label
    log_("sub-UMAP refinement: %d cells removed, %d CD8/NK/ILC/T_gd labels set (T_gd now %d cells, %d of them S1)",
         length(subumap_removed), nrow(subumap_labels),
         cell_table[label == "T_gd", .N], cell_table[label == "T_gd" & dataset == "s1", .N])
  }

  # Step 0c. QC tables: non-S1 cells whose published label differs from the final label
  relabelled <- cell_table[dataset != "s1" & !is.na(label_published) & label_published != label,
                           .N, by = .(dataset, label_published, label)]
  relabelled <- relabelled[order(-N)]
  fwrite(relabelled, file.path(OUT, "qc", "published_vs_cluster_label.tsv"), sep = "\t")
  log_("non-S1 cells relabelled published -> cluster type: %d of %d",
       sum(relabelled$N), cell_table[dataset != "s1" & !is.na(label_published), .N])
  print(head(relabelled, 15))
  fwrite(ctype, file.path(OUT, "qc", "final_cluster_types.tsv"), sep = "\t")
  print(ctype)
  # S1: its ex vivo cell type vs the final label
  s1_table <- cell_table[dataset == "s1", .N, by = .(s1_cell_type, label)]
  s1_table <- s1_table[order(s1_cell_type, -N)]
  fwrite(s1_table, file.path(OUT, "qc", "S1_exvivo_type_vs_cluster_label.tsv"), sep = "\t")
  print(s1_table)
  cell_table[, c("s1_cell_type", "loop_type", "final_cluster", "cluster_type", "label_published") := NULL]

  # Step 0d. Raw counts of the kept cells (cells removed above are dropped here too)
  cnt <- LayerData(obj, assay = "RNA", layer = "counts")[, cell_table$cell]
  rm(obj)
  invisible(gc())
  extra         <- setdiff(universe, rownames(cnt))               # universe genes missing from the object (the IG/TR genes)
  extra_symbols <- v47$symbol[match(extra, v47$ensg)]
  log_("r10: %d genes x %d cells; adding %d genes: %s", nrow(cnt), ncol(cnt),
       length(extra), paste(extra_symbols, collapse = " "))

  # Helper: counts of the 'extra' genes for one dataset (features x cells -> extra x cells).
  # Each input feature is mapped to its v47 gene; features mapping to the same gene are summed.
  extra_rows <- function(counts_matrix, dataset_name) {
    feature_map <- fread(file.path(REF, sprintf("feature_map_%s.tsv", dataset_name)))
    feature_map <- feature_map[match(rownames(counts_matrix), feature)]    # one row per matrix feature, same order
    is_extra    <- !is.na(feature_map$v47) & feature_map$v47 %in% extra
    # indicator matrix: input feature (row) -> extra gene (column)
    to_gene <- sparseMatrix(i = which(is_extra), j = match(feature_map$v47[is_extra], extra), x = 1,
                            dims = c(nrow(counts_matrix), length(extra)))
    result <- as(crossprod(to_gene, counts_matrix), "CsparseMatrix")
    rownames(result) <- extra
    result
  }

  # Step 0e. Extra-gene counts from each dataset's input object (cell names get the dataset prefix used in the object)
  extra_by_dataset <- list()

  # lab
  input_obj <- readRDS(file.path(W, "input/lab_scRNA_seurat.rds"))
  input_counts <- LayerData(input_obj, assay = "RNA", layer = "counts")
  rm(input_obj)
  colnames(input_counts) <- paste0("lab_", colnames(input_counts))
  extra_by_dataset$lab <- extra_rows(input_counts, "lab")
  rm(input_counts)
  invisible(gc())

  # Aoki (counts are in the "originalexp" assay)
  input_obj <- readRDS(file.path(W, "input/aoki_filtered.rds"))
  input_counts <- LayerData(input_obj, assay = "originalexp", layer = "counts")
  rm(input_obj)
  colnames(input_counts) <- paste0("aoki_", colnames(input_counts))
  extra_by_dataset$aoki <- extra_rows(input_counts, "aoki")
  rm(input_counts)
  invisible(gc())

  # Stewart: the stored "counts" are log-normalised. Back to raw UMI (B0, as step B):
  #   expm1 -> scale each cell so it sums to its total_counts -> round -> drop zeros
  input_obj <- readRDS(file.path(W, "input/stewart_filtered.rds"))
  stewart_lognorm <- LayerData(input_obj, assay = "RNA", layer = "counts")
  total_counts    <- input_obj$total_counts
  rm(input_obj)
  unlogged <- stewart_lognorm
  unlogged@x <- expm1(unlogged@x)
  stewart_raw <- unlogged %*% Diagonal(x = total_counts / colSums(unlogged))
  stewart_raw@x <- round(stewart_raw@x)
  stewart_raw <- drop0(stewart_raw)
  dimnames(stewart_raw) <- list(rownames(stewart_lognorm), paste0("stewart_", colnames(stewart_lognorm)))
  extra_by_dataset$stewart <- extra_rows(stewart_raw, "stewart")
  rm(stewart_lognorm, unlogged, stewart_raw)
  invisible(gc())

  # 10x
  input_obj <- readRDS(file.path(W, "input/tenx_filtered.rds"))
  input_counts <- LayerData(input_obj, assay = "RNA", layer = "counts")
  rm(input_obj)
  colnames(input_counts) <- paste0("tenx_", colnames(input_counts))
  extra_by_dataset$tenx <- extra_rows(input_counts, "tenx")
  rm(input_counts)
  invisible(gc())

  # S1 (only Sample1 of the ex vivo object)
  input_obj <- readRDS(S1_OBJ)
  input_obj <- subset(input_obj, subset = sample == "Sample1")
  input_counts <- LayerData(input_obj, assay = "RNA", layer = "counts")
  rm(input_obj)
  colnames(input_counts) <- paste0("s1_", colnames(input_counts))
  extra_by_dataset$s1 <- extra_rows(input_counts, "s1")
  rm(input_counts)
  invisible(gc())

  # Step 0f. Add the extra genes to the counts, put genes in universe order, rename rows to v47 symbols
  extra_counts <- do.call(cbind, extra_by_dataset)
  rm(extra_by_dataset)
  stopifnot(all(colnames(cnt) %in% colnames(extra_counts)))
  cnt <- rbind(cnt, extra_counts[, colnames(cnt)])
  cnt <- cnt[universe, ]
  rm(extra_counts)
  invisible(gc())
  rownames(cnt) <- sym
  # % of each dataset's UMI that falls in the added IG/TR genes
  print(cell_table[, .(cells = .N,
                       IG_TR_UMI_pct = round(100 * sum(cnt[extra_symbols, cell]) / sum(cnt[, cell]), 2)),
                   by = dataset])
  saveRDS(list(counts = cnt, meta = cell_table), file.path(OUT, "cHL_v2_r10_counts_universe.rds"))
  log_("counts: %d genes x %d cells saved", nrow(cnt), ncol(cnt))
}


## ---- Step 1. Reference cells (C1) ---------------------------------------------------------------------------------
## Default (CAP_SAMPLE >= 1, author 2026-10-21): per type, keep every cell if the type has <= CAP_TYPE cells,
## otherwise a simple random sample of CAP_TYPE cells (seed = 42 + the type's position in the alphabetical type list).
##
## Optional per-sample cap (CAP_SAMPLE < 1; the plan's rule): per type, T = min(CAP_TYPE, available) cells, sampled in
## proportion to each sample's size but with no sample above 25 % of T (capped proportional allocation; the excess of
## capped samples is redistributed in proportion to the others). When the cap cannot be met (fewer than 4 samples, or
## one sample dominating), T is lowered to the largest feasible total; with fewer than 4 samples the cap is
## 1 / (number of samples). Fixed 2026-10-14: the first version did not limit the total.

# Helper (per-sample-cap mode only): how many cells to take from each sample.
#   sample_sizes = cells per sample, total = cells wanted, cap = most cells any one sample may give
allocate_quota <- function(sample_sizes, total, cap) {
  limit <- pmin(sample_sizes, floor(cap))
  quota <- numeric(length(sample_sizes))
  fixed <- rep(FALSE, length(sample_sizes))              # samples already set to their limit
  # share the remaining total among the not-yet-fixed samples in proportion to their size;
  # any sample that goes over its limit is fixed at the limit, then repeat
  repeat {
    remaining <- total - sum(quota[fixed])
    weight <- ifelse(fixed, 0, sample_sizes)
    quota[!fixed] <- remaining * weight[!fixed] / sum(weight[!fixed])
    over <- !fixed & quota > limit
    if (!any(over)) break
    quota[over] <- limit[over]
    fixed <- fixed | over
  }
  # round down, then hand out the cells lost to rounding, samples with the most room first
  quota <- floor(quota)
  short <- total - sum(quota)
  room <- limit - quota
  for (i in order(-room)) {
    if (short <= 0) break
    add <- min(room[i], short)
    quota[i] <- quota[i] + add
    short <- short - add
  }
  quota
}

# Helper: pick the reference cells of ONE cell type. type_rows = that type's rows of the eligible-cell table,
# in object order; seed = 42 + the type's position in the list.
pick_type <- function(type_rows, seed) {
  if (CAP_SAMPLE >= 1) {
    # default: simple random sample, no per-sample cap
    set.seed(seed)                                   # seed is set for every type (also when no sampling is needed)
    type_cells <- type_rows$cell
    if (length(type_cells) <= CAP_TYPE) {
      return(type_cells)                             # small type: keep every cell
    }
    return(sample(type_cells, CAP_TYPE))             # large type: random CAP_TYPE cells
  }

  # per-sample cap mode (CAP_SAMPLE < 1)
  cells_per_sample <- type_rows[, .N, by = sample]                      # samples in order of first appearance
  cap_fraction <- max(CAP_SAMPLE, 1 / nrow(cells_per_sample))           # with few samples the cap is 1 / #samples
  # a total is feasible when the samples, each capped at cap_fraction * total, can supply it
  is_feasible <- function(total) sum(pmin(cells_per_sample$N, floor(cap_fraction * total))) >= total
  total <- min(CAP_TYPE, sum(cells_per_sample$N))
  while (total > 0 && !is_feasible(total)) {
    total <- total - 1                       # largest feasible total (rounding makes feasibility non-monotone)
  }
  quota <- allocate_quota(cells_per_sample$N, total, cap_fraction * total)
  stopifnot(sum(quota) == total, all(quota <= cells_per_sample$N))
  set.seed(seed)
  picked <- c()
  for (k in seq_len(nrow(cells_per_sample))) {
    this_sample <- cells_per_sample$sample[k]
    sample_cells <- type_rows[sample == this_sample, cell]
    if (quota[k] >= length(sample_cells)) {
      picked <- c(picked, sample_cells)                       # quota covers the whole sample: keep all
    } else {
      picked <- c(picked, sample(sample_cells, quota[k]))     # random quota[k] cells of this sample
    }
  }
  picked
}

# Choose the reference cells of one arm and write the per-type summary table
select_ref <- function(arm_cells, tag) {
  # every cell with a final label, EXCEPT cells with a ref_exclude flag (e.g. the 10x cells published as generic "B cell";
  # S1's flag was cleared in Step 0)
  eligible <- arm_cells[!is.na(label) & is.na(ref_exclude)]
  cell_types <- sort(unique(eligible$label))           # alphabetical order
  chosen <- c()
  for (i in seq_along(cell_types)) {
    type_rows <- eligible[label == cell_types[i]]
    picked <- pick_type(type_rows, seed = 42 + i)      # seed = 42 + the type's position in the list
    chosen <- c(chosen, picked)
  }
  selected <- eligible[cell %in% chosen]                # selected cells, in object order

  # per type: available cells / samples, and selected cells / largest sample share / cells per dataset
  available_per_type <- eligible[, .(available = .N, samples = uniqueN(sample)), by = label]
  selected_per_type  <- selected[, .(selected       = .N,
                                     max_sample_pct = round(100 * max(table(sample)) / .N, 1),
                                     datasets       = paste(names(table(dataset)), table(dataset), sep = "=", collapse = " ")),
                                 by = label]
  summary_table <- merge(available_per_type, selected_per_type, by = "label")
  fwrite(summary_table, file.path(OUT, "qc", sprintf("reference_cells_%s.tsv", tag)), sep = "\t")
  cat("\n==", tag, "\n")
  print(summary_table)
  selected
}


## ---- Step 2. Batch-confounded genes (C2) --------------------------------------------------------------------------
## Per cell type, pseudo-bulk log2(CPM + 1) per dataset (datasets with >= 30 reference cells of that type). A gene is
## evaluable in a type when its highest dataset CPM >= 10; it is batch-confounded when it is evaluable in >= 3 types and
## its between-dataset range is > 3 (8-fold) in >= 50 % of them. Cell-type markers are evaluable in few types, so they
## are not removed by this rule; genes systematically shifted by dataset / chemistry are.
batch_genes <- function(ref_cells, tag) {
  # cell type x dataset groups with at least 30 reference cells
  groups <- ref_cells[, .N, by = .(label, dataset)]
  groups <- groups[N >= 30]

  # one CPM pseudo-bulk column per group
  pseudobulk_list <- list()
  for (g in seq_len(nrow(groups))) {
    group_type    <- groups$label[g]
    group_dataset <- groups$dataset[g]
    group_cells <- ref_cells[label == group_type & dataset == group_dataset, cell]
    gene_sums <- Matrix::rowSums(cnt[, group_cells, drop = FALSE])
    pseudobulk_list[[g]] <- 1e6 * gene_sums / sum(gene_sums)
  }
  cpm <- do.call(cbind, pseudobulk_list)

  # per type with >= 2 datasets: each gene's highest CPM and its log2 range between datasets
  per_type_list <- list()
  for (this_type in unique(groups$label)) {
    type_columns <- which(groups$label == this_type)
    if (length(type_columns) < 2) next                                  # one dataset only: nothing to compare
    type_cpm <- cpm[, type_columns, drop = FALSE]
    per_type_list[[length(per_type_list) + 1]] <- data.table(
      gene       = rownames(cpm),
      type       = this_type,
      max_cpm    = apply(type_cpm, 1, max),
      range_log2 = apply(log2(type_cpm + 1), 1, function(v) diff(range(v))))     # max - min of log2(CPM + 1)
  }
  per_type <- rbindlist(per_type_list)
  if (!nrow(per_type)) {
    log_("%s: one dataset only - no batch-gene filter", tag)
    return(character(0))
  }

  # per gene, over the types where it is evaluable (max CPM >= 10)
  gene_scores <- per_type[max_cpm >= 10, .(types_evaluable = .N,
                                           frac_shifted    = mean(range_log2 > 3),
                                           median_range    = median(range_log2)),
                          by = gene]
  removed_genes <- gene_scores[types_evaluable >= 3 & frac_shifted >= 0.5, gene]
  fwrite(gene_scores[order(-frac_shifted, -median_range)],
         file.path(OUT, "qc", sprintf("batch_gene_scores_%s.tsv", tag)), sep = "\t")
  log_("%s: %d batch-confounded genes removed before HVG (of %d evaluable; types with >= 2 datasets: %d)", tag,
       length(removed_genes), nrow(gene_scores), length(unique(per_type$type)))
  removed_genes
}


## ---- Steps 3-4. Variable genes + reference files (C3) ---------------------------------------------------------------
build_reference <- function(arm_cells, tag) {
  ref_cells <- select_ref(arm_cells, tag)

  # Step 2 (optional): batch-confounded genes
  if (BATCH_FILTER) {
    removed_genes <- batch_genes(ref_cells, tag)
  } else {
    log_("%s: batch-gene filter OFF", tag)
    removed_genes <- character(0)
  }
  keep_genes <- setdiff(rownames(cnt), removed_genes)

  # Step 3: 4,000 vst variable genes, vst run per dataset (one layer per dataset)
  vf <- CreateSeuratObject(counts = cnt[keep_genes, ref_cells$cell],
                           meta.data = data.frame(dataset = ref_cells$dataset, row.names = ref_cells$cell))
  if (uniqueN(ref_cells$dataset) > 1) {
    vf[["RNA"]] <- split(vf[["RNA"]], f = vf$dataset)
  }
  vf <- FindVariableFeatures(vf, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
  hvg <- VariableFeatures(vf)
  rm(vf)
  invisible(gc())
  stopifnot(length(hvg) == N_HVG, all(hvg %in% sym))
  fwrite(data.table(gene = hvg, ensg = universe[match(hvg, sym)]),
         file.path(OUT, "qc", sprintf("hvg_%d_%s.tsv", N_HVG, tag)), sep = "\t")

  # Step 4: reference files (genes x cells; column header = each cell's label)
  ref_counts <- cnt[hvg, ref_cells$cell]
  write_dense_tsv(ref_counts, file.path(OUT, "reference", sprintf("cHL_v2_scRNA_reference_%d_%s.txt", N_HVG, tag)),
                  header = c("gene", ref_cells$label), writer = "fwrite")
  # log-normalised version = NormalizeData(LogNormalize, scale 1e4): log1p(count / cell total over ALL genes * 1e4)
  ref_lognorm <- ref_counts
  cell_totals <- Matrix::colSums(cnt[, ref_cells$cell])
  ref_lognorm@x <- log1p(ref_lognorm@x / rep(cell_totals, diff(ref_lognorm@p)) * 1e4)   # each non-zero / its cell's total
  write_dense_tsv(ref_lognorm, file.path(OUT, "reference", sprintf("cHL_v2_scRNA_reference_%d_%s_lognorm.txt", N_HVG, tag)),
                  header = c("gene", ref_cells$label), writer = "fwrite")
  fwrite(ref_cells[, .(cell, label, dataset, sample)],
         file.path(OUT, "qc", sprintf("reference_cell_ids_%s.tsv.gz", tag)), sep = "\t")
  log_("%s reference: %d genes x %d cells, %d types", tag, length(hvg), nrow(ref_cells), uniqueN(ref_cells$label))
}

build_reference(cell_table, "merged")                          # primary
# Author 2026-10-23: only the merged reference is used ("remove lab_only and merged_noS1").
#   The round-13 run also built lab_only (cell_table[dataset == "lab"]) and merged_noS1 (cell_table[dataset != "s1"]).


## ---- Step 5. Pseudo-bulk benchmarks (C5) ----------------------------------------------------------------------------------
if (SKIP_PB) {
  log_("pseudo-bulk skipped (SKIP_PSEUDOBULK=1)")
  quit(save = "no")
}
# 10x cells published as generic "B cell" (flagged in ref_exclude) are kept as their own type "B_cell" (S2D convention)
cell_table[, pb_label := fifelse(dataset == "tenx" & !is.na(ref_exclude) & grepl("^10x generic B_cell", ref_exclude),
                                 "B_cell", label)]
## Five benchmark sets (author 2026-10-20; supersedes "cHL = lab + S1" of 2026-10-13): cHL = lab only (scRNA + snRNA),
## S1 alone (new dataset, Flex chemistry), Aoki, 10X, Stewart (v1 names). Types are not filtered by size (as v1);
## the cell count per type is logged.
## Author 2026-10-23: pseudo-bulk benchmark on the MERGED set only (every cell of every dataset), not per source.
##   Until then: cHL = "lab", S1 = "s1", Aoki = "aoki", 10X = "tenx", Stewart = "stewart" (+ Merged from v2_C9).
SETS <- list(Merged = c("lab", "s1", "aoki", "tenx", "stewart"))
for (set_name in names(SETS)) {
  set_cells <- cell_table[dataset %in% SETS[[set_name]] & !is.na(pb_label)]
  set_dir <- file.path(OUT, "pseudobulk", paste0("pseudobulk_", set_name))
  dir.create(set_dir, showWarnings = FALSE)
  set_counts <- cnt[, set_cells$cell]
  type_sizes <- table(set_cells$pb_label)
  # one mixture file per cell type (spike-in of that type into all other cells of the set)
  for (cell_type in names(type_sizes)) {
    mixtures <- make_pseudobulk(set_counts, set_cells$pb_label, cell_type)
    write_pseudobulk(mixtures, file.path(set_dir, sprintf("pseudobulk_%s_%s_cpm.txt", set_name, cell_type)))
  }
  log_("pseudo-bulk %s: %d cells, %d types (%s)", set_name, nrow(set_cells), uniqueN(set_cells$pb_label),
       paste(names(type_sizes), type_sizes, sep = "=", collapse = " "))
  rm(set_counts)
  invisible(gc())
}
log_("step C done")
