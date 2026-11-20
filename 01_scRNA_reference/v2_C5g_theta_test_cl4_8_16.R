# =============================================================================
# v2_C5g_theta_test_cl4_8_16.R
#
# Purpose (author 2026-10-15): Harmony theta test on the round-3 cluster 4 / 8 / 16 sub-UMAP.
#   S1 integrates poorly with the default theta 2; test a higher theta here before using it elsewhere.
#   Same cells and PCA as round 3 (14,207 cells); only Harmony changes. Four settings:
#     A  dataset + chemistry + sample, theta 2 / 2 / 2   (current)
#     B  dataset + chemistry + sample, theta 4 / 4 / 2
#     C  dataset + chemistry + sample, theta 8 / 8 / 4
#     D  chemistry + sample,           theta 12 / 6      (author's NLPHL setting)
#   Then, per setting: FindNeighbors (harmony 1:30) -> Louvain res 1.2 (seed 42) -> UMAP (seed 42).
#   (The author chose B on 2026-10-16; used from round 3B on, see v2_C5f.)
#
# Metrics per setting:
#   - S1 mixing: inverse Simpson index of dataset among each cell's 30 nearest neighbours (median over all cells and over
#     S1 cells), and the share of S1 cells' neighbours that are S1
#   - subclusters that are > 80 % S1 (how many, and how many cells they hold)
#   - best T_gd subcluster (purity, capture, S1 cells) and best ILC subcluster (purity, capture)
#
# Usage: Rscript v2_C5g_... <round-3 obj.rds> <Box table dir> <Box figure dir> [point size, default 0.6]
#   (runs locally; nothing overwritten)
# 2026-10-16: bigger dots (author: 0.15 was too small) and UMAP coordinates saved per setting (no Aoki counts, so OK in Box).
#
# Output: <Box table dir>/2026-10-15_r11rolling_sub_cl4_8_16_r3_theta_test.tsv
#         <Box table dir>/2026-10-15_r11rolling_sub_cl4_8_16_r3_theta_test_umap_coords.tsv.gz
#         <Box figure dir>/2026-10-16_r11rolling_sub_cl4_8_16_r3_theta_test_bigdots.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})
args       <- commandArgs(TRUE)
sub        <- readRDS(args[1])       # round-3 sub-UMAP object (has the PCA)
table_dir  <- args[2]
figure_dir <- args[3]
stamp      <- "2026-10-15"           # date prefix of the two tables (the PDF name below is fixed to 2026-10-16)
point_size <- if (length(args) >= 4) as.numeric(args[4]) else 0.6

# The four Harmony settings: covariates and theta (one theta per covariate)
settings <- list(
  A = list(covariates = c("dataset", "chemistry", "sample"), theta = c(2, 2, 2)),
  B = list(covariates = c("dataset", "chemistry", "sample"), theta = c(4, 4, 2)),
  C = list(covariates = c("dataset", "chemistry", "sample"), theta = c(8, 8, 4)),
  D = list(covariates = c("chemistry", "sample"),            theta = c(12, 6))
)

published_label <- ifelse(sub$dataset == "s1", NA, sub$label)        # S1 has no published label
is_s1           <- sub$dataset == "s1"
label_colors    <- c(ILC = "#00C19A", NK = "#FDBF6F", T_CD8 = "#999999", T_gd = "#A65628")

# Helper: the subcluster with the highest share of one published label (purity), and the share of all cells with that
#         label it holds (capture), as a text line
best_subcluster <- function(cell_info, label_name) {
  per_cluster <- cell_info[, .(n = .N, k = sum(published_label %in% label_name), s1 = sum(is_s1)), by = cluster]
  per_cluster <- per_cluster[k > 0]                                  # subclusters with at least one such cell
  per_cluster[, purity := k / n]
  per_cluster[, capture := k / sum(k)]
  per_cluster <- per_cluster[order(-purity)]
  best <- per_cluster[1]
  sprintf("cl%s: %.0f%% %s, captures %.0f%%, %d cells, %d S1",
          best$cluster, 100 * best$purity, label_name, 100 * best$capture, best$n, best$s1)
}

metrics     <- list()
coords      <- list()
plots       <- list()
for (setting in names(settings)) {
  covariates <- settings[[setting]]$covariates
  theta      <- settings[[setting]]$theta

  # Step 1. Harmony with this setting -> neighbours -> Louvain res 1.2 -> UMAP (seed 42 throughout)
  set.seed(42)
  obj <- RunHarmony(sub, group.by.vars = covariates, theta = theta,
                    reduction.use = "pca", reduction.save = "harmony_t", verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony_t", dims = 1:30, verbose = FALSE, graph.name = c("t_nn", "t_snn"))
  obj <- FindClusters(obj, graph.name = "t_snn", resolution = 1.2, random.seed = 42, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony_t", dims = 1:30, seed.use = 42, reduction.name = "umap_t", verbose = FALSE)
  cluster <- as.character(obj$seurat_clusters)

  # Step 2. S1 mixing: 30 nearest neighbours of every cell in the Harmony space
  neighbours <- FindNeighbors(obj, reduction = "harmony_t", dims = 1:30, k.param = 30,
                              return.neighbor = TRUE, verbose = FALSE)@neighbors[[1]]
  neighbour_index <- Indices(neighbours)                             # one row per cell: positions of its 30 neighbours
  cell_dataset    <- obj$dataset
  # inverse Simpson index of dataset among each cell's neighbours (1 = all one dataset; higher = better mixed)
  inverse_simpson <- numeric(nrow(neighbour_index))
  for (i in seq_len(nrow(neighbour_index))) {
    neighbours_i <- neighbour_index[i, ]
    proportions  <- table(cell_dataset[neighbours_i]) / length(neighbours_i)
    inverse_simpson[i] <- 1 / sum(proportions^2)
  }
  # for each S1 cell: the share of its neighbours that are also S1
  s1_rows  <- neighbour_index[is_s1, , drop = FALSE]
  s1_share <- numeric(nrow(s1_rows))
  for (i in seq_len(nrow(s1_rows))) {
    s1_share[i] <- mean(cell_dataset[s1_rows[i, ]] == "s1")
  }

  # Step 3. Subclusters that are > 80 % S1, and the best T_gd / ILC subclusters
  cell_info <- data.table(cluster, published_label, is_s1)
  s1_fraction        <- cell_info[, .(n = .N, s1 = mean(is_s1)), by = cluster]
  mostly_s1_clusters <- s1_fraction[s1 > 0.8]

  metrics[[setting]] <- data.table(
    setting                   = setting,
    covariates                = paste(covariates, collapse = "+"),
    theta                     = paste(theta, collapse = "/"),
    clusters                  = uniqueN(cluster),
    median_iLISI_all          = round(median(inverse_simpson), 2),
    median_iLISI_S1           = round(median(inverse_simpson[is_s1]), 2),
    S1_neighbours_that_are_S1 = round(median(s1_share), 2),
    clusters_over80pct_S1     = nrow(mostly_s1_clusters),
    cells_in_them             = sum(mostly_s1_clusters$n),
    best_T_gd                 = best_subcluster(cell_info, "T_gd"),
    best_ILC                  = best_subcluster(cell_info, "ILC")
  )
  print(metrics[[setting]])

  # Step 4. UMAP coordinates for this setting (saved below) and the plot label: S1 / published label / "other"
  umap_table <- data.table(Embeddings(obj, "umap_t"))
  umap_table[, dataset := obj$dataset]
  umap_table[, lab := fifelse(is_s1, "S1", fifelse(obj$label %in% names(label_colors), obj$label, "other"))]
  setnames(umap_table, 1:2, c("U1", "U2"))
  coords[[setting]] <- cbind(setting = setting, cell = colnames(obj), cluster = cluster, umap_table)

  # Step 5. Plot: UMAP coloured by dataset | by published label (S1 grey; T_gd drawn on top)
  plot_title <- sprintf("%s: %s theta %s", setting, paste(covariates, collapse = "+"), paste(theta, collapse = "/"))
  shuffled   <- umap_table[sample(.N)]                               # random drawing order (no set.seed of its own)
  dataset_plot <- ggplot(shuffled, aes(U1, U2, colour = dataset)) +
                  geom_point(size = point_size, stroke = 0, alpha = 0.7) +
                  coord_equal() +
                  theme_classic(12) +
                  ggtitle(plot_title) +
                  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))
  tgd_on_top <- umap_table[order(lab == "T_gd")]
  label_plot <- ggplot(tgd_on_top, aes(U1, U2, colour = lab)) +
                geom_point(size = point_size, stroke = 0, alpha = 0.7) +
                coord_equal() +
                theme_classic(12) +
                scale_colour_manual(values = c(label_colors, S1 = "#D9D9D9", other = "#CCCCCC")) +
                ggtitle("published label (S1 grey)") +
                guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))
  plots[[setting]] <- dataset_plot | label_plot

  rm(obj, neighbours, neighbour_index)
  invisible(gc())
}


# Step 6. Save the metrics table, the UMAP coordinates of every setting, and the plots (one page per setting)
metrics_table <- rbindlist(metrics)
fwrite(metrics_table, file.path(table_dir, paste0(stamp, "_r11rolling_sub_cl4_8_16_r3_theta_test.tsv")), sep = "\t")
print(metrics_table)
fwrite(rbindlist(coords), file.path(table_dir, paste0(stamp, "_r11rolling_sub_cl4_8_16_r3_theta_test_umap_coords.tsv.gz")), sep = "\t")

pdf(file.path(figure_dir, "2026-10-16_r11rolling_sub_cl4_8_16_r3_theta_test_bigdots.pdf"), width = 20, height = 9)
for (setting in names(plots)) {
  print(plots[[setting]])
}
dev.off()
cat("done\n")
