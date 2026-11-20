# =============================================================================
# v2_D11_r14_configs.R — list every distinct per-cell-type state configuration that EcoTyper's own step-4 rule can
# produce on the ROUND-14 discovery (v2_D10, 100 NMF restarts, ranks 2-6), for the 2- / 3-ecotype search (author 2026-10-24).
#
# Same method as the round-13 search (ECOTYPER_SEARCH_PLAN_v2.md, FROZEN RULES 2):
#   - settings = cophenetic cutoff 0.900-1.000 (step 0.0005) x maximum states per cell type;
#   - for each setting, replay EcoTyper's step 4 (pipeline/state_discovery_rank_selection.R): the LAST crossing of the
#     cutoff, then whichever of rank k-1 / k is closer to it; no crossing -> the rank before the biggest drop;
#   - a configuration = one distinct vector of chosen ranks (cell types in alphabetical order, as rank_data.txt);
#     it is run at its median setting, and the number of settings producing it is its robustness.
# Difference from round 13: the round-14 NMF has ranks 2-6 only (author: max states 6), so max states is 4, 5 or 6
#   (max 3 crashes EcoTyper step 4; max 7-8 would need new NMF ranks).
#
# Check: the replay must reproduce the real 0.95 / max 6 run (Custom/Analysis/rank_selection/rank_data.txt of v2_D10).
#
# Run on Sherlock:  module load R/4.2.0; Rscript v2_D11_r14_configs.R
# Output: output/stepD/ecotyper_r14_merged/r14_configs.tsv
#   columns: vec, n_settings, maxes, rep_max, rep_cut, total_states (same as the round-13 r100_configs.tsv)
# =============================================================================

# Step 1. Folders. Optional arguments (author 2026-10-25, DC exploration): <discovery name> <output file>;
#   without arguments = the main round-14 discovery, as first run.
ECOTYPER_DIR   <- "${PROJECT_ROOT}/Ecotyper/ecotyper-master"
arguments      <- commandArgs(TRUE)
DISCOVERY_NAME <- "2026-10-24_cHL_v2_r14_merged_discovery"
OUTPUT_FILE    <- "${PROJECT_ROOT}/cHL_v2/output/stepD/ecotyper_r14_merged/r14_configs.tsv"
if (length(arguments) == 2) {
  DISCOVERY_NAME <- arguments[1]
  OUTPUT_FILE    <- arguments[2]
}
STATES_DIR     <- file.path(ECOTYPER_DIR, "EcoTyper", DISCOVERY_NAME, "Custom", "Cell_States", "discovery")
REAL_RANK_FILE <- file.path(ECOTYPER_DIR, "EcoTyper", DISCOVERY_NAME, "Custom", "Analysis", "rank_selection", "rank_data.txt")

# Step 2. Search grid (frozen-plan cutoffs; max states limited by the NMF ranks that exist)
CUTOFFS    <- round(seq(0.900, 1.000, by = 0.0005), 4)
MAX_STATES <- 4:6
NMF_RANKS  <- 2:6

# Step 3. Cophenetic coefficient per cell type and rank (one rank_data.txt per <cell type>/<rank> folder, as step 4 reads)
cell_types <- sort(list.files(STATES_DIR))
cophenetic <- data.frame()
for (cell_type in cell_types) {
  for (rank in NMF_RANKS) {
    rank_file <- file.path(STATES_DIR, cell_type, rank, "rank_data.txt")
    if (!file.exists(rank_file)) {
      next
    }
    rank_data <- read.delim(rank_file)
    cophenetic <- rbind(cophenetic,
                        data.frame(CellType = cell_type, n_clusters = rank_data$n_clusters, Cophenetic = rank_data$Cophenetic))
  }
}
cell_types <- sort(unique(cophenetic$CellType))       # cell types with NMF results (pDC is skipped by EcoTyper if too few genes)
cat("cell types:", paste(cell_types, collapse = " "), "\n")

# Step 4. EcoTyper's step-4 rule for one cell type (rows = ranks 2..max_states, in rank order)
choose_rank <- function(one_type, cutoff, max_states) {
  one_type <- one_type[one_type$n_clusters <= max_states, ]
  one_type <- one_type[order(one_type$n_clusters), ]
  coph <- one_type$Cophenetic

  # 4a. a "crossing" at row i: coefficient i below the cutoff, the two rows before it at or above it
  crossing <- rep(FALSE, length(coph))
  for (i in 3:length(coph)) {
    crossing[i] <- coph[i] < cutoff & coph[i - 1] >= cutoff & coph[i - 2] >= cutoff
  }

  # 4b. no crossing: EcoTyper falls back to the rank before the biggest drop in the coefficient
  if (!any(crossing)) {
    drops <- c(NA, coph[-length(coph)] - coph[-1])
    return(one_type$n_clusters[which.max(drops) - 1])
  }

  # 4c. last crossing: rank k-1 or k, whichever coefficient is closer to the cutoff (k-1 if tied)
  last_cross <- max(which(crossing))
  distances <- c(abs(cutoff - coph[last_cross - 1]), abs(cutoff - coph[last_cross]))
  one_type$n_clusters[last_cross - 2 + which.min(distances)]
}

# chosen ranks of all cell types for one setting, as one "3-2-4-..." text in alphabetical cell-type order
choose_all_ranks <- function(cutoff, max_states) {
  chosen <- c()
  for (cell_type in cell_types) {
    chosen <- c(chosen, choose_rank(cophenetic[cophenetic$CellType == cell_type, ], cutoff, max_states))
  }
  paste(chosen, collapse = "-")
}

# Step 5. Check the replay against the real 0.95 / max 6 run
real_ranks <- read.delim(REAL_RANK_FILE)
real_ranks <- real_ranks[order(real_ranks$CellType), ]
real_vec   <- paste(real_ranks$Chosen_Rank, collapse = "-")
replay_vec <- choose_all_ranks(0.95, 6)
cat("real 0.95 / max 6:", real_vec, "\nreplay           :", replay_vec, "\n")
if (!identical(real_ranks$CellType, cell_types) || real_vec != replay_vec) {
  stop("the replay does not reproduce the real run; configurations not written")
}

# Step 6. Chosen ranks for every setting (max states x cutoff)
settings <- data.frame()
for (max_states in MAX_STATES) {
  for (cutoff in CUTOFFS) {
    settings <- rbind(settings, data.frame(max = max_states, cut = cutoff, vec = choose_all_ranks(cutoff, max_states)))
  }
}

# Step 7. One row per distinct configuration; its representative setting = the median one when its settings are
#   ordered by max states, then cutoff (as round 13)
configs <- data.frame()
for (vec in unique(settings$vec)) {
  its_settings <- settings[settings$vec == vec, ]
  its_settings <- its_settings[order(its_settings$max, its_settings$cut), ]
  median_row <- ceiling(nrow(its_settings) / 2)
  configs <- rbind(configs, data.frame(
    vec          = vec,
    n_settings   = nrow(its_settings),
    maxes        = paste(sort(unique(its_settings$max)), collapse = ","),
    rep_max      = its_settings$max[median_row],
    rep_cut      = its_settings$cut[median_row],
    total_states = sum(as.integer(strsplit(vec, "-")[[1]]))))
}
configs <- configs[order(-configs$n_settings), ]

# Step 8. Write the table and a summary
write.table(configs, OUTPUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
cat("distinct configurations:", nrow(configs), "from", nrow(settings), "settings ->", OUTPUT_FILE, "\n")
print(head(configs, 10), row.names = FALSE)
