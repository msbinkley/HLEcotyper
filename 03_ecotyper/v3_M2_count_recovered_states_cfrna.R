# =============================================================================
# v3_M2_count_recovered_states_cfrna.R — how many of the 30 cell states of the final model are detected in the plasma
# cfRNA recovery (v3_M1_cfrna_recovery_12samples.sbatch). Rule of the manuscript: a cell state is detected when at
# least one plasma sample is assigned to it.
#
# Inputs:  <model folder>     the final model (resources/ecotyper_model: <CellType>/gene_info.txt)
#          <recovery folder>  the recovery output (<CellType>/state_assignment.txt)
# Output:  printed: the samples of every detected state, the number of detected states, the states not detected
# How to run:  Rscript --vanilla v3_M2_count_recovered_states_cfrna.R <model folder> <recovery folder>
# Randomness: none.
# =============================================================================
arguments <- commandArgs(trailingOnly = TRUE)
MODEL_DIR <- arguments[1]
RECOVERY_DIR <- arguments[2]
CELL_TYPES <- c("B_preGC", "B_postGC", "CAF", "DC_Conventional", "DC_Plasmacytoid", "HRS", "Macrophage",
                "NK", "Plasma", "T_CD4", "T_CD8", "T_gd", "T_Treg")

# Step 1. All states of the model (from each cell type's gene_info.txt)
model_states <- character(0)
for (cell_type in CELL_TYPES) {
  gene_info <- read.delim(file.path(MODEL_DIR, cell_type, "gene_info.txt"), check.names = FALSE)
  model_states <- c(model_states, paste(cell_type, sort(unique(gene_info$State))))
}
cat("states in the model:", length(model_states), "\n")

# Step 2. The states with at least one assigned sample, and the samples of each state
detected <- character(0)
for (cell_type in CELL_TYPES) {
  assignment <- read.delim(file.path(RECOVERY_DIR, cell_type, "state_assignment.txt"), check.names = FALSE)
  for (state in sort(unique(assignment$State))) {
    detected <- c(detected, paste(cell_type, state))
    samples <- assignment$ID[assignment$State == state]
    cat(sprintf("  %-22s %d samples: %s\n", paste(cell_type, state), length(samples), paste(samples, collapse = " ")))
  }
}

# Step 3. The count (run of the paper: 28 of 30; not detected: DC_Conventional S01, DC_Conventional S02)
cat(sprintf("%d of %d states detected; not detected: %s\n", length(detected), length(model_states),
            paste(setdiff(model_states, detected), collapse = ", ")))
