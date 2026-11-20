# =============================================================================
# v2_D14_r14_search_eval.R — score every run of the round-14 EcoTyper search (v2_D12) and pick the representative
# 2-ecotype and 3-ecotype models (author 2026-10-24: "do not worry about Steidl survival so far").
#
# Same scoring as the round-13 search (v2_D9_search_eval.R, ECOTYPER_SEARCH_PLAN_v2.md), minus Steidl:
#   - number of ecotypes; hard rule = no ecotype with two states of the same cell type
#   - collagen ecotype = the ecotype with the highest mean number of collagen/ECM genes per state among the top-50
#     genes (MaxFC) of its states; Th2 and IFN-gamma hits per state reported separately
#   - % of the 187 tumours assigned; ARI vs v1 over all 187 (Unassigned = its own class)
#   - max Guatemala share; robustness (n settings = number of cutoff x max-states settings giving the configuration)
#   - training Cox of the collagen-ecotype abundance (+ StageNum + Gender + age; duplicates removed; 10-year censoring;
#     5-year for the round-16 / 17 sets, author 2026-10-26):
#     REPORTED ONLY, not used to choose (author 2026-10-24: no survival in the selection for now)
#
# Representatives (the plan's selection rule without survival):
#   valid = 2 or 3 ecotypes + hard rule passed + a collagen ecotype (>= 1 collagen/ECM gene per state on average);
#   among valid runs with the same number of ecotypes: most settings first, then higher % assigned, then higher min
#   states (6 = the v1 value).
#
# Runs locally, on files pulled by v2_D13_r14_pull.sh into data_v2/02_intermediate/r14_search/min<M>/<tag>/.
# Usage:  Rscript v2_D14_r14_search_eval.R            (scores every pulled pass: min 6, 5, 4, 3)
# Output: figures_v2/EcoTyper_search/<date>_r14_search_ALL.tsv  and  <date>_r14_search_representatives.tsv
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival) })

# Step 1. Project folder and the shared v2 config (gives V2, p(), p2(), DUP_PATIENTS, CENSOR_YEARS)
MANUSCRIPT_DIR <- "${BOX_ROOT}/paper/Manuscript"
setwd(file.path(MANUSCRIPT_DIR, "script"))
source("../script_v2/figures/00_config_v2.R")

# Step 2. Passes (min states per ecotype) and the configuration table
#   tag of each configuration = c<representative cutoff>_m<representative max>_e<min>, as v2_D12 names its folders
#   (format() of the whole cutoff column, trailing zeros dropped)
#   Optional argument (author 2026-10-25): "dc" = the cDC + pDC pooled exploration search (r14_DC_search/; outputs
#   named r14_DC_search_*); default = the main round-14 search, as first run.
SET <- commandArgs(TRUE)[1]
if (is.na(SET)) SET <- "main"
MIN_PASSES <- c(6, 5, 4, 3)
config_name <- "r14_configs.tsv"
if (SET == "dc") {
  search_root <- file.path(V2, "r14_DC_search")
  output_label <- "r14_DC_search"
} else if (SET == "r16main") {           # round 16 (author 2026-10-26; v2_H1 / v2_H2), 13 types
  search_root <- file.path(V2, "r16_search")
  output_label <- "r16_search"
  config_name <- "r16_configs.tsv"
} else if (SET == "r16dc") {             # round 16, 12 types with cDC + pDC = DC
  search_root <- file.path(V2, "r16_DC_search")
  output_label <- "r16_DC_search"
  config_name <- "r16_configs.tsv"
} else if (SET == "r17main") {           # round 17 (author 2026-10-26; v2_H1_r17 / v2_H2_r17), 13 types
  search_root <- file.path(V2, "r17_search")
  output_label <- "r17_search"
  config_name <- "r17_configs.tsv"
} else if (SET == "r17dc") {             # round 17, 12 types with cDC + pDC = DC
  search_root <- file.path(V2, "r17_DC_search")
  output_label <- "r17_DC_search"
  config_name <- "r17_configs.tsv"
} else {
  search_root <- file.path(V2, "r14_search")
  output_label <- "r14_search"
}

# Survival endpoint (author 2026-10-26: "use 5-year PFS for the r16 and r17 search too"): the round-16 / 17 sets censor
# PFS at 5 years, as Fig 4 (fig4.R). The round-14 sets keep the config value (10 years) so their earlier tables reproduce.
if (SET %in% c("r16main", "r16dc", "r17main", "r17dc")) {
  CENSOR_YEARS <- 5
}
cat("PFS censored at", CENSOR_YEARS, "years\n")
configs <- fread(file.path(search_root, config_name))
cutoff_text <- format(configs$rep_cut, trim = TRUE, drop0trailing = TRUE)
MIN_ECM_PER_STATE <- 1            # "a collagen ecotype exists": at least 1 collagen/ECM gene per state on average

# Step 3. Gene lists (as v2_D9)
#   collagen / ECM genes: every COL<number> gene plus this list
ECM_GENES <- c("LUM", "DCN", "BGN", "SPARC", "THBS1", "THBS2", "MMP2", "MMP14", "LOX", "LOXL1", "LOXL2", "FN1", "POSTN", "CDH11",
               "PDGFRA", "PDGFRB", "FAP", "ACTA2", "TAGLN", "CALD1", "NID1", "NID2", "FBN1", "ELN", "TIMP3", "CTHRC1", "VCAN", "AEBP1",
               "EFEMP1", "EFEMP2", "MXRA5", "FSTL1", "CCDC80", "SFRP2", "SFRP4", "OLFML2B", "ADAM12", "PLOD2", "EMILIN1", "MFAP5",
               "SERPINF1", "TGFBI", "SNAI2", "ANGPTL2", "MMP11", "COMP", "SPON2", "PCOLCE", "SERPINH1", "P4HA2", "COLGALT1", "P3H1")
is_ecm <- function(genes) grepl("^COL[0-9]", genes) | genes %in% ECM_GENES      # TRUE for each collagen / ECM gene
TH2_GENES <- c("CCL17", "CCL22", "IL13RA2", "CD86", "IL5RA", "IL22RA2", "IL13", "CCR4")
IFN_GENES <- c("STAT1", "GBP1", "GBP2", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "IRF1", "WARS1", "UBD")

# Step 4. Helper functions (as v2_D9)
# adjusted Rand index of two labelings of the same tumours
adjusted_rand_index <- function(labels_a, labels_b) {
  counts <- table(labels_a, labels_b)
  pairs <- function(x) sum(choose(x, 2))
  expected <- pairs(rowSums(counts)) * pairs(colSums(counts)) / choose(sum(counts), 2)
  maximum  <- (pairs(rowSums(counts)) + pairs(colSums(counts))) / 2
  (pairs(counts) - expected) / (maximum - expected)
}

# ecotype per tumour, as a vector named by tumour ID
read_assignment <- function(assignment_file) {
  assignment <- read.delim(assignment_file, check.names = FALSE)
  setNames(as.character(assignment$Ecotype), assignment$ID)
}

# Cox model -> hazard ratio, 95% CI and p of the FIRST term (the collagen-ecotype abundance); all NA if the fit fails
cox_first_term <- function(formula, data) {
  fit_summary <- tryCatch(summary(coxph(formula, data = data)), error = function(e) NULL)
  if (is.null(fit_summary)) return(c(HR = NA, lo = NA, hi = NA, p = NA))
  c(HR = fit_summary$conf.int[1, 1],
    lo = fit_summary$conf.int[1, 3],
    hi = fit_summary$conf.int[1, 4],
    p  = fit_summary$coefficients[1, 5])
}

# Step 5. Inputs shared by all runs
# 5a. The 187 training tumours (round-14 13-type fractions, the EcoTyper input) and their v1 ecotype ("Unassigned" if none)
ids187 <- fread(p2("fractions13"))$Mixture
v1_assignment <- read_assignment(file.path(p("ecotyper_model_dir"), "Ecotypes", "ecotype_assignment.txt"))
v1_labels <- ifelse(ids187 %in% names(v1_assignment), v1_assignment[ids187], "Unassigned")

# 5b. Training clinical table: numeric age, duplicate patients removed, PFS in years censored at CENSOR_YEARS
clinical <- read.delim(p("clinical"))
clinical$age <- suppressWarnings(as.numeric(clinical$age))
clinical <- clinical[!clinical$ID %in% DUP_PATIENTS, ]
clinical$yrs <- clinical$PFS_Days / 365.25
clinical$ev  <- clinical$Progression_Censor
clinical$ev[clinical$yrs > CENSOR_YEARS] <- 0
clinical$yrs <- pmin(clinical$yrs, CENSOR_YEARS)

# 5c. Cohort of every training tumour (all tumours, duplicates included), named by ID
all_clinical <- read.delim(p("clinical"))
cohort <- setNames(all_clinical$Cohort, all_clinical$ID)

# Step 6. Score each run (one row per pass x configuration)
result_rows <- list()
for (min_states in MIN_PASSES) {
  search_dir <- file.path(search_root, paste0("min", min_states))
  if (!dir.exists(search_dir)) {
    cat("min", min_states, ": not pulled, skipped\n")
    next
  }
  for (i in seq_len(nrow(configs))) {
    tag <- sprintf("c%s_m%d_e%d", cutoff_text[i], configs$rep_max[i], min_states)
    config_dir <- file.path(search_dir, tag)
    base <- data.table(tag = tag, min = min_states, vec = configs$vec[i], n_settings = configs$n_settings[i],
                       max = configs$rep_max[i], cut = configs$rep_cut[i])

    # 6a. No ecotypes: say whether the run gave no ecotypes or was not run / not pulled
    ecotypes_file <- file.path(config_dir, "Ecotypes", "ecotypes.txt")
    if (!file.exists(ecotypes_file)) {
      if (dir.exists(config_dir)) {
        note <- "no ecotypes"
      } else {
        note <- "not run / not pulled"
      }
      result_rows[[length(result_rows) + 1]] <- cbind(base, n_eco = 0L, note = note)
      next
    }

    # 6b. Ecotype states; hard rule = no (ecotype, cell type) pair with more than one state
    ecotypes <- fread(ecotypes_file)
    n_eco <- uniqueN(ecotypes$Ecotype)
    states_per_type <- ecotypes[, .N, by = .(Ecotype, CellType)]
    duplicates <- states_per_type[N > 1]

    # 6c. Top-50 genes (highest MaxFC) of every ecotype state
    top_genes <- list()
    for (j in seq_len(nrow(ecotypes))) {
      gene_file  <- file.path(config_dir, ecotypes$CellType[j], "gene_info.txt")
      this_state <- ecotypes$State[j]
      if (!file.exists(gene_file)) {
        top_genes[[j]] <- character(0)
      } else {
        gene_info <- fread(gene_file)
        gene_info <- gene_info[gene_info$State == this_state, ]
        top_genes[[j]] <- head(gene_info$Gene[order(-gene_info$MaxFC)], 50)
      }
    }

    # 6d. Per ecotype (alphabetical): number of states, its states, and the mean number of ECM / Th2 / IFN genes per state
    per_ecotype_list <- list()
    for (ecotype in sort(unique(ecotypes$Ecotype))) {
      state_rows <- which(ecotypes$Ecotype == ecotype)
      ecm_counts <- c()
      th2_counts <- c()
      ifn_counts <- c()
      for (j in state_rows) {
        ecm_counts <- c(ecm_counts, sum(is_ecm(top_genes[[j]])))
        th2_counts <- c(th2_counts, sum(top_genes[[j]] %in% TH2_GENES))
        ifn_counts <- c(ifn_counts, sum(top_genes[[j]] %in% IFN_GENES))
      }
      per_ecotype_list[[length(per_ecotype_list) + 1]] <- data.table(
        Ecotype = ecotype, n_states = length(state_rows),
        states = paste(paste0(ecotypes$CellType[state_rows], "_", ecotypes$State[state_rows]), collapse = ","),
        ecm_per_state = mean(ecm_counts), th2_per_state = mean(th2_counts), ifn_per_state = mean(ifn_counts))
    }
    per_ecotype <- rbindlist(per_ecotype_list)

    # 6e. Collagen ecotype = highest mean ECM genes per state (first one if tied)
    collagen_eco <- per_ecotype[which.max(ecm_per_state), Ecotype]

    # 6f. Tumour assignment over the 187 ("Unassigned" if not listed)
    assignment <- read_assignment(file.path(config_dir, "Ecotypes", "ecotype_assignment.txt"))
    labels187 <- ifelse(ids187 %in% names(assignment), assignment[ids187], "Unassigned")

    # 6g. Largest share of Guatemala tumours in any ecotype
    guatemala_shares <- c()
    for (ecotype in unique(assignment)) {
      eco_ids <- names(assignment)[assignment == ecotype]
      guatemala_shares <- c(guatemala_shares, mean(cohort[eco_ids] == "Guatemala", na.rm = TRUE))
    }
    max_guatemala <- max(guatemala_shares)

    # 6h. Training Cox (reported only): collagen-ecotype abundance + StageNum + Gender + age
    abundance <- read.delim(file.path(config_dir, "Ecotypes", "ecotype_abundance.txt"), check.names = FALSE)
    collagen_abundance <- data.frame(ID = colnames(abundance), col = unlist(abundance[collagen_eco, ]))
    training <- merge(clinical, collagen_abundance, by = "ID")
    cox_training <- cox_first_term(Surv(yrs, ev) ~ col + StageNum + Gender + age, training)

    # 6i. Text columns
    if (nrow(duplicates)) {
      duplicates_text <- paste(duplicates$Ecotype, duplicates$CellType, collapse = ";")
    } else {
      duplicates_text <- ""
    }
    eco_summary <- paste(sprintf("%s:%d states, ECM %.1f, Th2 %.1f, IFN %.1f", per_ecotype$Ecotype, per_ecotype$n_states,
                                 per_ecotype$ecm_per_state, per_ecotype$th2_per_state, per_ecotype$ifn_per_state),
                         collapse = " | ")
    eco_states <- paste(paste0(per_ecotype$Ecotype, "=", per_ecotype$states), collapse = " | ")
    eco_sizes  <- paste(paste0(names(table(labels187)), ":", as.integer(table(labels187))), collapse = " ")

    # 6j. The run's row (HR / CI / p to 3 significant digits)
    result_rows[[length(result_rows) + 1]] <- cbind(base,
      n_eco                  = n_eco,
      hard_rule_ok           = nrow(duplicates) == 0,
      dup                    = duplicates_text,
      collagen_eco           = collagen_eco,
      collagen_ecm_per_state = round(per_ecotype[Ecotype == collagen_eco, ecm_per_state], 2),
      eco_summary            = eco_summary,
      eco_states             = eco_states,
      tumours_per_ecotype    = eco_sizes,
      pct_assigned           = round(100 * mean(labels187 != "Unassigned"), 1),
      ARI_v1_all187          = round(adjusted_rand_index(labels187, v1_labels), 3),
      max_guatemala_share    = round(max_guatemala, 2),
      train_HR = signif(cox_training[["HR"]], 3),
      train_lo = signif(cox_training[["lo"]], 3),
      train_hi = signif(cox_training[["hi"]], 3),
      train_p  = signif(cox_training[["p"]], 3),
      note     = "")
  }
}
results <- rbindlist(result_rows, fill = TRUE)

# Step 7. Valid runs and the representative 2- and 3-ecotype models
results[, valid := n_eco %in% 2:3 & hard_rule_ok %in% TRUE & collagen_ecm_per_state >= MIN_ECM_PER_STATE]
valid_runs <- results[valid == TRUE]
setorder(valid_runs, n_eco, -n_settings, -pct_assigned, -min)
valid_runs[, rank_within_n_eco := seq_len(.N), by = n_eco]
representatives <- valid_runs[rank_within_n_eco <= 5]          # the winner (rank 1) plus up to 4 runner-ups per n_eco

# Step 8. Write the date-prefixed tables and print a summary
out_dir <- file.path(MANUSCRIPT_DIR, "figures_v2", "EcoTyper_search")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(results, file.path(out_dir, sprintf("%s_%s_ALL.tsv", format(Sys.Date()), output_label)), sep = "\t")
fwrite(representatives, file.path(out_dir, sprintf("%s_%s_representatives.tsv", format(Sys.Date()), output_label)), sep = "\t")
for (min_states in unique(results$min)) {
  pass <- results[min == min_states]
  n_eco_counts <- table(pass$n_eco)
  cat(sprintf("min %d: %d runs; n_eco: %s; hard rule ok: %d; valid 2-3 HLE: %d\n", min_states, nrow(pass),
              paste(names(n_eco_counts), n_eco_counts, sep = "=", collapse = " "),
              sum(pass$hard_rule_ok %in% TRUE), sum(pass$valid)))
}
print(representatives[, .(n_eco, rank_within_n_eco, tag, n_settings, pct_assigned, ARI_v1_all187, max_guatemala_share,
                          collagen_eco, collagen_ecm_per_state, train_HR, train_p)])
