#!/usr/bin/env Rscript
# =============================================================================
# 05_aggregate_filter.R
# Merge Kraken2/Bracken species calls, remove host & contaminants, apply
# abundance/minimizer/prevalence filters, and (optionally) run decontam against
# negative controls.
#
# Called by 05_downstream.sbatch, which passes config values as positional args:
#  1 RESULTS_DIR          2 BRACKEN_BIOM(species)   3 MINIMIZER_STATS
#  4 TOTAL_READS_LOOKUP   5 READCOUNTS_DIR          6 BLOCKLIST
#  7 MIN_READS            8 MIN_RPM                 9 MIN_PREVALENCE_FRAC
# 10 MIN_MINIMIZER_COVERAGE 11 MIN_DISTINCT_MINIMIZERS  12 CONTROLS_LIST (path or "NONE")
# =============================================================================
suppressPackageStartupMessages({
  # explicit core packages (not the tidyverse meta-package)
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(purrr)
  library(phyloseq); library(biomformat)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 12) stop(sprintf("Expected exactly 12 arguments, got %d; see header.", length(a)))
RESULTS_DIR <- a[1]; BIOM <- a[2]; MMZ <- a[3]; TOTALS <- a[4]
READCOUNTS <- a[5]; BLOCKLIST <- a[6]
MIN_READS <- as.numeric(a[7]); MIN_RPM <- as.numeric(a[8])
MIN_PREV <- as.numeric(a[9]); MIN_COV <- as.numeric(a[10])
MIN_DISTINCT <- as.numeric(a[11]); CONTROLS <- a[12]

# Fail-closed safety toggles (env-var overrides; default = safe/fail-closed). These
# guard against silently bypassing a gate when a step-03 input is missing/partial.
#   ALLOW_MINIMIZER_FAILOPEN: taxa with no minimizer row auto-PASS the distinct-
#       minimizer gate (legacy behavior). Default FALSE = such taxa FAIL the gate,
#       and a missing minimizer_stats.tsv is fatal instead of skipping the gate.
#   ALLOW_RPM_FALLBACK: when a sample has no total_reads denominator, substitute
#       classified_reads (legacy). Default FALSE = abort, since it silently changes
#       what MIN_RPM means.
istrue <- function(s) tolower(trimws(s)) %in% c("true","1","yes","t")
ALLOW_MINIMIZER_FAILOPEN <- istrue(Sys.getenv("ALLOW_MINIMIZER_FAILOPEN", "false"))
ALLOW_RPM_FALLBACK       <- istrue(Sys.getenv("ALLOW_RPM_FALLBACK", "false"))

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
log <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---- helpers ----------------------------------------------------------------
clean_sample <- function(x) {
  # strip extensions outside-in so e.g. "S1.species_bracken.kreport" -> "S1"
  x %>% basename() %>%
    str_remove("\\.kreport\\.minimizer$") %>%
    str_remove("\\.kreport$") %>%
    str_remove("\\.species_bracken$") %>% str_remove("\\.genus_bracken$") %>%
    str_remove("_bracken$") %>% str_remove("\\.bracken$") %>%
    str_remove("\\.species$") %>% str_remove("\\.genus$") %>%
    str_remove("_microbe$") %>%
    str_trim()
}
clean_tax <- function(x) x %>% str_remove_all("^[a-z]__") %>%
  str_replace_all(";", "") %>% str_trim()
norm_name <- function(x) x %>% tolower() %>% str_trim()
# Genus-qualified species key so bare epithets (e.g. "pneumoniae") cannot
# false-match across genera (Klebsiella vs Streptococcus pneumoniae). NA-safe.
sp_key_of  <- function(g, s) {
  g <- ifelse(is.na(g), "", g); s <- ifelse(is.na(s), "", s)
  norm_name(paste(g, s))
}
# Human-readable label: prepend genus UNLESS the species string already begins with
# it. This genus-qualifies bare epithets AND genus-less placeholder names (e.g.
# "sp. (in: enterobacteria)"), so the label cannot collide ACROSS genera (Serratia
# vs Sodalis), while true binomials ("Staphylococcus aureus") are left intact. NA-safe.
sp_display <- function(g, s) {
  g <- ifelse(is.na(g), "", g); s <- ifelse(is.na(s), "", s)
  starts_with_genus <- nzchar(g) & str_starts(norm_name(s), stringr::fixed(norm_name(g)))
  ifelse(starts_with_genus, str_trim(s), str_trim(paste(g, s)))
}

# ---- 1. Kraken2/Bracken (species) from BIOM ---------------------------------
log("Loading Bracken BIOM:", BIOM)
kb <- tryCatch({
  ps <- import_biom(BIOM)
  m  <- psmelt(ps)
  rk <- paste0("Rank", 1:7)
  m %>%
    rename(reads = Abundance, sample = Sample) %>%
    mutate(taxid = as.character(OTU),
           superkingdom = clean_tax(.data[[rk[1]]]),
           phylum = clean_tax(.data[[rk[2]]]), class = clean_tax(.data[[rk[3]]]),
           order = clean_tax(.data[[rk[4]]]), family = clean_tax(.data[[rk[5]]]),
           genus = clean_tax(.data[[rk[6]]]), species = clean_tax(.data[[rk[7]]]),
           sample = clean_sample(sample)) %>%
    filter(reads > 0) %>%
    select(sample, taxid, superkingdom, phylum, class, order, family, genus, species, reads)
}, error = function(e) { log("WARN: could not read BIOM:", conditionMessage(e)); NULL })

if (is.null(kb) || nrow(kb) == 0) stop("No Kraken2/Bracken data parsed; cannot continue.")

# ---- 2. Total reads (host-input) for RPM ------------------------------------
totals <- tryCatch(
  read_tsv(TOTALS, col_names = c("sample","total_reads"), show_col_types = FALSE) %>%
    mutate(sample = clean_sample(sample)),
  error = function(e) tibble(sample = character(), total_reads = numeric()))

kb <- kb %>% left_join(totals, by = "sample") %>%
  group_by(sample) %>% mutate(classified_reads = sum(reads)) %>% ungroup()
# RPM is defined against the FULL input library (total_reads_lookup.tsv). If a
# sample has no/zero denominator, RPM would silently become classified-read-
# normalized and MIN_RPM would mean something different -> abort by default.
missing_tot <- kb %>% filter(is.na(total_reads) | total_reads <= 0) %>%
  distinct(sample) %>% pull(sample)
if (length(missing_tot) > 0) {
  if (ALLOW_RPM_FALLBACK) {
    log("WARN:", length(missing_tot), "sample(s) lack a positive total_reads denominator;",
        "falling back to classified_reads for RPM (ALLOW_RPM_FALLBACK=true). Samples:",
        paste(head(missing_tot, 10), collapse = ","))
    kb <- kb %>% mutate(total_reads = ifelse(is.na(total_reads) | total_reads <= 0,
                                             classified_reads, total_reads))
  } else {
    stop(sprintf(paste0("Missing/zero total_reads denominator for %d sample(s); RPM/MIN_RPM would be wrong. ",
                        "Provide a complete total_reads_lookup.tsv, or set ALLOW_RPM_FALLBACK=true to use ",
                        "classified_reads. Samples: %s"),
                 length(missing_tot), paste(head(missing_tot, 20), collapse = ",")))
  }
}
kb <- kb %>% mutate(RPM = reads / total_reads * 1e6,
                    percent = reads / classified_reads * 100)

# ---- 3. Distinct-minimizer (KrakenUniq-style) FP filter ---------------------
if (file.exists(MMZ)) {
  mmz <- read_tsv(MMZ, show_col_types = FALSE) %>%
    mutate(sample = clean_sample(sample), taxid = as.character(taxid),
           coverage = ifelse(total_minimizers > 0, distinct_minimizers/total_minimizers, 0))
  kb <- kb %>% left_join(select(mmz, sample, taxid, distinct_minimizers, coverage),
                         by = c("sample","taxid"))
  kb <- kb %>% mutate(
    distinct_minimizers = ifelse(is.na(distinct_minimizers), NA_real_, distinct_minimizers),
    # "fail-open" = no minimizer row at all (e.g. Bracken redistributed reads to a
    # species the raw Kraken2 report never assigned directly, so it carries no
    # distinct-minimizer evidence). FAIL-CLOSED by default: such a taxon does NOT
    # pass the KrakenUniq-style gate, so the FP filter can never be silently
    # bypassed. Set ALLOW_MINIMIZER_FAILOPEN=true to restore the legacy auto-pass.
    minimizer_failopen = is.na(distinct_minimizers),
    pass_minimizer = dplyr::coalesce(distinct_minimizers >= MIN_DISTINCT & coverage >= MIN_COV, FALSE) |
      (ALLOW_MINIMIZER_FAILOPEN & minimizer_failopen))
  n_na_mmz <- sum(is.na(kb$distinct_minimizers))
  log("Applied minimizer filter (distinct>=", MIN_DISTINCT, ", coverage>=", MIN_COV,
      "; fail-open=", ALLOW_MINIMIZER_FAILOPEN, ")")
  if (n_na_mmz > 0)
    log("NOTE:", n_na_mmz, "taxon-calls had no minimizer row ->",
        if (ALLOW_MINIMIZER_FAILOPEN) "auto-PASSED (ALLOW_MINIMIZER_FAILOPEN=true);"
        else "FAILED the gate (fail-closed default);",
        "verify step 02 wrote a *.kreport.minimizer for every sample.")
} else {
  if (!ALLOW_MINIMIZER_FAILOPEN)
    stop("minimizer_stats.tsv not found (", MMZ, "); the distinct-minimizer FP gate ",
         "would be skipped for EVERY taxon. Re-run step 03 so it is produced, or set ",
         "ALLOW_MINIMIZER_FAILOPEN=true to intentionally skip the gate.")
  kb <- kb %>% mutate(distinct_minimizers = NA_real_, coverage = NA_real_,
                      minimizer_failopen = TRUE, pass_minimizer = TRUE)
  log("WARN: no minimizer_stats.tsv; ALLOW_MINIMIZER_FAILOPEN=true -> skipping KrakenUniq-style filter (all taxa pass).")
}

# ---- 4. Host + contaminant (kitome) flags -----------------------------------
block <- character(0)
if (file.exists(BLOCKLIST)) {
  block <- read_lines(BLOCKLIST) %>% str_remove("#.*$") %>% str_trim()
  block <- norm_name(block[block != ""])
}
is_host <- function(g, s) norm_name(g) == "homo" | norm_name(s) == "homo sapiens"
kb <- kb %>% mutate(
  # coalesce: a row with NA genus/species would make `==` return NA, which then
  # poisons high_confidence, the rejected complement, and the funnel sum()s
  # (NA stage counts). Treat unknown taxonomy as not-host (fail-open to "keep").
  is_host = dplyr::coalesce(is_host(genus, species), FALSE),
  is_blocklisted = norm_name(genus) %in% block | norm_name(species) %in% block,
  pass_abundance = reads >= MIN_READS & RPM >= MIN_RPM)

# ---- 5. Build the species-keyed table ---------------------------------------
# Genus-qualified species key + human-readable label; drop genus-only rows with an
# empty species, exactly as before.
xval <- kb %>%
  mutate(sp_key = sp_key_of(genus, species), ge_key = norm_name(genus),
         species_label = sp_display(genus, species)) %>%
  filter(!is.na(species), norm_name(species) != "", norm_name(species) != "na")

# prevalence across samples
prev <- xval %>% filter(pass_abundance, !is_host, !is_blocklisted, pass_minimizer) %>%
  group_by(sp_key) %>% summarise(n_samples = n_distinct(sample), .groups = "drop")
# Cohort denominator = every sample that entered classification (one
# readcount.tsv per host-removed sample), NOT just samples that happened to
# yield a surviving Bracken call -- otherwise prevalence is inflated and every
# "X / Y samples" figure reports the wrong Y.
cohort_n <- if (dir.exists(READCOUNTS))
  length(list.files(READCOUNTS, pattern = "\\.readcount\\.tsv$")) else 0L
# When readcounts/ is absent, fall back to the total_reads_lookup universe (one row
# per classified sample) rather than n_distinct(xval$sample): the latter is taken
# AFTER reads>0 filtering, so zero-call samples would drop out and INFLATE prevalence.
totals_n <- dplyr::n_distinct(totals$sample)
n_total  <- max(cohort_n, totals_n, n_distinct(xval$sample))
if (n_total <= 0)
  stop("Cannot establish a cohort denominator for prevalence (no readcounts/ dir and ",
       "empty total_reads_lookup.tsv). Provide one so MIN_PREVALENCE_FRAC is computed ",
       "against the true sample count.")
xval <- xval %>% left_join(prev, by = "sp_key") %>%
  mutate(n_samples = replace_na(n_samples, 0),
         pass_prevalence = (n_samples / pmax(n_total,1)) >= MIN_PREV)

# final high-confidence flag
xval <- xval %>% mutate(
  high_confidence = pass_abundance & pass_minimizer & pass_prevalence &
    !is_host & !is_blocklisted)

# ---- 6. Optional decontam against negative controls -------------------------
decontam_tab <- NULL
contam_ids   <- character(0)   # decontam-flagged taxids; applied to high_confidence below
if (CONTROLS != "NONE" && file.exists(CONTROLS) && requireNamespace("decontam", quietly = TRUE)) {
  log("Running decontam (prevalence) against controls in:", CONTROLS)
  ctl <- read_tsv(CONTROLS, col_names = c("sample","type"), show_col_types = FALSE) %>%
    mutate(sample = clean_sample(sample))
  # Prevalence matrix from the Kraken2/Bracken read counts.
  mat <- xval %>% group_by(sample, taxid) %>% summarise(reads = sum(reads), .groups="drop") %>%
    pivot_wider(names_from = sample, values_from = reads, values_fill = 0) %>%
    column_to_rownames("taxid") %>% as.matrix()
  smp <- intersect(colnames(mat), ctl$sample)
  if (length(smp) > 2 && any(ctl$type[match(smp, ctl$sample)] == "control")) {
    is_ctrl <- ctl$type[match(smp, ctl$sample)] == "control"
    dc <- decontam::isContaminant(t(mat[, smp, drop=FALSE]), neg = is_ctrl, method = "prevalence")
    decontam_tab <- dc %>% rownames_to_column("taxid")
    contam_ids <- decontam_tab %>% filter(contaminant) %>% pull(taxid)
    write_csv(decontam_tab, file.path(RESULTS_DIR, "decontam_contaminants.csv"))
    log("decontam flagged", length(contam_ids), "contaminant taxa.")
  } else log("WARN: controls list lacks usable control samples; skipping decontam.")
} else if (CONTROLS != "NONE") {
  log("NOTE: decontam not run (no controls file or 'decontam' package missing).")
}

# Apply the decontam exclusion. No-op when decontam did not run (contam_ids is
# empty -> taxid %in% character(0) is all FALSE).
xval <- xval %>% mutate(is_decontam_contaminant = taxid %in% contam_ids,
                        high_confidence = high_confidence & !is_decontam_contaminant)

log("Kraken2/Bracken sample-taxon calls:", nrow(xval))

# ---- 7. Write outputs -------------------------------------------------------
out_cols <- c("sample","taxid","superkingdom","phylum","class","order","family",
              "genus","species","species_label","reads","total_reads","RPM","percent",
              "distinct_minimizers","coverage","minimizer_failopen","n_samples",
              "is_host","is_blocklisted","pass_minimizer","pass_abundance",
              "pass_prevalence","high_confidence")
out_cols <- intersect(out_cols, names(xval))
write_csv(xval %>% select(all_of(out_cols)),
          file.path(RESULTS_DIR, "cross_validated_species.csv"))

hc <- xval %>% filter(high_confidence)
write_csv(hc %>% select(all_of(out_cols)),
          file.path(RESULTS_DIR, "high_confidence_taxa.csv"))

# Calls that did NOT pass every high-confidence gate (the complement of
# high_confidence within cross_validated). `fail_reason` names which gate(s) each
# call failed, so the table is filterable. Gzipped because this audit table can be large.
rejected <- xval %>% filter(!high_confidence) %>%
  mutate(fail_reason = str_remove(paste0(
    ifelse(is_host, "host;", ""),
    ifelse(is_blocklisted, "kitome_contaminant;", ""),
    ifelse(!dplyr::coalesce(pass_abundance,  FALSE), "low_abundance;", ""),
    ifelse(!dplyr::coalesce(pass_minimizer,  TRUE),  "low_distinct_minimizers;", ""),
    ifelse(!dplyr::coalesce(pass_prevalence, TRUE),  "low_prevalence;", ""),
    ifelse( dplyr::coalesce(is_decontam_contaminant, FALSE), "decontam_contaminant;", "")),
    ";$"))
write_csv(rejected %>% select(any_of(c(out_cols, "fail_reason"))),
          file.path(RESULTS_DIR, "rejected_calls.csv.gz"))
log("Wrote", nrow(rejected), "calls that failed >=1 gate to rejected_calls.csv.gz")

# Group on species_label (genus + de-doubled species), which is consistent whether a
# source stored the bare epithet or the full binomial -- so the same organism is not
# split into two cohort rows.
cohort <- hc %>% group_by(superkingdom, genus, species_label) %>%
  summarise(species = dplyr::first(species),
            n_samples = n_distinct(sample),
            mean_RPM   = mean(RPM),
            median_RPM = stats::median(RPM),
            max_RPM    = max(RPM),
            .groups = "drop") %>%
  relocate(species, .after = genus) %>%
  arrange(desc(n_samples), desc(mean_RPM))
write_csv(cohort, file.path(RESULTS_DIR, "cohort_high_confidence_summary.csv"))

# --- Filter funnel: how many sample-taxon calls survive each successive gate ---
g_host  <- !xval$is_host
g_block <- g_host  & !xval$is_blocklisted
g_minz  <- g_block & xval$pass_minimizer
g_abund <- g_minz  & xval$pass_abundance
g_prev  <- g_abund & xval$pass_prevalence
# How often did the distinct-minimizer guard actually bite among the FINAL calls?
# A high fail-open share means most high-confidence species had no minimizer row
# (Bracken-redistributed) and the KrakenUniq-style gate never screened them --
# worth knowing before trusting it as an FP control.
n_hc          <- nrow(hc)
n_hc_failopen <- sum(hc$minimizer_failopen, na.rm = TRUE)
log(sprintf("Minimizer gate among high-confidence: %d/%d on real distinct-minimizer evidence, %d fail-open (no minimizer row).",
            n_hc - n_hc_failopen, n_hc, n_hc_failopen))
funnel <- tibble(
  stage = c("0. raw sample-taxon calls (Kraken2/Bracken)",
            "1. after host (Homo sapiens) removal",
            "2. after kitome blocklist",
            "3. after distinct-minimizer filter",
            "4. after abundance filter (reads/RPM)",
            "5. after prevalence filter",
            "6. HIGH-CONFIDENCE (final)"),
  calls = c(nrow(xval), sum(g_host), sum(g_block), sum(g_minz),
            sum(g_abund), sum(g_prev), nrow(hc)),
  distinct_species = c(n_distinct(xval$sp_key), n_distinct(xval$sp_key[g_host]),
            n_distinct(xval$sp_key[g_block]), n_distinct(xval$sp_key[g_minz]),
            n_distinct(xval$sp_key[g_abund]), n_distinct(xval$sp_key[g_prev]),
            n_distinct(hc$sp_key)))
funnel <- funnel %>% mutate(calls_dropped_here = dplyr::lag(calls) - calls)
write_csv(funnel, file.path(RESULTS_DIR, "filter_funnel_summary.csv"))

# --- QC read accounting (per-stage removals + Kraken2 classified/unclassified) ---
qc <- tryCatch({
  files <- list.files(READCOUNTS, pattern="\\.readcount\\.tsv$", full.names=TRUE)
  if (length(files)==0) NULL else map_dfr(files, ~read_tsv(.x, show_col_types = FALSE))
}, error=function(e) NULL)
if (!is.null(qc)) {
  if (!"after_t2t" %in% names(qc)) qc$after_t2t <- qc$after_dedup   # back-compat
  qc <- qc %>% mutate(
    host_removed_STAR         = input_pairs - star_unmapped,
    rRNA_removed              = star_unmapped - after_rrna,
    lowcomplexity_removed     = after_rrna - after_entropy,
    duplicates_removed        = after_entropy - after_dedup,
    t2t_removed               = after_dedup - after_t2t,
    pct_host_STAR             = round(100*(input_pairs - star_unmapped)/pmax(input_pairs,1), 2),
    pct_rRNA_removed          = round(100*(star_unmapped - after_rrna)/pmax(star_unmapped,1), 2),
    pct_lowcomplexity_removed = round(100*(after_rrna - after_entropy)/pmax(after_rrna,1), 2),
    pct_duplicates_removed    = round(100*(after_entropy - after_dedup)/pmax(after_entropy,1), 2),
    pct_t2t_removed           = round(100*(after_dedup - after_t2t)/pmax(after_dedup,1), 2),
    pct_microbial_final       = round(100*final_microbial/pmax(input_pairs,1), 4))
  # fold in Kraken2 per-sample classified/unclassified (written by step 03)
  kqc_path <- file.path(dirname(TOTALS), "kraken2_classification_QC.tsv")
  if (file.exists(kqc_path)) {
    kqc <- tryCatch(read_tsv(kqc_path, show_col_types = FALSE) %>%
                      mutate(sample = clean_sample(sample)), error = function(e) NULL)
    if (!is.null(kqc)) qc <- qc %>% left_join(kqc, by = "sample")
  }
  # per-sample count of high-confidence species (believable microbes per sample)
  qc <- qc %>%
    left_join(hc %>% count(sample, name = "n_high_confidence_species"), by = "sample") %>%
    mutate(n_high_confidence_species = tidyr::replace_na(n_high_confidence_species, 0L))
  write_csv(qc, file.path(RESULTS_DIR, "QC_read_accounting.csv"))
}

# --- One-glance narrative summary (pipeline_summary.txt) ---------------------
med <- function(x) if (length(x)) round(stats::median(suppressWarnings(as.numeric(x)), na.rm = TRUE), 2) else NA
n_samp <- n_total   # cohort size (samples that entered classification), matches the prevalence denominator
con <- file(file.path(RESULTS_DIR, "pipeline_summary.txt"), "w")
wl <- function(...) writeLines(paste0(...), con)
wl("Tumor-microbiome pipeline - run summary")
wl("=======================================")
wl("Samples analyzed: ", n_samp)
wl("")
wl("READ FUNNEL  (median %, each stage relative to the reads ENTERING that stage):")
if (exists("qc") && !is.null(qc)) {
  wl(sprintf("  host removed (STAR) ........ %s%%", med(qc$pct_host_STAR)))
  wl(sprintf("  rRNA removed ............... %s%%", med(qc$pct_rRNA_removed)))
  wl(sprintf("  low-complexity removed ..... %s%%", med(qc$pct_lowcomplexity_removed)))
  wl(sprintf("  duplicates removed ......... %s%%", med(qc$pct_duplicates_removed)))
  wl(sprintf("  T2T residual host removed .. %s%%", med(qc$pct_t2t_removed)))
  wl(sprintf("  -> final candidate-microbial %s%% of input", med(qc$pct_microbial_final)))
  if ("pct_unclassified" %in% names(qc)) wl(sprintf("  Kraken2 unclassified (median) %s%%", med(qc$pct_unclassified)))
} else wl("  (no read-accounting files found)")
wl("")
wl("TAXON-CALL FUNNEL  (Kraken2/Bracken; one call = one species in one sample):")
for (i in seq_len(nrow(funnel)))
  wl(sprintf("  %-46s calls=%-7d species=%d", funnel$stage[i], funnel$calls[i], funnel$distinct_species[i]))
wl(sprintf("  minimizer gate among high-confidence: %d/%d on real distinct-minimizer evidence, %d fail-open (Bracken-redistributed, no minimizer row)",
           n_hc - n_hc_failopen, n_hc, n_hc_failopen))
wl("")
wl("TOP HIGH-CONFIDENCE SPECIES (by prevalence):")
top <- utils::head(cohort, 15)
if (nrow(top) > 0) for (i in seq_len(nrow(top)))
  wl(sprintf("  %2d. %-34s %d/%d samples  median RPM %.1f",
             i, top$species_label[i], top$n_samples[i], n_samp, top$median_RPM[i]))
wl("")
wl("Key files: high_confidence_taxa.csv | cohort_high_confidence_summary.csv |")
wl("           cross_validated_species.csv | QC_read_accounting.csv | filter_funnel_summary.csv |")
wl("           rejected_calls.csv.gz (failed >=1 gate, with fail_reason)")
close(con)

log("Done. Wrote results to:", RESULTS_DIR)
log("Wrote pipeline_summary.txt (human-readable run summary).")
log(sprintf("High-confidence: %d sample-taxon calls, %d distinct species.",
            nrow(hc), n_distinct(hc$species)))
