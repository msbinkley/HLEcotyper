# Differential tumor-microbiome abundance — method & caveats

**Data:** step-05 `cross_validated_species.csv` (raw Kraken2/Bracken per-sample species calls,
depth-normalized RPM) for cHL (206) + NLPHL (192). Host (Homo sapiens) and kitome-blocklisted
taxa removed. Absent call = 0 RPM (structural/sampling zero, not interpreted as true absence).

**Data-integrity cleaning (added after first pass — duplicates + cross-cohort misclassification):**
Samples are filtered to `sample_design_clean.csv` `keep==TRUE`. Excluded (13 total):
- cHL run (206→199): 4 NLPHL-misclassified (`subtype_final`=NLPHL: cHL129/134/135/140) + 3 duplicate-
  patient extras (`is_multi_sample_patient`=TRUE; kept higher-depth: dropped cHL09/cHL20/cHL151).
- NLPHL run (192→186): 1 cHL-misclassified (`TPV_Subtype`="…cHL": LP143) + 5 duplicate extras
  (`dup_role`=duplicate, kept MAIN: dropped LP142/LP153/LP66/LP55/LP53). Join needed zero-pad
  normalization (microbiome `LP1` vs metadata `LP01`).
NLPHL EBV now from `TPV_EBV_Status_EBER` where available (1 kept NLPHL is EBER+, LP19; rest neg/unknown).
Residual unknowns: 17 cHL `subtype=NA` + 4 cHL unmatched to metadata could hide further misclassification.

**Grouping:** cHL EBV status from clinical metadata (`2026-06-21_cHL_metadata.csv`, join on `id`).
Pre-cleaning join: 202/206 matched, 76 EBV+ / 97 EBV− / 33 unknown. **After data-integrity cleaning
(below), the analysis arms are 75 EBV+ / 92 EBV− / 32 unknown.** NLPHL is treated as EBV− *by
clinical definition* (NLPHL is ~90–95% EBV-negative) and as a default where per-sample EBER is
missing — but this is a CAVEAT, not a validated fact: 1 kept NLPHL is EBER+ (LP19) and most kept
NLPHL have unknown EBER status.

**Comparisons (cleaned-design arm sizes):**
1. cHL (all 199) vs NLPHL (186) — CROSS-COHORT, batch-confounded → discovery only
2. EBV+ cHL (75) vs EBV− cHL (92) — INTERNAL, clean. **Primary.**
3. EBV− cHL (92) vs NLPHL (186) — cross-cohort; NLPHL = *mostly assumed/unknown* EBV− (not verified) → discovery only
4. EBV+ cHL (75) vs NLPHL (186) — cross-cohort, EBV-dominated → discovery only

**Test (codex-agreed):** per comparison, restrict to species with prevalence ≥10% in either arm
(EBV force-retained). Primary = zero-inclusive **Wilcoxon rank-sum on RPM**; companion = **Fisher
exact** on presence/absence. log2FC = log2((mean RPM_A + pc)/(mean RPM_B + pc)), pc = ½ a one-read
RPM (= 0.5·1e6/median total_reads) — used for effect size/plots **only**, never in the test.
Per-comparison Benjamini-Hochberg FDR. Species/taxid level primary (classifier + EBV control are
taxid-based); genus is the natural sensitivity level.

**Controls:**
- POSITIVE: EBV/HHV-4 (taxid 3050299 in step-05 species table; 10376 = HHV-4 S1 node). Must be up
  in EBV+ arms. Observed: EBV+ cHL prev 98.7% / EBV− cHL 49.5% / NLPHL 23.4%; log2FC 4.6 (24×) in
  the clean internal comparison, padj≈0. ✓ Supports EBV signal recovery + cHL metadata-join alignment.
  (Does NOT validate the "NLPHL is EBV−" assumption — NLPHL microbial-EBV prevalence is ~22–23% and
  1 kept NLPHL is clinically EBER+; NLPHL EBV status remains a caveat.)
- NEGATIVE / process: reagent & environmental contaminants — plant fungi (Pyricularia/rice blast,
  Zymoseptoria, Cercospora, Remersonia, Akanthomyces…), algae (Cyanidioschyzon, Hemiselmis = DB
  artifacts), kitome bacteria, and "viruses" like Betaretrovirus maspfimon (Mason-Pfizer monkey
  virus), Gorganvirus, Pahexavirus (Propionibacterium phages). Flagged via `contaminant` column;
  these should NOT be read as biology.

**Hard caveats:**
- Comparisons 1/3/4 are **fully confounded** by cohort = sequencing/processing batch (cHL_fastp vs
  NovaSeq_trimmed). No statistical adjustment is valid for a complete confound. EBV is the sole
  trustworthy cross-cohort signal because it is genuine tumor biology, not batch.
- Low-biomass FFPE → many zeros, low counts; non-EBV "hits" are dominated by contaminants/batch and
  should be treated as hypothesis-generating at best.

**Outputs:** `results/<comparison>_results.csv` (taxid, species, prevA/B, mean/med RPM, log2FC,
wilcox_p/padj, fisher_p/padj, contaminant, is_EBV); `results/figures/volcano_all_comparisons.pdf`;
`results/figures/EBV_positive_control.pdf` (written by `diff_abundance.R`). **Reproduce in order:**
`Rscript diff_abundance.R` (result CSVs + EBV control fig) → `Rscript plot_volcano_broken.R`
(the broken-axis `volcano_all_comparisons.pdf`; diff_abundance.R no longer writes it, to avoid
clobbering) → `Rscript ebv_concordance.R` (concordance metrics + figures).
