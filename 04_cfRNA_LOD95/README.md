# 4. cfRNA LOD95 — limit of detection of HRS cells in plasma cell-free RNA

Simulation-based LOD (Fig. 6c–d): healthy-volunteer cfRNA profiles are spiked in silico with
HRS cell-line transcripts, deconvoluted with the cHL signature matrix, and the detection
probability is modelled as a function of spike-in fraction.

```bash
# 1. build the simulated mixtures (R)
Rscript 01_make_spikein_mixtures.R cfRNA_coding_symbol_TPM.txt LL100_TPM.txt LOD95
#    -> LOD95_HRS_background_spikein_TPM.txt  (copy to $INPUT_DIR)   + LOD95_design.txt
# 2. deconvolute (CIBERSORTx HiRes S-mode, SLURM)
bash 02_cibersortx_hires_spikein.sh
# 3. threshold + logistic fit
Rscript 03_fit_LOD95.R $CIBERSORT_ROOT/cfRNA_LOD95_spikein/CIBERSORTxGEP_NA_Fractions-Adjusted.txt LOD95_design.txt LOD95
```

| item | value |
|---|---|
| cfRNA input | healthy-volunteer cfRNA TPM (STAR + RSEM), protein-coding genes; volunteers `CD01`–`CD20` |
| HRS source | L-428 cell line RNA-seq TPM from the LL-100 blood-cancer cell-line panel (Quentmeier et al. 2019; ArrayExpress E-MTAB-7721) |
| background set | volunteers CD11–CD20 × 10 noisy replicates, no spike-in |
| spike-in set | volunteers CD01–CD10 × 16 fractions (100 % → 0.001 %) × 10 replicates |
| noise | Gaussian, sd = 0.1 % of the profile mean, floored at 0; columns renormalised to 1e6 |
| detection threshold | mean + 3 SD of background HRS fraction |
| model | `glm(Detected ~ log10(SpikeIn), family = binomial)` |
| result | LOD95 = 1.17 % (manuscript) |

Both simulation steps use `set.seed(1)`; the replicate counts are constants at the top of
`01_make_spikein_mixtures.R`.
