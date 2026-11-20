#!/usr/bin/env python3
# =============================================================================
# v3_I10a_mean_sorted_hrs16.py — the HRS spike-in profile for the Fig 6d LOD rebuild: the mean of 16 sorted primary HRS
# samples of GSE301492 (author 2026-11-03: "yes drop sample 10 and 19, get mean from remaining 16 HRS samples as the
# input for spike in").
#
# Why 2 of the 18 samples are left out (screen v3_I9, 100 % spike-in in CD01; marker CPM vs the median of the other 16):
#   Sample10-HRS  HRS 0.057, fit r 0.64: CCL17 247 (15,575), MS4A1 1,236 (19), MKI67 1,001 (179); looks like proliferating B
#   Sample19-HRS  HRS 0.164, fit r 0.78, T_CD4 0.32: CD2 539 (9), CD3E 238 (10); looks contaminated with T cells
# Processing (the same as v3_I9a, so the 16 samples are treated as in the screen):
#   * Counts -> TPM: counts / gene length (kb), rescaled to 1e6. Gene length = union of the gene's exons in GENCODE
#     v24lift37 (GRCh37, as the GSE alignment); a symbol not in v24 takes its Ensembl 91 length. Genes with no length are
#     dropped (count printed). The 2 duplicated symbols are summed.
#   * Mean = plain mean of the 16 TPM profiles, gene by gene.
#   * Both gene names are kept: the GSE symbol and the name matched to the round-16 signature, because the mixtures must
#     follow the gene names of whichever clinical run they are compared with. Values are never changed by a rename.
#   * No rescaling to shared genes here; the step that builds the mixtures does that (as v1 2025-cHL_misc.Rmd).
# This script only writes the profile. It builds no mixtures and runs no CIBERSORTx.
# Output (OAK cHL_v2/output/r16/stepI/lod_rebuild_input/):
#   2026-11-03_mean_sorted_HRS16_TPM.tsv   gene_gse, gene_signature, mean_HRS16, then the 16 samples
# Run on Sherlock (in a job; the GTFs are large):  sbatch v3_I10a_mean_sorted_hrs16.sh
# =============================================================================
import csv
import gzip
import math
import os
import re

# Step 1. Files and settings
BRIAN = "${PROJECT_ROOT}"
R16_INPUT = BRIAN + "/cHL_v2/output/r16/stepD/input"              # round-16 sigmtx_merged.txt
OUT_DIR = BRIAN + "/cHL_v2/output/r16/stepI/lod_rebuild_input"
OUT_FILE = OUT_DIR + "/2026-11-03_mean_sorted_HRS16_TPM.tsv"
SIGNATURE_FILE = R16_INPUT + "/sigmtx_merged.txt"
GSE_COUNTS = BRIAN + "/cHL_v2/data/external/GSE301492/GSE301492_Roshal_Kong_readCounts.txt.gz"   # md5 964306b5...
EXCLUDED_SAMPLES = ["Sample10-HRS", "Sample19-HRS"]               # author's decision, reasons in the header
# gene annotations (the same four as v3_I9a)
TENX_GTF = "${ANNOTATION_DIR}/genes.gtf.gz"
V36_GTF = "${ANNOTATION_DIR}/gencode.v36.annotation.gtf"
E91_GTF = "${ANNOTATION_DIR}/Homo_sapiens.GRCh38.91.gtf"
V24_GTF = "${LAB_ROOT}/sequencing/pipeline/cfRNA_indices/gencode.v24lift37.annotation.plus.ERCC.gtf"

if os.path.exists(OUT_FILE):
    raise SystemExit(OUT_FILE + " already exists; refusing to overwrite")
os.makedirs(OUT_DIR, exist_ok=True)


# Step 2. Read a GTF: symbol -> Ensembl ID (first seen), Ensembl ID -> symbol, and symbol -> exon-union length (bp)
#         (unchanged from v3_I9a)
def read_gtf(path, with_lengths=False):
    symbol_to_ensembl = {}
    ensembl_to_symbol = {}
    exons_by_symbol = {}
    if path.endswith(".gz"):
        handle = gzip.open(path, "rt")
    else:
        handle = open(path)
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        if fields[2] != "gene" and not (with_lengths and fields[2] == "exon"):
            continue
        gene_id = re.search('gene_id "([^"]+)"', fields[8]).group(1)
        if gene_id.startswith("ENSGR") or "_PAR_Y" in gene_id:
            continue                                             # pseudo-autosomal Y copies
        name_match = re.search('gene_name "([^"]+)"', fields[8])
        ensembl = gene_id.split(".")[0]
        if name_match:
            symbol = name_match.group(1)
        else:
            symbol = ensembl
        if fields[2] == "gene":
            symbol_to_ensembl.setdefault(symbol, ensembl)
            ensembl_to_symbol[ensembl] = symbol
        else:
            exons_by_symbol.setdefault(symbol, []).append((int(fields[3]), int(fields[4])))
    handle.close()
    lengths = {}
    for symbol, exons in exons_by_symbol.items():
        exons.sort()
        total = 0
        current_start, current_end = exons[0]
        for start, end in exons[1:]:
            if start > current_end:
                total += current_end - current_start + 1
                current_start, current_end = start, end
            elif end > current_end:
                current_end = end
        total += current_end - current_start + 1
        lengths[symbol] = total
    return symbol_to_ensembl, ensembl_to_symbol, lengths


tenx_symbol_to_ensembl, _, _ = read_gtf(TENX_GTF)
v36_symbol_to_ensembl, _, _ = read_gtf(V36_GTF)
e91_symbol_to_ensembl, _, e91_lengths = read_gtf(E91_GTF, with_lengths=True)
v24_symbol_to_ensembl, _, v24_lengths = read_gtf(V24_GTF, with_lengths=True)


# signature_ensembl: Ensembl ID of a signature symbol (10x first, then GENCODE v36, then Ensembl 91) — as v3_I9a
def signature_ensembl(symbol):
    for lookup in (tenx_symbol_to_ensembl, v36_symbol_to_ensembl, e91_symbol_to_ensembl):
        if symbol in lookup:
            return lookup[symbol]
    return None


# older_ensembl: Ensembl ID of an older (GRCh37-era) symbol: v24 first, then Ensembl 91, then v36 — as v3_I9a
def older_ensembl(symbol):
    for lookup in (v24_symbol_to_ensembl, e91_symbol_to_ensembl, v36_symbol_to_ensembl):
        if symbol in lookup:
            return lookup[symbol]
    return None


# Step 3. Signature genes
signature_genes = []
with open(SIGNATURE_FILE) as handle:
    next(handle)
    for line in handle:
        signature_genes.append(line.split("\t", 1)[0].strip('"'))
signature_set = set(signature_genes)
signature_by_ensembl = {}
for symbol in signature_genes:
    ensembl = signature_ensembl(symbol)
    if ensembl is not None:
        signature_by_ensembl[ensembl] = symbol
print("signature genes:", len(signature_set))

# Step 4. Read the counts and sum duplicated symbols
with gzip.open(GSE_COUNTS, "rt") as handle:
    rows = list(csv.reader(handle, delimiter="\t"))
header = rows[0]
if header[0] == "":
    sample_names = header[1:]
else:
    sample_names = header
if len(sample_names) != len(rows[1]) - 1:
    raise SystemExit("GSE header / row width mismatch")
counts_by_symbol = {}
n_duplicates = 0
for row in rows[1:]:
    counts = [float(value) for value in row[1:]]
    if row[0] in counts_by_symbol:
        n_duplicates += 1
        counts_by_symbol[row[0]] = [a + b for a, b in zip(counts_by_symbol[row[0]], counts)]
    else:
        counts_by_symbol[row[0]] = counts

# Step 5. Which samples go into the mean: every sorted HRS sample except the two excluded ones (must leave 16)
all_hrs_samples = [name for name in sample_names if name.endswith("-HRS")]
kept_samples = [name for name in all_hrs_samples if name not in EXCLUDED_SAMPLES]
print("sorted HRS samples:", len(all_hrs_samples), "| excluded:", EXCLUDED_SAMPLES, "| kept:", len(kept_samples))
if len(all_hrs_samples) != 18 or len(kept_samples) != 16:
    raise SystemExit("expected 18 sorted HRS samples and 16 kept")

# Step 6. Gene lengths; a symbol with no length cannot be converted to TPM and is dropped
gene_lengths = {}
n_v24 = 0
n_e91 = 0
for symbol in counts_by_symbol:
    if symbol in v24_lengths:
        gene_lengths[symbol] = v24_lengths[symbol]
        n_v24 += 1
    elif symbol in e91_lengths:
        gene_lengths[symbol] = e91_lengths[symbol]
        n_e91 += 1
print(f"symbols {len(counts_by_symbol)} ({n_duplicates} duplicates summed); length from v24 {n_v24}, "
      f"from Ensembl 91 {n_e91}, none (dropped) {len(counts_by_symbol) - n_v24 - n_e91}")
symbols = [symbol for symbol in counts_by_symbol if symbol in gene_lengths]

# Step 7. Counts -> TPM for each of the 18 samples (the 18-sample mean is only used for the check in Step 10)
tpm_by_sample = {}
for sample_name in all_hrs_samples:
    column = sample_names.index(sample_name)
    rates = []
    for symbol in symbols:
        rates.append(counts_by_symbol[symbol][column] / (gene_lengths[symbol] / 1000))
    total_rate = sum(rates)
    tpm_by_sample[sample_name] = [rate / total_rate * 1e6 for rate in rates]

# Step 8. Mean of the 16 kept samples, and of all 18 for comparison
mean_16 = []
mean_18 = []
for i in range(len(symbols)):
    mean_16.append(sum(tpm_by_sample[name][i] for name in kept_samples) / len(kept_samples))
    mean_18.append(sum(tpm_by_sample[name][i] for name in all_hrs_samples) / len(all_hrs_samples))

# Step 9. Name matched to the round-16 signature: a symbol is renamed to the signature symbol of the same Ensembl gene,
#         if that symbol is not already a row (the v3_I9a rule; two old names never take the same signature symbol)
present_symbols = set(symbols)
used_targets = set()
signature_names = []
for symbol in symbols:
    new_name = symbol
    if symbol not in signature_set:
        target = signature_by_ensembl.get(older_ensembl(symbol))
        if target is not None and target not in present_symbols and target not in used_targets:
            new_name = target
            used_targets.add(target)
    signature_names.append(new_name)
if len(signature_names) != len(set(signature_names)):
    raise SystemExit("duplicate gene names after matching to the signature")
print(f"signature genes present: {len(signature_set & set(symbols))} under the GSE symbols, "
      f"{len(signature_set & set(signature_names))} after matching (of {len(signature_set)})")

# Step 10. Checks: the mean sums to 1e6, and how much the two excluded samples changed it (signature genes, log2 TPM + 1)
print("mean_HRS16 sum:", round(sum(mean_16), 3))
log_16 = []
log_18 = []
for i, name in enumerate(signature_names):
    if name in signature_set:
        log_16.append(math.log2(mean_16[i] + 1))
        log_18.append(math.log2(mean_18[i] + 1))
average_16 = sum(log_16) / len(log_16)
average_18 = sum(log_18) / len(log_18)
covariance = sum((a - average_16) * (b - average_18) for a, b in zip(log_16, log_18))
spread_16 = math.sqrt(sum((a - average_16) ** 2 for a in log_16))
spread_18 = math.sqrt(sum((b - average_18) ** 2 for b in log_18))
print("Pearson r, mean of 16 vs mean of 18, signature genes, log2(TPM + 1):", round(covariance / (spread_16 * spread_18), 4))
for symbol in ["CCL17", "TNFRSF8", "MS4A1", "CD3E", "CD2", "MKI67"]:
    if symbol in symbols:
        i = symbols.index(symbol)
        print(f"  {symbol}: mean of 16 = {mean_16[i]:.1f} TPM, mean of 18 = {mean_18[i]:.1f} TPM")

# Step 11. Write the profile
with open(OUT_FILE, "w") as out:
    out.write("\t".join(["gene_gse", "gene_signature", "mean_HRS16"] + kept_samples) + "\n")
    for i, symbol in enumerate(symbols):
        sample_values = [repr(tpm_by_sample[name][i]) for name in kept_samples]
        out.write("\t".join([symbol, signature_names[i], repr(mean_16[i])] + sample_values) + "\n")
print("written:", OUT_FILE)
