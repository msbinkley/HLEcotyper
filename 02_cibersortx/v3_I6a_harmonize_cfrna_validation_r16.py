#!/usr/bin/env python3
# =============================================================================
# v3_I6a_harmonize_cfrna_validation_r16.py — the cfRNA mixture of the CardioGuardian VALIDATION re-run: the 112 earlier
# cfRNA samples + the 29 NovaSeq B35 samples (CC69-CC97), with gene names matched to the ROUND-16 signature matrix
# (author 2026-11-02: "rerun deconvolution with the validation cohort, and compare cHL vs NLPHL vs other cancers vs
# healthy").
#
# Same logic as the 2026-09-25 B35 harmonization (Cibersort/2026-09-25_cHL_cfRNA_B35/gene_harmonization/
# harmonize_genes.py, made for the v1 signature), pointed at the round-16 signature:
#   * the cfRNA matrices use GENCODE v24 symbols; the round-16 signature uses the newer 10x GRCh38-2020-A symbols
#     (e.g. ADA2, CCN1, H1-0). A mixture row is renamed to the signature symbol when both are the same Ensembl gene.
#   * if a mixture row already carries a signature symbol but belongs to a DIFFERENT Ensembl gene, and the
#     signature's Ensembl gene is in the mixture under another v24 name, the row is set aside as
#     '<symbol>_v24_otherENSG' so the right gene can take the symbol.
#   * values are never changed; only row names.
# Not needed for round 16: the v1 step that restored R-mangled names ('-' -> '.') in the signature. Instead every
#   signature gene is checked to be in the round-16 reference.
#
# Run on Sherlock:  python3 v3_I6a_harmonize_cfrna_validation_r16.py
# Outputs (OAK cHL_v2/output/r16/stepI/cfrna_validation_input/):
#   cfRNA_TPM_141_harmonized_r16.txt   genes in both matrices x 141 samples (112 earlier, then 29 B35)
#   mixture_rename_log.tsv             every renamed / set-aside row
# =============================================================================
import csv
import gzip
import os
import re

# Step 1. Files
BRIAN = "${PROJECT_ROOT}"
R16_INPUT = BRIAN + "/cHL_v2/output/r16/stepD/input"              # round-16 sigmtx_merged.txt and ref_merged.txt
B35_INPUTS = BRIAN + "/Cibersort/2026-09-25_cHL_cfRNA_B35/inputs"  # the two raw cfRNA TPM matrices
OUT_DIR = BRIAN + "/cHL_v2/output/r16/stepI/cfrna_validation_input"
SIGNATURE_FILE = R16_INPUT + "/sigmtx_merged.txt"
REFERENCE_FILE = R16_INPUT + "/ref_merged.txt"
EARLIER_TPM = B35_INPUTS + "/2025-09-23_cfRNA_coding_symbol_TPM.txt"     # 112 samples (the v3 Fig 6b run's mixture)
B35_TPM = B35_INPUTS + "/2026-09-25_NovaSeqB35_genes_symbol_TPM.txt"      # 29 samples, columns Sample_CCxx_cfrna
# gene annotations (the same four as the 2026-09-25 harmonization)
TENX_GTF = "${ANNOTATION_DIR}/genes.gtf.gz"
V36_GTF = "${ANNOTATION_DIR}/gencode.v36.annotation.gtf"
E91_GTF = "${ANNOTATION_DIR}/Homo_sapiens.GRCh38.91.gtf"
V24_GTF = "${LAB_ROOT}/sequencing/pipeline/cfRNA_indices/gencode.v24lift37.annotation.plus.ERCC.gtf"

if os.path.exists(OUT_DIR + "/cfRNA_TPM_141_harmonized_r16.txt"):
    raise SystemExit(OUT_DIR + " already has the mixture; refusing to overwrite")
os.makedirs(OUT_DIR, exist_ok=True)


# Step 2. Read a GTF into two lookups: symbol -> Ensembl ID (first one seen) and Ensembl ID -> symbol
def read_gtf_genes(path):
    symbol_to_ensembl = {}
    ensembl_to_symbol = {}
    if path.endswith(".gz"):
        handle = gzip.open(path, "rt")
    else:
        handle = open(path)
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "gene":
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
        symbol_to_ensembl.setdefault(symbol, ensembl)
        ensembl_to_symbol[ensembl] = symbol
    handle.close()
    return symbol_to_ensembl, ensembl_to_symbol


tenx_symbol_to_ensembl, _ = read_gtf_genes(TENX_GTF)
v36_symbol_to_ensembl, _ = read_gtf_genes(V36_GTF)
e91_symbol_to_ensembl, _ = read_gtf_genes(E91_GTF)
v24_symbol_to_ensembl, v24_ensembl_to_symbol = read_gtf_genes(V24_GTF)


# signature_ensembl: Ensembl ID of a signature symbol (10x first, then GENCODE v36, then Ensembl 91)
def signature_ensembl(symbol):
    for lookup in (tenx_symbol_to_ensembl, v36_symbol_to_ensembl, e91_symbol_to_ensembl):
        if symbol in lookup:
            return lookup[symbol]
    return None


# Step 3. Signature genes; each must be in the round-16 reference (first column of ref_merged.txt)
signature_genes = []
with open(SIGNATURE_FILE) as handle:
    next(handle)
    for line in handle:
        signature_genes.append(line.split("\t", 1)[0].strip('"'))
signature_set = set(signature_genes)
reference_genes = set()
with open(REFERENCE_FILE) as handle:
    next(handle)
    for line in handle:
        reference_genes.add(line.split("\t", 1)[0].strip('"'))
not_in_reference = signature_set - reference_genes
print("signature genes:", len(signature_set), "| not in the round-16 reference:", len(not_in_reference))
if not_in_reference:
    raise SystemExit("signature genes missing from the reference: " + ", ".join(sorted(not_in_reference)[:20]))

# Step 4. Rename rules: v24 symbol -> signature symbol, when both are the same Ensembl gene.
#   A v24 symbol that maps to two signature symbols is ambiguous and is not renamed.
candidates = {}
for symbol in signature_genes:
    ensembl = signature_ensembl(symbol)
    if ensembl is None:
        continue
    v24_symbol = v24_ensembl_to_symbol.get(ensembl)
    if v24_symbol is not None and v24_symbol != symbol:
        candidates.setdefault(v24_symbol, []).append(symbol)
rename_rules = {}
ambiguous = {}
for v24_symbol, targets in candidates.items():
    if len(targets) == 1:
        rename_rules[v24_symbol] = targets[0]
    else:
        ambiguous[v24_symbol] = targets
print("rename rules:", len(rename_rules), "| ambiguous (not renamed):", ambiguous)


# Step 5. Harmonize one mixture's row names (see the header for the two rules)
def belongs_to_other_gene(name, present_names):
    # True when 'name' is a signature symbol, but the signature's Ensembl gene is in the mixture under another
    # v24 name and 'name' itself is a different Ensembl gene in v24 (a bare ID change keeps the symbol match)
    ensembl = signature_ensembl(name)
    if ensembl is None:
        return False
    other_name = v24_ensembl_to_symbol.get(ensembl)
    if other_name is None or other_name == name or other_name not in present_names:
        return False
    return v24_symbol_to_ensembl.get(name) != ensembl


def harmonize(rows, tag):
    names = [row[0] for row in rows]
    present_names = set(names)
    set_aside = set()
    for name in names:
        if name in signature_set and belongs_to_other_gene(name, present_names):
            set_aside.add(name)
    names_kept = present_names - set_aside
    new_names = []
    log = []
    for name in names:
        if name in set_aside:
            new_name = name + "_v24_otherENSG"
            log.append((name, new_name, "collision_moved_aside", tag))
        elif name in rename_rules and rename_rules[name] not in names_kept:
            new_name = rename_rules[name]
            log.append((name, new_name, "renamed_to_signature_symbol", tag))
        else:
            new_name = name
        new_names.append(new_name)
    if len(new_names) != len(set(new_names)):
        raise SystemExit(tag + ": duplicate gene names after harmonization")
    new_rows = []
    for new_name, row in zip(new_names, rows):
        new_rows.append([new_name] + row[1:])
    n_renamed = sum(1 for entry in log if entry[2] == "renamed_to_signature_symbol")
    n_aside = sum(1 for entry in log if entry[2] == "collision_moved_aside")
    print(f"{tag}: renamed {n_renamed}, set aside {n_aside}; signature genes present "
          f"{len(signature_set & set(new_names))}/{len(signature_set)} (before: {len(signature_set & present_names)})")
    return new_rows, log


def read_matrix(path):
    with open(path) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0], rows[1:]


earlier_header, earlier_rows = read_matrix(EARLIER_TPM)
b35_header, b35_rows = read_matrix(B35_TPM)
b35_samples = []
for column in b35_header[1:]:
    b35_samples.append(re.sub(r"^Sample_(CC\d+)_cfrna$", r"\1", column))   # Sample_CC69_cfrna -> CC69
print("samples: earlier", len(earlier_header) - 1, "| B35", len(b35_samples))

earlier_rows, earlier_log = harmonize(earlier_rows, "earlier cfRNA (112)")
b35_rows, b35_log = harmonize(b35_rows, "B35 (29)")

# Step 6. Combine: genes in both matrices, in the order of the earlier matrix
b35_by_gene = {}
for row in b35_rows:
    b35_by_gene[row[0]] = row[1:]
combined_header = ["Gene"] + earlier_header[1:] + b35_samples
combined_rows = []
for row in earlier_rows:
    if row[0] in b35_by_gene:
        combined_rows.append(row + b35_by_gene[row[0]])
combined_genes = set(row[0] for row in combined_rows)
print("combined:", len(combined_rows), "genes x", len(combined_header) - 1, "samples; signature genes present",
      len(signature_set & combined_genes), "/", len(signature_set))

with open(OUT_DIR + "/cfRNA_TPM_141_harmonized_r16.txt", "w") as out:
    out.write("\t".join(combined_header) + "\n")
    for row in combined_rows:
        out.write("\t".join(row) + "\n")
with open(OUT_DIR + "/mixture_rename_log.tsv", "w") as out:
    out.write("v24_symbol\tnew_name\taction\tfile\n")
    for entry in earlier_log + b35_log:
        out.write("\t".join(entry) + "\n")
print("written to", OUT_DIR)
