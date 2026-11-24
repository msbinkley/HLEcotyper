# =============================================================================
# v2_D_undot_sigmatrix.py — CIBERSORTx writes signature-matrix gene names through R make.names
# (HLA-A -> HLA.A, 7SK -> X7SK). Map every name back to the reference's own gene symbol, so it matches the
# bulk / pseudo-bulk mixtures.
#
# Usage: python3 v2_D_undot_sigmatrix.py <refsample.txt> <sigmatrix.txt> <out.txt>
#   (called by job A of v2_D1_cibersortx.sh)
# Output: the signature matrix with the gene names restored; every other column unchanged.
# Exit code 1 if any signature gene could not be matched to the reference.
# =============================================================================
import re
import sys

# Step 1. File names from the command line
ref_file, sig_file, out_file = sys.argv[1:4]

# Step 2. Reference gene symbols = first column of the reference, header line skipped
reference_genes = [line.split("\t", 1)[0] for line in open(ref_file)][1:]


# Step 3. R make.names, as CIBERSORTx applies it:
#   every character other than a letter, digit, "." or "_" becomes ".",
#   and "X" is put in front if the name starts with a digit, "_", or "." followed by a digit
def make_names(gene):
    name = re.sub(r"[^A-Za-z0-9._]", ".", gene)
    if re.match(r"[0-9_]|\.[0-9]", name):
        return "X" + name
    return name


# Step 4. Lookup table: R name -> reference symbol (if two symbols give the same R name, the first one wins)
r_name_to_symbol = {}
for gene in reference_genes:
    r_name_to_symbol.setdefault(make_names(gene), gene)

# Step 5. Rewrite the signature matrix line by line: header unchanged; for every gene line,
#   keep the name if it is already a reference symbol, otherwise translate it back, otherwise keep it and count it
sig_lines = open(sig_file).read().splitlines()
n_changed = 0
n_missing = 0
with open(out_file, "w") as out:
    out.write(sig_lines[0] + "\n")
    for line in sig_lines[1:]:
        gene, rest = line.split("\t", 1)
        if gene in reference_genes:
            new_name = gene
        elif gene in r_name_to_symbol:
            new_name = r_name_to_symbol[gene]
            n_changed += 1
        else:
            new_name = gene
            n_missing += 1
        out.write(new_name + "\t" + rest + "\n")

# Step 6. Report; fail (exit code 1) if any gene is not in the reference
print(f"{len(sig_lines) - 1} signature genes; renamed back {n_changed}; not in reference {n_missing}")
if n_missing:
    sys.exit(1)
