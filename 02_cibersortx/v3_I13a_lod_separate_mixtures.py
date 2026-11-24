#!/usr/bin/env python3
# =============================================================================
# v3_I13a_lod_separate_mixtures.py — the two mixtures of the Fig 6d LOD with background and spike-in in SEPARATE S-mode runs.
# Author 2026-11-05: "for each 47 non-cancer samples, generate 4 gausian replicates. two replicate goes to background, the
# other two got spiked in at various concentrations (add back 2 to 1,5,10). run s mode on background group and spike in
# group seperately."
#
# Contents:
#   background file:  47 donors x replicates 1, 2, no spike-in         =    94 samples   Frac_0_Sample<donor>_Rep_<1 or 2>
#   spiked file:      47 donors x replicates 3, 4 x 16 levels          = 1,504 samples   Frac_<f>_Sample<donor>_Rep_<3 or 4>
# Choices:
#   * Every donor is in both groups, each group with its own replicates. Why not one group per donor: author's design.
#   * Levels: the 1-2-5-10 series, 100 % down to 0.001 %, the 16 levels of v1. "add back 2" restores the 20 %, 0.2 %, 0.02 % and
#     0.002 % dropped on 2026-11-03.
#   * Spike-in profile: mean of the 16 sorted primary HRS samples of GSE301492 (v3_I10a).
#   * Donors: the 47 non-cancer samples of the v3_I6a matrix (gene names matched to the round-16 signature).
#   * Shared genes = v3_I6a matrix n HRS profile, the same genes as v3_I10b. Every profile is rescaled to 1e6 on them (as v1).
#   * Noise size as v1: Gaussian, mean 0, SD = 0.001 x mean TPM, drawn per gene and per replicate, negative values set to 0.
#     Seeded here (v1 had no seed).
#   * Order, as the author described it: the 4 replicates are made first, then replicates 3 and 4 are spiked:
#     mixed = f x HRS + (1 - f) x noisy donor replicate. So a replicate carries the same noise through all levels. (v1 and
#     v3_I10b drew new noise after mixing.) One result: at f = 1 the donor term is 0, so all 94 columns at 100 % are identical.
# Outputs (OAK cHL_v2/output/r16/stepI/lod_separate_input/):
#   2026-11-05_lod_background.txt    18,8xx genes x 94 samples
#   2026-11-05_lod_spiked.txt        18,8xx genes x 1,504 samples
# Run on the Mac:   python3 v3_I13a_lod_separate_mixtures.py <local folder holding cfrna_validation_input/ and lod_rebuild_input/>
#                   then copy lod_separate_input/ to OAK stepI/. (Or on Sherlock:  sbatch v3_I13a_lod_separate_mixtures.sh)
# =============================================================================
import array
import os
import random
import sys

# Step 1. Files and settings. An optional argument replaces STEP_DIR, to build on the Mac from local copies of the two inputs
#         (author 2026-11-05: "you can run it locally"), keeping the same folder names inside it.
STEP_DIR = "${PROJECT_ROOT}/cHL_v2/output/r16/stepI"
if len(sys.argv) > 1:
    STEP_DIR = sys.argv[1]
REAL_FILE = STEP_DIR + "/cfrna_validation_input/cfRNA_TPM_141_harmonized_r16.txt"
PROFILE_FILE = STEP_DIR + "/lod_rebuild_input/2026-11-03_mean_sorted_HRS16_TPM.tsv"
OUT_DIR = STEP_DIR + "/lod_separate_input"
OUT_BACKGROUND = OUT_DIR + "/2026-11-05_lod_background.txt"
OUT_SPIKED = OUT_DIR + "/2026-11-05_lod_spiked.txt"
NOISE_SEED = 6                      # the project's seed number; makes the noise reproducible
NOISE_FRACTION = 0.001              # noise SD = this x mean TPM (v1)
BACKGROUND_REPLICATES = [1, 2]      # author: two replicates go to background
SPIKED_REPLICATES = [3, 4]          # author: the other two are spiked
SPIKE_LEVELS = [1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001]
N_DONORS = 47

if os.path.exists(OUT_DIR):
    raise SystemExit(OUT_DIR + " already exists; refusing to overwrite")
os.makedirs(OUT_DIR)

# Step 2. Read the HRS profile, under the gene names matched to the signature
hrs_by_gene = {}
with open(PROFILE_FILE) as handle:
    header = handle.readline().rstrip("\n").split("\t")
    if header[1] != "gene_signature" or header[2] != "mean_HRS16":
        raise SystemExit("unexpected columns in the HRS profile")
    for line in handle:
        fields = line.split("\t", 3)
        hrs_by_gene[fields[1]] = float(fields[2])

# Step 3. Read the non-cancer donors (names start with CB or CD), keeping only genes also in the HRS profile
with open(REAL_FILE) as handle:
    header = handle.readline().rstrip("\n").split("\t")
    real_samples = header[1:]
    donors = []
    for sample in real_samples:
        if sample.startswith("CB") or sample.startswith("CD"):
            donors.append(sample)
    donors.sort()
    donor_columns = {}
    for donor in donors:
        donor_columns[donor] = array.array("d")
    shared_genes = []
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if fields[0] not in hrs_by_gene:
            continue
        shared_genes.append(fields[0])
        for sample, value in zip(real_samples, fields[1:]):
            if sample in donor_columns:
                donor_columns[sample].append(float(value))
if len(donors) != N_DONORS or len(shared_genes) != len(set(shared_genes)):
    raise SystemExit("expected 47 donors and unique gene names")
n_genes = len(shared_genes)
print("donors:", len(donors), "| shared genes:", n_genes)

# Step 4. Rescale every profile to 1e6 on the shared genes (as v1)
for donor in donors:
    total = sum(donor_columns[donor])
    rescaled = array.array("d")
    for value in donor_columns[donor]:
        rescaled.append(value / total * 1e6)
    donor_columns[donor] = rescaled
hrs_values = [hrs_by_gene[gene] for gene in shared_genes]
hrs_total = sum(hrs_values)
hrs_profile = array.array("d", [value / hrs_total * 1e6 for value in hrs_values])

# Step 5. Make the 4 noisy replicates of each donor. The mean TPM is 1e6 / n_genes after rescaling, so the noise SD is the
#         same for every donor.
noise_sd = NOISE_FRACTION * (1e6 / n_genes)
print("noise SD (TPM):", round(noise_sd, 4))
noise_generator = random.Random(NOISE_SEED)
replicates = {}
for donor in donors:
    for replicate_number in BACKGROUND_REPLICATES + SPIKED_REPLICATES:
        noisy = array.array("d")
        for value in donor_columns[donor]:
            noisy.append(max(0.0, value + noise_generator.gauss(0, noise_sd)))
        replicates[(donor, replicate_number)] = noisy

# Step 6. Background group: replicates 1 and 2, unchanged
background_names = []
background_columns = []
for donor in donors:
    for replicate_number in BACKGROUND_REPLICATES:
        background_names.append(f"Frac_0_Sample{donor}_Rep_{replicate_number}")
        background_columns.append(replicates[(donor, replicate_number)])

# Step 7. Spiked group: replicates 3 and 4 at every level, mixed = f x HRS + (1 - f) x noisy replicate
spiked_names = []
spiked_columns = []
for fraction in SPIKE_LEVELS:
    for donor in donors:
        for replicate_number in SPIKED_REPLICATES:
            noisy = replicates[(donor, replicate_number)]
            mixed = array.array("d")
            for i in range(n_genes):
                mixed.append(fraction * hrs_profile[i] + (1 - fraction) * noisy[i])
            spiked_names.append(f"Frac_{fraction}_Sample{donor}_Rep_{replicate_number}")
            spiked_columns.append(mixed)

# Step 8. Checks before writing: counts, unique names, no negatives, sums near 1e6
if len(background_names) != 94 or len(spiked_names) != 1504:
    raise SystemExit("wrong number of samples")
if len(set(background_names + spiked_names)) != 94 + 1504:
    raise SystemExit("duplicate sample names")
for label, columns in (("background", background_columns), ("spiked", spiked_columns)):
    column_sums = [sum(column) for column in columns]
    smallest_value = min(min(column) for column in columns)
    print(label, "samples:", len(columns), "| column sums: min", round(min(column_sums), 1), "max", round(max(column_sums), 1),
          "| smallest value:", smallest_value)
    if smallest_value < 0:
        raise SystemExit("negative value in " + label)


# Step 9. Write one mixture file
def write_mixture(path, names, columns):
    with open(path, "w") as out:
        out.write("\t".join(["Gene"] + names) + "\n")
        for i, gene in enumerate(shared_genes):
            values = [repr(column[i]) for column in columns]
            out.write(gene + "\t" + "\t".join(values) + "\n")
    print("written:", path)


write_mixture(OUT_BACKGROUND, background_names, background_columns)
write_mixture(OUT_SPIKED, spiked_names, spiked_columns)
