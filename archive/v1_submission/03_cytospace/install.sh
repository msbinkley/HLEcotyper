#!/bin/bash
# =============================================================================
# install.sh — one-time CytoSPACE installation (conda), run on the cluster.
# CytoSPACE: Vahid et al., Nat Biotechnol 2023; https://github.com/digitalcytometry/cytospace
# =============================================================================
CYTOSPACE_ROOT="${CYTOSPACE_ROOT:-/path/to/Cytospace}"
cd "$CYTOSPACE_ROOT"
git clone https://github.com/digitalcytometry/cytospace
cd cytospace
conda env create -f environment.yml          # creates env "cytospace"
conda activate cytospace
pip install .
pip install lapjv==1.3.14                    # lap_CSPR solver backend used in run_cytospace.sh
