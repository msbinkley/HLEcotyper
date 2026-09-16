# 3. CytoSPACE — single-cell to Visium spot assignment

CytoSPACE v1.1 (https://github.com/digitalcytometry/cytospace) maps cells from the annotated
cHL scRNA/snRNA-seq reference onto each Visium spot of the four spatially profiled tumors
(Fig. 5a; Supplementary Fig.). Solver: `lap_CSPR` (`lapjv`).

```bash
bash install.sh              # once: conda env "cytospace" + lapjv
# build inputs (R, Seurat):
Rscript -e 'source("generate_cytospace_from_seurat_object.R");
            obj <- readRDS("cHL_annotated.rds"); Idents(obj) <- obj$pooled_celltype;
            generate_cytospace_from_scRNA_seurat_object(obj, dir_out="cHL_scRNA_Input")'
Rscript generate_cytospace_input_from_spaceranger_output.R /path/to/spaceranger/<tumor>/outs cHL_Visium_Input/<tumor>
bash run_cytospace.sh        # submits one SLURM job per tumor
```

| file | produced by | content |
|---|---|---|
| `cHL_scRNA_Input/scRNA_data.txt` | `generate_cytospace_from_seurat_object.R` | genes × cells raw counts (`GENES` first column) |
| `cHL_scRNA_Input/cell_type_labels.txt` | same | `Cell IDs`, `CellType` (13 pooled labels) |
| `cHL_Visium_Input/ST_data_<tumor>.csv` | `generate_cytospace_input_from_spaceranger_output.R` | genes × spots counts |
| `cHL_Visium_Input/Coordinates_<tumor>.csv` | same | `Spot ID`, row, col |
| `output/<tumor>_results/` | `run_cytospace.sh` | `assigned_locations.csv`, `cell_type_assignments_by_spot.csv`, `fractional_abundances_by_spot.csv`, plots |

The two R helpers are the unmodified utilities distributed with CytoSPACE.
`fractional_abundances_by_spot.csv`, renormalised so each spot sums to 1, is the
`Recovery cell type fractions` input for the `*_Cyto` variant of EcoTyper Visium recovery
(`02_ecotyper/06_recovery_visium.sh`); spatial autocorrelation of the resulting cell-type
and cell-state maps was summarised per tumor with Moran's I (Seurat `RunMoransI`).
