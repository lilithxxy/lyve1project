# 05-chapter6-donor-variability — thesis Chapter 6.2.3, Supplementary S7.x

Donor-to-donor variability in LYVE-1⁺ differentiation capacity. PBMCs from four
donors were hashtag-labelled and pooled (hash IDs `240106-Donor1-PBMC-L1-15`,
`240120-Donor1-PBMC-L1-64`, `240229-Donor1-PBMC-L1`, `240309-Donor2-PBMC-L1` →
`donor_15`, `donor_64`, `donor_6`, `donor_90`) and grouped by their LYVE-1⁺
yield (`efficacy` / `l1_capacity`): **low** = donor_15, donor_6; **high** =
donor_64, donor_90. The analysis is a Seurat (R) workflow with a CellTypist
(Python) QC/annotation notebook.

## Scripts (run order)

| # | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| 1 | `conversion.R` | Reads the per-donor `PBMC*.rds` Seurat objects (output of the demultiplexing step), drops the HTO assay, converts each to `.h5ad` via `sceasy` for the Python side | `~/Desktop/bioinfo/sc_R/unfiltered/PBMC*.rds` | `<sample>.h5ad` per donor |
| 2 | `momac.R` | Merges the four PBMC objects; QC `nFeature_RNA > 200 & < 8000`, `nCount_RNA < 80000`, `percent.mt < 5`; LogNormalize, 3000 HVGs, PCA, UMAP, Leiden-style clustering (res 0.5); donor → efficacy mapping; `FindMarkers(high vs low)`; pseudobulk `AggregateExpression` by efficacy × donor | `PBMC*.rds` | `integrated.rds` (save line commented), DE tables |
| 3 | `monos_.R` | Same merge with **SCTransform**, cell-cycle scoring, `percent.mt < 10`, res 0.4; adds the Day-0 MoMac sample (`day0M`); PBMC marker dictionary; `FindAllMarkers` top-10 heatmap; UMAPs by cluster/efficacy/phase; selected transporter-gene plots per donor | `~/Desktop/bioinfo/sc_R/unfiltered/toy/*.rds` | `monos_withm.rds`, `umap.png`, `ggsave1.png` |
| 4 | `integration_annotation.R` | Same lineage (SCTransform, `percent.mt < 10`, res 0.4) with manual cluster annotation → CD14⁺ Mono, CD16⁺ Mono, CD4⁺ T, CD8⁺ T, B, CD56±NK, DC, Platelet, Debris; clustered dot plots and heatmaps | `/home/yc/Desktop/monos/PBMC*.rds` | annotated object (in session), plots |
| 5 | `qc_celltypist.ipynb` | Python side: downloads CellTypist models, loads `Healthy_COVID19_PBMC.pkl`, reads the four `donor_*.h5ad`, labels `sample` and `l1_capacity` (low/high), computes scanpy QC metrics (mt/ribo/hb) and violin plots per donor | `donor_6/15/64/90.h5ad` | QC figures, CellTypist labels |

**Recovered file:** `qc_celltypist.ipynb` was copied from a Jupyter checkpoint
(`monos_final/pre_processing/.ipynb_checkpoints/qc-checkpoint.ipynb`,
2025-06-09). No live copy of this notebook existed anywhere; the checkpoint is
the only surviving version.

**Known issue at source:** `integration_annotation.R` calls
`SCTransform(integrated, verbose = FALSE, vars.to.regress = )` with an empty
`vars.to.regress` argument — this line errors as written and must have been
edited interactively when run.

## Not copied from this project

`monos_final/pbulk_ontology_cellchat.R` — not copied (its `#### cellchat`
section is an empty heading; the rest duplicates the pseudobulk/edgeR steps in
`momac.R`).

## Inputs (not in repo)

Per-donor Seurat objects produced by the HTODemux step
(`03-chapter4-scrna-timecourse/00_demultiplex_htodemux_full.R`), with obs
`hash.ID`, `HTO_classification`, `nCount_RNA`, `nFeature_RNA`.
