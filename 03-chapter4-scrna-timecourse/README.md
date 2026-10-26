# 03-chapter4-scrna-timecourse — thesis Chapter 4

Time-resolved scRNA-seq of human monocyte-derived macrophages differentiated
with MDI vs MDIX, Days 1–5 (one donor, hashtag-pooled per timepoint;
Figure 4.4.1). Clusters `leiden_res_0_8` correspond to thesis Mac1–Mac7.

## Main scripts (run order)

| # | Script | What it does | Inputs | Outputs |
|---|---|---|---|---|
| 00 | `00_demultiplex_htodemux_full.R` | **Entry point of the chapter (R/Seurat).** Hashtag demultiplexing: `NormalizeData(assay="HTO", normalization.method="CLR")` → `HTODemux(assay="HTO", positive.quantile = 0.99)` → keep `HTO_classification.global == "Singlet"` & not `Negative` → split by hashtag → per-hashtag `.rds`. Then Seurat-side QC (`nFeature_RNA 500–8000, nCount_RNA < 50000, percent.mt < 15`; second sample `nFeature_RNA 300–8000, nCount_RNA < 60000, percent.mt < 25`), merge of day objects with `timepoint`/`condition`, CCA integration, marker plots, `AddModuleScore` of the bulk LYVE-1⁺/⁻ top-500 signatures (`l1pos`/`l1neg`; Fig 4.4.12 forward projection), cell-cycle scoring. Interactive session script — several objects are created outside the file. | cellranger `outs/filtered_feature_bc_matrix` per run (Gene Expression + Antibody Capture) under `~/Desktop/bioinfo/singlecell_human/` | `<hash.ID>.rds` per sample, merged Seurat object |
| 00b | `00b_demultiplex_htodemux_single.R` | Compact single-run version of the same demultiplexing (same CLR + `positive.quantile = 0.99`, Singlet filter), merges per-day `.rds` with `timepoint`, converts to `final.h5ad` via `sceasy` | `R1_DMo_D0M/outs/filtered_feature_bc_matrix` | `final.h5ad` |
| 01 | `01_day1_qc_cluster.ipynb` | Per-day processing (see note below): QC plots, **mt filter `< 30%`**, cell-cycle scoring, normalise/HVG/PCA, neighbours, UMAP, Leiden (res 1 and 0.6), PAGA, diffusion map, apoptosis scoring, marker DE, pre-annotation | `~/Desktop/scanpy/3_1/3_1_day1/3_1_day1.h5ad` | `day1final.h5ad`, `~/Desktop/scanpy/4/day1_30.h5ad`, DE CSVs, QC/UMAP figures |
| 02 | `02_day2_qc_cluster.ipynb` | same, Day 2, **mt `< 30%`** | `.../3_1_day2/3_1_day2.h5ad` | `day2_30.h5ad`, DE CSVs |
| 03 | `03_day3_qc_cluster.ipynb` | same, Day 3, **mt `< 30%`** (+ LYVE1hi/lo comparison) | `.../3_1_day3/...` | `day3_30.h5ad`, DE CSVs |
| 04 | `04_day4_qc_cluster.ipynb` | same, Day 4, **mt `< 40%`** (+ LYVE1hi/lo comparison) | `.../3_1_day4/...` | `day4_40.h5ad`, DE CSVs |
| 05 | `05_day5_qc_cluster.ipynb` | same, Day 5, **mt `< 50%`** (+ LYVE1hi/lo comparison) | `.../3_1_day5/...` | `day5_50.h5ad`, DE CSVs |
| 06 | `06_integration_benchmark.ipynb` | Benchmarks integration across the five per-day h5ads: BBKNN, Harmony, Seurat CCA, STACAS, uncorrected (+ scIB/iLISI metrics) | `~/Desktop/scanpy/4/*.h5ad` | inline plots only |
| 07 | `07_cell_counts.py` | Cell numbers before/after QC per day × condition (Figure 4.4.4); counts hard-coded from the QC notebooks | none | `cell_counts_qc_facet.png` |
| 08 | `08_scenic_regulons.ipynb` | SCENIC MAF-regulon analysis, 3 figures (thesis Figures 4.5.1–4.5.3) | `inte_scenic.loom`, `inte_palantir.h5ad`, `reg_day1..4.csv` | `fig1a/1b_MAF_*`, `fig2_MAF_network_day4`, `fig3_MAF_target_trajectory` |
| 09 | `09_temporal_gsea.R` | clusterProfiler GO/GSEA/KEGG loop over per-cluster DE tables | `sc_deg.csv` (see below) | `sc_deg/`, `sc_go/`, `sc_gsea/` CSVs + RDS, gseaplot2 figures |
| 10 | `10_slingshot_pseudotime.R` | **Slingshot pseudotime (R)** on the UMAP embedding of the Seurat object (`slingshot()`, `getLineages`, start cluster set manually), extracts `slingPseudotime_1/2/3`, beeswarm plot of cells ordered by pseudotime (one of the three methods in Fig 4.4.20) | `scored_final.rds` / the MDIX Seurat object | plots (in session) |
| 11 | `11_monocle3_pseudotime.R` | **Monocle3 pseudotime (R)**: Seurat preprocessing of `final.rds` (PCA/clusters res 0.6/UMAP) → `SeuratWrappers::as.cell_data_set` → `cluster_cells`, `learn_graph(use_partition = F)`, `order_cells`, pseudotime plots; `AddModuleScore` LYVE-1/MHC scores; manual cluster naming (immature / intermediate / LYVE1hi / postLYVE1hi / MDI / MDIX) | `final.rds` | plots (in session) |
| — | `exploratory/celloracle_pax8_exploratory.ipynb` | Exploratory CellOracle run on the condition-level GRN — **PAX8 knockout only, not reported in the thesis** (moved out of the main pipeline). | `inte_palantir.h5ad` | figures (in session) |

**Naming note (misnamed at source):** `01–05` were originally `qc_dayN_XX.ipynb`.
They are **not** QC-only — each contains neighbours, Leiden clustering,
diffusion maps and gene scoring. The trailing number in the original filename
(`30/30/30/40/50`) encodes the mitochondrial-percent threshold applied in that
notebook, listed above.

**Integration decision (deliberate):** the final Chapter-4 analysis uses **no
batch correction** — per-day objects are concatenated and processed jointly
(see `momac/uncorrected_new.ipynb`). The benchmark in `06` informed this
choice; note its metric-scoring cells were not fully executed.

**Ambient RNA / doublets (in-vitro data):** `00_demultiplex_htodemux_full.R`
loads `library(decontX)` under a "decontamination" comment but **never calls
`decontX()`**; doublet removal relies on the HTODemux `Singlet` classification
(hashtag doublets) only — no transcriptome-based doublet tool is run on this
dataset (`sc.pp.scrublet` is commented out in `momac/scvelo.ipynb`).

### Input formats

- `sc_deg.csv` (for `09`): concatenated per-cluster DE results in 5-column
  chunks — `<cluster>_names, <cluster>_scores, <cluster>_pvals,
  <cluster>_pvals_adj, <cluster>_logfoldchanges` for each cluster.
- Per-day h5ads carry obs columns `condition` (`MDI` / `MDIX`), `timepoint`,
  and (upstream of `momac/scvelo.ipynb`) a hashtag-assignment column `hash.ID`.

## `momac/` — final in-vitro analysis chain

Ran on a Linux server (`/home/...` paths); reads the per-day data and produces
the final integrated object (`inte_palantir.h5ad`) and most Chapter-4 figures.

| Script | Role |
|---|---|
| `scvelo.ipynb` | Entry point: reads upstream QC'd h5ads, merges spliced/unspliced looms, assigns `condition` (`MDIX`/`MDI`/`Momac`) from hashtag IDs, QC metrics (mt/ribo/hb), filters `min_genes=100`, `min_cells=3` (mito filter present but commented out), writes the MDIX-subset h5ads |
| `uncorrected_new.ipynb` | Concatenates per-day h5ads, HVG/PCA/neighbours/UMAP/Leiden **without batch correction**, Palantir pseudotime; writes `integrated.h5ad` / `inte_palantir.h5ad` |
| `Preprocessing.ipynb` | DE by condition on the integrated object, cluster renaming, re-export |
| `reclustering.ipynb` | Reload/re-export of the Palantir object; `forrds.h5ad` export for R |
| `pyscenic.ipynb` | Full pySCENIC pipeline (grn/ctx/aucell CLI calls) on the timecourse object |
| `dynamo.ipynb` | RNA velocity with dynamo (scanpy 1.10.4 / anndata 0.11.1 recorded in output) |
| `decon.ipynb` | Pseudobulk per cluster + Spearman/Pearson correlation against the bulk profile (a correlation-based concordance check, not a deconvolution method) |
| `Figures_plotting.ipynb` | General figure notebook (2024 iteration); panels overlap figs 4.4.5–4.4.8 and 4.4.17–4.4.18 |
| `Figure3K.ipynb` | C1QA/C1QC/LYVE1 expression vs Palantir pseudotime per cluster (provenance unresolved) |
| `Figure_velo_cellrank.ipynb` | PAGA + scVelo velocity (UMAP & diffusion map), velocity pseudotime, spliced/unspliced proportions (Fig 4.4.16), CellRank (VelocityKernel, RealTimeKernel + moscot, GPCCA) |

### Figure notebooks renamed from manuscript `Figure3*` names

| Original name | New name | Thesis figure(s) |
|---|---|---|
| `Figure3A_umaps.ipynb` | `fig4_4_5_to_8_umaps_dotplots.ipynb` | 4.4.5–4.4.8 (UMAP by Leiden, markers, per-condition density, LYVE1 dot plot) |
| `Figure3HI_heatmap_j.ipynb` | `fig4_4_9_to_11_marker_module_heatmaps.ipynb` | 4.4.9–4.4.11 (marker dot plot/heatmap, functional-module heatmap incl. endocytosis) |
| `Figure3N_pygsea.ipynb` | `fig4_4_14_mac4_go_enrichment.ipynb` | 4.4.14 (gseapy prerank, GO_Biological_Process_2021, Mac4 markers) |
| `Figure3G_Freq.ipynb` | `fig4_4_15_cluster_proportions.ipynb` | 4.4.15 (cluster proportions per condition over time) |
| `Figure3BCDEF_traj_velo.ipynb` | `fig4_4_17_to_19_velocity_trajectory.ipynb` | 4.4.17–4.4.19 (velocity streams UMAP/diffmap, pseudotime) |

Mapping established against the thesis figure captions; unresolved cases keep
their original names.

## `exploratory/`

- `integration.ipynb` — scVI/scANVI label-transfer sandbox (imports scrublet
  but never calls it — annotated in-file; harmony/bbknn trials). Not part of
  the final pipeline.
- `Figure3A_original.ipynb` — the 2026-02 first version of
  `fig4_4_5_to_8_umaps_dotplots.ipynb` (superseded 2026-07 version kept in
  `momac/`).
- `celloracle_pax8_exploratory.ipynb` — exploratory CellOracle run
  (condition-level GRN, PAX8 knockout only). CellOracle is **not** part of the
  thesis (see "Exploratory analyses" in the root README).

## Missing from this chapter

Figure 4.4.13 (MoMacVerse reverse projection) is implemented in
`../02-chapter3-bulk-rnaseq/03_momacverse_projection.R`; the forward
projection (4.4.12) is the `AddModuleScore(l1pos/l1neg)` block in
`00_demultiplex_htodemux_full.R`. Slingshot and Monocle3 (Fig 4.4.20) are
`10_`/`11_` above. CellOracle is **not reported in the thesis**; the only
executed perturbation, a PAX8 knockout, is kept in `exploratory/` (see
"Exploratory analyses" in the root README).
