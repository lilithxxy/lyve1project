# Analysis code for "LYVE1⁺ macrophage" PhD thesis

This repository contains the analysis code behind the PhD thesis on LYVE1⁺
macrophage biology. It covers four analyses: bulk RNA-seq of sorted
macrophage populations, an in-vitro single-cell RNA-seq differentiation
timecourse, a human heart single-cell/single-nucleus macrophage atlas, and a
four-donor PBMC study of donor-to-donor variability.

## Repository layout (folders follow thesis chapters)

| Folder | Thesis chapter | Content |
|---|---|---|
| `00-upstream/` | Methods | Alignment & quantification scripts (reconstructed — see its README) |
| `01-chapter2-flow-qpcr/` | Chapter 2 | qRT-PCR figure (Figure 2.5.1) |
| `02-chapter3-bulk-rnaseq/` | Chapter 3 | Bulk RNA-seq DESeq2 + TF activity |
| `03-chapter4-scrna-timecourse/` | Chapter 4 | Time-resolved scRNA-seq of MDI vs MDIX differentiation |
| `04-chapter5-heart-atlas/` | Chapter 5 | Human heart sc/sn macrophage atlas |
| `05-chapter6-donor-variability/` | Chapter 6.2.3, Suppl. S7.x | Donor-to-donor variability in LYVE-1⁺ yield (four-donor PBMC pool) |
| `env/` | — | Environment specifications |

Each numbered folder has its own `README.md` listing the scripts in run order,
what each does, its inputs and outputs.

## Nomenclature (important)

Certain components of the differentiation cocktail and the antagonists used
against them are referred to by placeholder names throughout this repository
for intellectual-property reasons. Numeric concentrations are retained.

## Known reproducibility limitation

The analyses were run across several machines and conda environments over
2024–2026: notebooks record Python 3.10–3.12 kernels and both scvi-tools 1.2.0
and 1.3.0 (see `env/`). `env/environment.yml` records a single workable
environment with the only versions that are recorded in the sources pinned;
exact per-notebook environments were not captured at run time.

## Setup

Python side:

```bash
conda env create -f env/environment.yml
conda activate lyve1-thesis
```

R side: packages are listed (without invented versions) in
`env/sessionInfo.txt`; install with `BiocManager::install(...)` /
`install.packages(...)`. renv was never used for this project.

## Data availability

**Raw and intermediate data are NOT part of this repository.** All scripts
read inputs from absolute paths on the machines where they were run; these
paths must be adapted before the scripts are re-run. The Chapter 5 heart atlas
is built on public data from the DISCO database. A public data deposit of the
remaining data is planned:

> Data DOI (Zenodo): **[to be added]**

**Before any dataset is deposited alongside this code, the identifier mapping
held outside this repository must be applied to the data — the code in this
repository assumes the renamed labels and will not run against the original
data files.**

## Run order (summary)

1. `00-upstream/` — alignment, quantification, spliced/unspliced looms.
2. `01-chapter2-flow-qpcr/fig2_5_1_qpcr.R` — qPCR figure.
3. `02-chapter3-bulk-rnaseq/01_deseq2_analysis.R` then `02_decoupler_tf_activity.R` (same R session);
   `03_momacverse_projection.R`, `04_music_deconvolution.R`, `05_celldistinguisher_deconvolution.R`
   and `helpers/plotpca.R` are standalone follow-ups.
4. `03-chapter4-scrna-timecourse/` — hashtag demultiplexing (`00`, `00b`; R/Seurat), per-day
   QC/clustering (01–05), integration benchmark (06), cell-count figure (07), SCENIC regulon
   figures (08), temporal GSEA (09), Slingshot/Monocle3 pseudotime (10, 11),
   then the `momac/` analysis chain and figure notebooks.
5. `04-chapter5-heart-atlas/` — atlas integration chain (`atlas/`, `sn/`), annotation,
   scoring, and R plotting scripts.
6. `05-chapter6-donor-variability/` — Seurat merge/QC/annotation of the four-donor PBMC pool
   and the CellTypist QC notebook.

Details, inputs and outputs per script are in the per-folder READMEs.

## Thesis figure → script index

Only figures whose producing script is established in the chapter READMEs
are listed.

| Thesis figure | Script |
|---|---|
| Fig 2.4b (flow t-SNE) | Removed from the thesis |
| Fig 2.5.1 (qRT-PCR) | `01-chapter2-flow-qpcr/fig2_5_1_qpcr.R` |
| Figs 4.4.1–4.4.3 | No producing script identified |
| Fig 4.4.4 (cell numbers) | `03-chapter4-scrna-timecourse/07_cell_counts.py` |
| Figs 4.4.5–4.4.8 | `03-chapter4-scrna-timecourse/momac/fig4_4_5_to_8_umaps_dotplots.ipynb` |
| Figs 4.4.9–4.4.11 | `03-chapter4-scrna-timecourse/momac/fig4_4_9_to_11_marker_module_heatmaps.ipynb` |
| Fig 4.4.12 (forward projection) | `AddModuleScore(l1pos/l1neg)` block in `03-chapter4-scrna-timecourse/00_demultiplex_htodemux_full.R` |
| Fig 4.4.13 (MoMacVerse projection) | `02-chapter3-bulk-rnaseq/03_momacverse_projection.R` |
| Fig 4.4.14 (Mac4 GO) | `03-chapter4-scrna-timecourse/momac/fig4_4_14_mac4_go_enrichment.ipynb` |
| Fig 4.4.15 (cluster proportions) | `03-chapter4-scrna-timecourse/momac/fig4_4_15_cluster_proportions.ipynb` |
| Fig 4.4.16 (spliced/unspliced QC) | `03-chapter4-scrna-timecourse/momac/Figure_velo_cellrank.ipynb` |
| Figs 4.4.17–4.4.19 | `03-chapter4-scrna-timecourse/momac/fig4_4_17_to_19_velocity_trajectory.ipynb` |
| Fig 4.4.20 (pseudotime, 3 methods) | `10_slingshot_pseudotime.R`, `11_monocle3_pseudotime.R`, Palantir in `momac/uncorrected_new.ipynb` |
| Figs 4.5.1–4.5.3 (SCENIC MAF) | `03-chapter4-scrna-timecourse/08_scenic_regulons.ipynb` (twin: `04-chapter5-heart-atlas/scenic_heart.ipynb`) |
| Chapter 5 figures | Per-figure script mapping not established |
| Chapter 6.2.3 / Suppl. S7.x | `05-chapter6-donor-variability/` (folder-level) |

## Exploratory analyses (attempted, not reported in the thesis)

The following analyses were attempted during the project but did not produce
results robust enough to report, and they are not described in the thesis.
The code is included so that the work done can be inspected.

- In silico transcription-factor perturbation (CellOracle) —
  `03-chapter4-scrna-timecourse/exploratory/celloracle_pax8_exploratory.ipynb`
- Fate mapping (CellRank) —
  `03-chapter4-scrna-timecourse/momac/Figure_velo_cellrank.ipynb`
- Cell–cell interaction analysis (LIANA) —
  `04-chapter5-heart-atlas/exploratory/heart_sn_liana.ipynb`

## License

License to be confirmed pending IP review.
