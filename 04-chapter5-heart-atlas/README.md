# 04-chapter5-heart-atlas — thesis Chapter 5

Human heart single-cell (sc) + single-nucleus (sn) macrophage atlas built from
public DISCO/HCA downloads, with scVI integration, myeloid→macrophage lineage
subsetting and annotation (clusters Mac0–Mac6).

## Top-level scripts (final analysis)

| Script | What it does | Inputs | Outputs |
|---|---|---|---|
| `annotation_sc.ipynb` | sc atlas annotation: Leiden→Mac0–Mac6 map, LYVE1-core scoring, marker/functional-module plots, metadata joins, disease/anatomy/age composition (Mann–Whitney + BH), pseudobulk export, decoupleR TF activity | `qc3_scvi.h5ad`, `heart_sc_clean_.csv` | `anno_layer1.csv`, pseudobulk TSVs, figures |
| `annotation_sn.ipynb` | sn atlas annotation (same workflow; drops one cluster, sn-specific module dict) | `sn_final.h5ad`, sample CSVs | `anno_layer1.csv` (sn), pseudobulk TSVs |
| `chapter5_sc.ipynb` | sc atlas figure factory: sample pies, UMAPs, dotplots, composition plots, adult-vs-feotus stats, HMOX1/heme panel | `qc3_scvi.h5ad`, `heart_sc_clean_.csv` | figures in `Chapter5_Figures/` paths |
| `scoring_and_annotation.ipynb` | UCell signature scoring (17 sc / 15 sn modules) + decoupleR TF ranking for both atlases | `qc3_scvi.h5ad`, `sn_final.h5ad` | `*_UCell` obs columns, TF matrixplots |
| `GO_heatmap.R` | curated GO dotplot + heatmap across Mac clusters | **`sc_go_curated.csv`** (configuration, see below) | `heart_go_dotplot.png` |
| `GSEA_loop.R` | clusterProfiler GO/GSEA/KEGG loop over per-cluster DE (twin of `03-chapter4.../09_temporal_gsea.R`) | `sc_deg.csv` | `sc_deg/`, `sc_go/`, `sc_gsea/` CSVs/RDS, gseaplot2 figures |
| `barplot_geneontology.R` | detox-term barplot (Mac4) + fibrosis barplot + wound-healing GSEA plot (Mac2) | `sc_go/heart_sc4_go.csv`, `sc_go/heart_sc2_go.csv`, `sc_gsea/heart_sc2_gsea.csv` | PNG figures |
| `sc_mac2.R` | standalone cleaned Mac2 fibrosis/wound-healing figures | same as above | PNG figures |
| `scenic_heart.ipynb` | markdown-narrated twin of `03-chapter4.../08_scenic_regulons.ipynb` (SCENIC MAF figures). **Provenance note:** despite living under `Chapter5_Figures/` at source, it analyses the Chapter-4 in-vitro timecourse (thesis Figs 4.5.1–4.5.3), not the heart atlas; kept here to preserve the source location. | `inte_scenic.loom`, `inte_palantir.h5ad` | 3 SCENIC MAF figures (PDF) |

**`sc_go_curated.csv` is hand-written configuration, not generated data** (72
rows; columns `Clusters, ID, Description, p.adjust, RichFactor`; the 12 GO
descriptions match `custom_order` in `GO_heatmap.R`). It is the only CSV in
this repository.

## `atlas/` — heart single-cell chain (run order)

1. `readh5.ipynb` — 10x `.h5` → per-sample `.h5ad` (`/home/yanjun3/atlas/heart_sc`).
2. `integration_preannotation_heartsc.ipynb` — per-sample QC loop → concat →
   scVI → Leiden → myeloid subsetting. Also contains a second-pass filter.
3. `integration_preannotation_group2/3/4.ipynb` — heart-**sn** sample groups,
   same QC loop + scVI + myeloid subsetting (these analyse sn groups; named
   "groupN" at source).
4. `integration_preannotation_myeloid.ipynb` — myeloid round: re-scVI, marker
   checks, macrophage subsetting → `qc3.h5ad`.
5. `integration_preannotation_macrophages.ipynb` — macrophage round: scVI →
   `qc3_scvi.h5ad`, Leiden 0.3/0.4/0.5, DE, final subsets.
6. `integration_sc_and_sn.ipynb` (+ `_cleanhvg`, `_rmrpl`, `_reintegrate`,
   `_labeltransfer` variants) — sc+sn combined scVI (continuous covariates
   mt/genes/counts), embedding-density by `rna_source`, DE, macrophage-marker
   dotplot → `snsc_final.h5ad`. `_cleanhvg` and `_rmrpl` vary the feature set
   fed to scVI (MT-/RPL/RPS genes dropped); `_labeltransfer` reads a different
   input; `_reintegrate` is a re-run writing to `figure_scsn_reintegrate`
   (see the rename note below).
   **Rename note:** `integration_sc_and_sn_reintegrate.ipynb` was called
   `integration_sc_and_sn_cellbender.ipynb` at source, but it contains no
   CellBender code — it is the same sc+sn scVI integration re-run, differing
   only in extra mt/ribo QC covariates, a different figure directory, and
   dropping the ribo/mts score plots. Its `adata_clean` (MT-/RPL-/RPS-removed)
   object is created but never used — scVI runs on the full object.
7. `celltypist.ipynb` — trains custom CellTypist models from `new_anno`,
   cross-annotates sc ↔ sn.
8. `celltypist_benchmark.ipynb` — QC + scrublet + CellTypist benchmark loop
   on an **adipose** reference dataset; scVI/scANVI/BBKNN scIB benchmarking.
9. `Decoupler.ipynb` — decoupleR TF activity on the day-4 in-vitro data
   (MDIX vs mean).
10. `annotation_sc.ipynb`, `annotation_sn.ipynb`, `scoring_and_annotation.ipynb`
    (top-level copies; the `atlas/` copies are the server-side originals of the
    same stages — `scoring_and_annotation` is byte-identical, the annotation
    pairs differ by run date).

## `sn/` — heart single-nucleus chain (run order)

1. `readh5.ipynb` — sn 10x → h5ad (differs from the atlas/ copy: different
   input directory and sample handling).
2. `integration_preannotation_myeloid.ipynb` — merges the per-group myeloid
   objects (`group*_myeloid.h5ad`) → `heart_myeloid_merged.h5ad` → scVI.
   (Differs substantially from the atlas/ sc-myeloid notebook.)
3. `integration_preannotation_v2.ipynb` → `v3.ipynb` → `v4.ipynb` — iterative
   macrophage refinement; v3 applies a mito filter (see cutoffs below) and
   cell-cycle scoring; v4 removes two clusters and writes the final
   `sn_final.h5ad`.

**Duplicates:** the `sn/` copies of `celltypist_benchmark.ipynb`,
`integration_preannotation_group2/3/4.ipynb`, `celltypist.ipynb` and
`annotation_sn.ipynb` were byte-identical to their `atlas/` counterparts
(sha256-verified) and are stored once, under `atlas/`. The two
`after_integration*` notebooks were exploratory work on an intestine reference
dataset — see `exploratory/`.

## QC cutoffs actually applied in this chapter

| Stage | Filter | File |
|---|---|---|
| Heart per-sample loop (sc and sn groups) | `sc.pp.filter_cells(adata, min_genes=300)`; `sc.pp.filter_genes(adata, min_cells=10)` | `atlas/integration_preannotation_heartsc.ipynb` (cell 4), group2/3/4 |
| Heart sc, second pass | `n_genes_by_counts` 500–10,000; `total_counts` 800–5,000; `pct_counts_mt < 5` | `atlas/integration_preannotation_heartsc.ipynb` (cell 17) |
| Heart sn pass | `pct_counts_mt < 2` | `sn/integration_preannotation_v3.ipynb` (cell 6) |
| Adipose benchmark loop | min_genes 300; `pct_counts_mt < 20`; min_cells 10; `sc.pp.scrublet` | `atlas/celltypist_benchmark.ipynb` (cell 4) |

No ambient-RNA removal (SoupX/CellBender/DecontX) and no doublet detection
were run on the heart data; `sc.pp.scrublet` runs only in the adipose
benchmark loop and its score columns are never filtered on.

## Input formats (not in repo)

- `qc3_scvi.h5ad` (6.5 GB) — final sc macrophage object; obs incl. `sample`,
  QC metrics, `leiden_res0_3/0_4/0_5`, scVI slots.
- `heart_sc_clean.csv` / `heart_sn.csv` — per-sample metadata joined on
  `sample_id`; columns incl. `project_id`, `tissue`, `anatomical_site`,
  `Disease`, `age_group`, `age`, `gender`, `race`, `rna_source`, `study`.
- `sc_deg.csv` — per-cluster DE chunks (5 columns per cluster, see Chapter-4 README).

## `exploratory/`

- `after_integration.ipynb` — post-integration analysis of an intestine
  reference dataset (scIB metrics, marker dicts); source of the
  `bulktop50`/`mac4top50`/`mac7top100` signature gene lists used for scoring
  elsewhere.
- `heart_sn_liana.ipynb` — LIANA cell–cell interaction analysis of the heart
  sn macrophage subclusters (reload from 4 h5ads → concat → LIANA →
  visualisations, 10 groups in `anno_layer1`). Not reported in the thesis.
