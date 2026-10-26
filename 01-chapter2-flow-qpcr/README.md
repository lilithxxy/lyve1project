# 01-chapter2-flow-qpcr — thesis Chapter 2

## Scripts

### `fig2_5_1_qpcr.R` — thesis Figure 2.5.1 (qRT-PCR)

- **Renamed at source.** The original file was `Chapter2_Figures/Figure2_4b.R`,
  but its contents are the qRT-PCR figure, not Figure 2.4b (the flow-cytometry
  t-SNE, for which no code exists).
- **What it does:** log2 2^-ΔΔCt bar plots with error bars and donor jitter
  points for three genes across three groups, with paired *t*-tests
  (`rstatix::pairwise_t_test`, `ggpubr::stat_pvalue_manual`), combined via
  `patchwork`.
- **Genes:** LYVE1, F13A1, MRC1.
- **Groups:** `Day0` (monocytes), `Day5_MDI`, `Day5_MDIX`.
- **Inputs:** none from disk — the ΔΔCt values are hard-coded matrices in the
  script (LYVE1: 5 donors × 3 groups; F13A1 and MRC1: 3 donors × 3 groups).
  Rows = donors, columns = Day0 / Day5_MDI / Day5_MDIX.
- **Outputs:** `qpcr_LYVE1.png`, `qpcr_F13A1.png`, `qpcr_MRC1.png` (300 dpi)
  and `qpcr_combined.png`.
- **Libraries:** ggplot2, rstatix, ggpubr, dplyr, patchwork.

## Missing from this chapter

Figure 2.4b (flow-cytometry t-SNE) and the flow-cytometry quantification of
Figure 4.2.1 have no code in this repository (FCS files and analysis were
external).
