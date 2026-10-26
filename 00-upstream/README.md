# 00-upstream — alignment and quantification

Upstream processing that produced the raw inputs for all downstream chapters.
All three scripts are **reconstructions**: they were rebuilt from the recorded
commands and pipeline outputs found on this machine, not copied from original
script files (the draft files in `~/Downloads` referenced by the author could
not be found; nothing was available to copy).

| Script | Role | Evidence |
|---|---|---|
| `00_preprocess_bulk.sh` | Bulk RNA-seq (human, paired-end Novogene): FastQC → Trimmomatic (PE, TruSeq3-PE, ILLUMINACLIP 2:30:7, LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15, MINLEN:36) → STAR paired-end (`--readFilesCommand zcat`, human toplevel FASTA + Ensembl 110 GTF; mouse FASTA unrecorded) → featureCounts (no `-p`) | candidate's recorded paired-end command; GTFs on disk; author's notes |
| `00_cellranger_count.sh` | 10x quantification per sample: `cellranger count` v7.0.0, GRCh38-3.0.0, per-sample GE + Antibody Capture libraries, 2-hashtag feature reference (reproduced, one id de-identified) | `_cmdline`/`_versions`/`_invocation`/feature-ref CSVs of the 8 pipestances under `bioinfo/singlecell_human/tmp/` |
| `02_velocyto.sh` | Spliced/unspliced looms: `velocyto run10x` per sample | loom paths in `momac/scvelo.ipynb` (`/home/yc/Desktop/loom/R1_*.loom`); looms not on this machine |

**No `cellranger multi` script is included**: every pipestance on disk is
`cellranger count` (`outs/filtered_feature_bc_matrix` directly under `outs/`,
no `per_sample_outs/`). Demultiplexing was done afterwards in Seurat
(`../03-chapter4-scrna-timecourse/00_demultiplex_htodemux_full.R`).

## Input formats on disk

Library CSV (one per sample; one Gene Expression and one Antibody Capture row), e.g. `R1_DMo_D0M_library.csv`:

```
fastqs,sample,library_type
/scdb/tmp/R1_DMo_D0M/,R1_DMo_D0M-SCI7T037-SCI5T037_HCJVTDSXC,Gene Expression
/scdb/tmp/R1_DMo_D0M_FB/,R1_DMo_D0M_FB-AK31823-AK31824_HCJVTDSXC,Antibody Capture
```

Feature reference (2 hashtags per sample; `R2`, pattern `5PNNNNNNNNNN(BC)`; the
cocktail-arm hashtag id/name is de-identified here — apply the identifier
mapping held outside this repository before depositing the originals):

```
id,name,read,pattern,sequence,feature_type
240316_Donor2_Mo_Day_1,240316_Donor2_Mo_Day_1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
240316_Donor2_MoMac_Day_0_M,240316_Donor2_MoMac_Day_0_M,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
240316_Donor2_MoMac_Day1_MDI,240316_Donor2_MoMac_Day1_MDI,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
240316_Donor2_MoMac_Day1_MDIX,240316_Donor2_MoMac_Day1_MDIX,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
240309_Donor2_PBMC_L1,240309_Donor2_PBMC_L1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
240229_Donor1_PBMC_L1,240229_Donor1_PBMC_L1,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
```

velocyto output: `run10x` writes `<sample>/velocyto/<sample>.loom` by default;
the looms actually read by the Chapter-4 notebook sit collected in
`/home/yc/Desktop/loom/` on the analysis server (moved there after generation).

Differences between these scripts and the submitted thesis Methods (6.1 / 7.1):
`cellranger count` vs multi, ILLUMINACLIP 2:30:7 vs 2:30:10, toplevel vs
primary assembly, featureCounts `-p` not recorded for the human bulk, and the
unverifiable velocyto repeat mask.
