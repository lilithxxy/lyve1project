# Cell Ranger 7.0.0 count per sample, reconstructed from recorded pipestance metadata.

# Count run (repeat per sample)
cellranger count \
  --id=R1_D1MDI_D1A \
  --libraries=./R1_D1MDI_D1A_library.csv \
  --transcriptome=<cellranger_ref_dir>/refdata-cellranger-GRCh38-3.0.0 \
  --feature-ref=./R1_D1MDI_D1A_feature_ref.csv \
  --localcores=32 \
  --localmem=128 \
  --chemistry=auto
