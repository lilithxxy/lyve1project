# velocyto run10x per sample, reconstructed; the looms are read by momac/scvelo.ipynb.

# Spliced and unspliced quantification
velocyto run10x \
  -o /home/yc/Desktop/loom \
  <cellranger_sample_dir>/outs \
  /path/to/refdata-cellranger-GRCh38-3.0.0/genes.gtf
  # -m <rmsk.gtf> (not recorded whether a repeat mask was used)
