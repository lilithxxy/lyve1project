# Bulk RNA-seq preprocessing (FastQC, Trimmomatic, STAR, featureCounts), reconstructed from recorded commands.

# FastQC
fastqc <sample>_1.clean.fq.gz <sample>_2.clean.fq.gz -threads 12 -o <qc_dir>

# Trimmomatic (paired-end)
java -jar <trimmomatic_dir>/trimmomatic.jar PE \
  -threads 12 \
  -phred33 \
  <sample>_1.clean.fq.gz <sample>_2.clean.fq.gz \
  trimmed_<sample>_1.fq.gz trimmed_<sample>_unpaired_1.fq.gz \
  trimmed_<sample>_2.fq.gz trimmed_<sample>_unpaired_2.fq.gz \
  ILLUMINACLIP:<trimmomatic_dir>/adapters/TruSeq3-PE.fa:2:30:7 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# STAR genome index (human)
STAR --runMode genomeGenerate \
  --genomeDir Hs_STAR_dir/ \
  --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.toplevel.fa \
  --sjdbGTFfile Homo_sapiens.GRCh38.110.gtf \
  --runThreadN 12

# STAR genome index (mouse; genome FASTA unrecorded)
STAR --runMode genomeGenerate \
  --genomeDir Mm_STAR_dir/ \
  --genomeFastaFiles <mouse_genome.fa> \
  --sjdbGTFfile Mus_musculus.GRCm39.110.gtf \
  --runThreadN 12

# STAR alignment (paired-end, gzipped input)
STAR --runMode alignReads \
  --genomeDir Hs_STAR_dir/ \
  --outSAMtype BAM SortedByCoordinate \
  --readFilesIn trimmed_<sample>_1.fq.gz trimmed_<sample>_2.fq.gz \
  --readFilesCommand zcat \
  --runThreadN 12 \
  --outFileNamePrefix trimmed_<sample>.bam

# featureCounts
featureCounts \
  -a Homo_sapiens.GRCh38.110.gtf \
  -T 12 \
  -o <sample>_counts.txt \
  trimmed_<sample>.bam

# Clean up counts table
cat <sample>_counts.txt | cut -f1,7 > cut_<sample>_counts.txt
