# ATAC_pipeline
The Transposase-Accessible Chromatin assay utilises the capability of Tn5 transposase and provides a way to identify the regions that are bound by transcription factors across the genome. The peaks which are the accessible regions can be compared among different samples of interest. To this end, several pipelines have been developed to identify the changes in molecular networks under different conditions. Here, we propose a pipeline for the analysis of ATAC-seq data. The key steps included in this pipeline include:
1. Data pre-processing
2. Alignment and filtering
3. Peak identification 
4. Downstream analysis

The downstream analysis includes the motif enrichment, footprinting, annotation and integrative analysis with RNA-seq and ChIP-Atlas data. We provide a set of scripts which the user can apply to any ATAC-seq data to perform similar analysis. The pipeline can be submitted as a job on a cluster to perform all the analysis as mentioned above.

The following packages are required to be installed by the user for running the pipeline:
  - bedtools
  - bowtie2
  - deeptools
  - fastqc
  - homer
  - macs2
  - pip
  - python>=3.9
  - samtools
  - trim-galore
  - tobias
  - R
    - ggplot2
    - cowplot
    - reshape2
    - dplyr
    - tidyr
    - RColorBrewer
    - DiffBind

All the required files and scripts can be installed by running the command:
    git clone ......

After cloning the repository, the user can set the working directory in the CONFIG_ATAC file. The user can also modify the sampleID.txt file with FASTQ file names of interest to user and the corresponding name for the treatment conditions. After making the required changes, the user can run the pipeline with the command
    sbatch atac.sh

