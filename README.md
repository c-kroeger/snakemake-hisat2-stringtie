Snakemake workflow: De novo transcriptome assembly with Hisat2 and StringTie. 

This pipeline comprises all steps of the 'new tuxedo suit' workflow published by Pertea et al. (1) and can be used to perform genome-guided de novo transcriptome assembly on bulk RNA-seq data with default parameters (without downstream analysis in R).
Additionally, the piepline comprises:
- adapter and quality trimmimg
- read quality control with FASTQC
- generation of a representative alingment file and bed file for the purpose of visualizing read coverage.



(1) Pertea, Mihaela; Kim, Daehwan; Pertea, Geo M.; Leek, Jeffrey T.; Salzberg, Steven L. (2016): Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. In: Nature Protocols 11, 1650 EP
