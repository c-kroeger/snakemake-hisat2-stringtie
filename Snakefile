# Author: Charlotte Kroeger
# Affiliation: Limes, Uni Bonn
# Aim: Snakefile for HISAT2 Alignment and StringTie Assembly to performe genome-guided de novo assembly of transcripts. 
#      The pipeline additionally includes quality control of sequencing and alignment, class code assignement with gffcompare, and generation of count_tables
# Date: 11 July 2019      
# Run: snakemake --use-conda --snakefile Snakefile


import pandas as pd

#### Load configfile and sample table ####

#Load configuration file
configfile: "config.yaml"

# Load sample_table from configuration file
sample_table = pd.read_csv(config["samples"], delimiter="\t").set_index("sample", drop=False)
samples = list(sample_table.index)



#### Defining input functions ####

def get_fastq_r1(wildcards):
	return sample_table.loc[(wildcards.sample), "fq1"]

def get_fastq_r2(wildcards):
	return sample_table.loc[(wildcards.sample), "fq2"]



#### target rule ####

rule all:
	input:
		expand("output/mapped_sorted/{sample}.sorted.bam.bai", sample=samples), 
		expand("output/ballgown/{sample}/{sample}.gtf", sample=samples),
		"output/gffcompare/GFFcompare.annotated.gtf",
		#"output/qc/multiqc_raw_reads.html",
		"output/qc/multiqc_trimmed_reads.html",
		"output/gene_count_matrix.csv",
		"output/mapped_sorted/merged_subsample.bw"



#### snakemake workflow ####

#### Trimming and Alignment ####

# Running Trimmomatic wrapper for PE reads:
rule trimmomatic_pe:
	input:
		r1 = get_fastq_r1,
		r2 = get_fastq_r2
	output:
		r1="output/trimmed/{sample}_1.fastq.gz",
		r2="output/trimmed/{sample}_2.fastq.gz",
		#reads whithout mate after trimming
		r1_unpaired="output/trimmed/{sample}_1.unpaired.fastq.gz",
		r2_unpaired="output/trimmed/{sample}_2.unpaired.fastq.gz"
	log:
		"output/logs/trimmomatic/{sample}.log"
	params:
		#list of trimmers:
		trimmer=["ILLUMINACLIP:{}:2:30:10".format(config["trimmomatic"]["adapters"]),
                         "LEADING:3",
                         "TRAILING:3",
                         "SLIDINGWINDOW:4:25",
                         "MINLEN:36"],
		#optional parameters
		extra="",
		compression_level="-9"
	threads: 8
	wrapper:
		"0.35.0/bio/trimmomatic/pe"


# Running the HISAT2 wrapper: Aligns reads and pipes output into samtools to convert into BAM file.
# Note: wrapper was split into environment and script to solve different python requirements of main workflow and wrapper.
rule hisat2:
	input:
		reads = ["output/trimmed/{sample}_1.fastq.gz", "output/trimmed/{sample}_2.fastq.gz"]
	output:
		temp("output/mapped/{sample}.bam")
	log:                                
		"output/logs/hisat2/{sample}.log"
	params:
    	# idx is required, extra is optional
		idx=config["reference"]["index"],
		extra="--dta --new-summary" # --dta improves de novo assembly with Stringtie, --new-summary required for multiqc of hisat statistics
	threads: 8
	conda:
		"envs/hisat2.yaml"
	script:
		"scripts/hisat2_wrapper0.34.0.py"


# Running samtool sort wrapper: sorts BAM files
rule samtools_sort:
	input:
		"output/mapped/{sample}.bam"
	output:
		"output/mapped_sorted/{sample}.sorted.bam"
	params:
		"-m 4G"
	threads: 8
	wrapper:
		"0.35.1/bio/samtools/sort"
		
rule samtools_index:
	input: "output/mapped_sorted/{sample}.sorted.bam"
	output: "output/mapped_sorted/{sample}.sorted.bam.bai"
	params:
		"" # optional params string
	wrapper:
		"0.35.1/bio/samtools/index"		
		
		
#### Assembly and Quantificaion ####

# Initial assembly of transcripts with StringTie for each sample 
rule stringtie_initial:
	input: 
		sbam="output/mapped_sorted/{sample}.sorted.bam",
		anno=config["reference"]["annotation"]
	output:
		"output/sample_gtf/{sample}.gtf"
	log:
		"output/logs/stringtie/{sample}.log"
	threads: 8 
	conda:
		"envs/stringtie.yaml"
	shell: 
		# -l adds a label to assembled transcripts. Is it necessary? 
		"stringtie -v -p {threads} -G {input.anno} -o {output} {input.sbam} 2> {log}"


# Merge assembled transcripts from all samples. This produces a single, consistent set of transcripts and should overcome low coverage (leading to fragmented transcripts in some samples.
rule stringtie_merge:
	input: 
		gtf=expand("output/sample_gtf/{sample}.gtf", sample=samples), # or prepare a mergelist.txt of files.
		anno=config["reference"]["annotation"]
	output:
		"output/stringtie_merged_assembly.gtf"
	log:
		"output/logs/stringtie_merge.log"
	threads: 8
	conda:
		"envs/stringtie.yaml"
	shell:
		"stringtie -v --merge -p {threads} -G {input.anno} -o {output} {input.gtf} 2> {log}"


# Comparison of reference annotation with all transcripts assembled by stringtie --merge. Thereby usefull class codes are assigned describing the relation betwenn reference transcripts and assembled transcripts.
rule gffcompare_transcripts:
	input:
		st_transcripts="output/stringtie_merged_assembly.gtf",
		anno=config["reference"]["annotation"]
	output:
		"output/gffcompare/GFFcompare.annotated.gtf"
	threads: 2
	conda:
		"envs/gffcompare.yaml"
	shell:
		"gffcompare -G -r {input.anno} -o output/gffcompare/GFFcompare {input.st_transcripts}"


# New stringtie assembly and quantification restricted to the transcript set produced by stringtie --merge. Additional input for ballgown analysis is prepared.
rule stringtie_ballgown:
	input: 
		anno="output/stringtie_merged_assembly.gtf", 	# important to use set of merged transcripts as annotation file!!
		sbam="output/mapped_sorted/{sample}.sorted.bam"
	output:
		"output/ballgown/{sample}/{sample}.gtf"
	threads: 8
	conda:
		"envs/stringtie.yaml"
	shell: "stringtie -e -B -p {threads} -G {input.anno} -o {output} {input.sbam}"


# generation of gene and transcript count matrices from stringtie quanitification using the python 2.7 script coming with stringtie. 
rule count_tables:
	input:
		expand("output/ballgown/{sample}/{sample}.gtf", sample=samples)
	output:
		"output/gene_count_matrix.csv",
		"output/transcript_count_matrix.csv"
	conda:
		"envs/count_tables.yaml"
	shell:
		"""
		cd output
		prepDE.py -i ballgown
		cd ..
		"""




#### Quality control of raw reads ####

# fastqc of read 1
rule fastqc_r1:
	input:
		get_fastq_r1
	output:
		html="output/qc/fastqc_raw/{sample}_1_fastqc.html",
		zip="output/qc/fastqc_raw/{sample}_1_fastqc.zip"
	params: ""
	log:
		"output/logs/fastqc_raw/{sample}_1.log"
	wrapper:
		"0.34.0/bio/fastqc"


# fastqc of read 2 
rule fastqc_r2:
	input:
		get_fastq_r2
	output:
		html="output/qc/fastqc_raw/{sample}_2_fastqc.html",
		zip="output/qc/fastqc_raw/{sample}_2_fastqc.zip"
	params: ""
	log:
		"output/logs/fastqc_raw/{sample}_2.log"
	wrapper:
		"0.34.0/bio/fastqc"


# multiqc summarising the fastqc files of raw reads
rule multiqc:
	input:
		expand("output/qc/fastqc_raw/{sample}_1_fastqc.zip", zip, sample=samples),
		expand("output/qc/fastqc_raw/{sample}_2_fastqc.zip", zip, sample=samples)
	output:
		"output/qc/multiqc_raw_reads.html"
	log:
		"output/logs/multiqc_raw_reads.log"
	wrapper:
		"0.35.1/bio/multiqc"




#### Quality control of trimmed reads ####

# fastqc of trimmed reads
rule fastqc_trimmed:
	input:
		"output/trimmed/{sample_read}.fastq.gz"
	output:
		html="output/qc/fastqc_trimmed/{sample_read}_fastqc.html",
		zip="output/qc/fastqc_trimmed/{sample_read}_fastqc.zip"
	params: ""
	log:
		"output/logs/fastqc_trimmed/{sample_read}.log"
	wrapper:
		"0.34.0/bio/fastqc"


# multiqc summarising the fastqc files after trimming and hisat alignment statistics
rule multiqc_trimmed:
	input:
		expand("output/qc/fastqc_trimmed/{sample}_1_fastqc.zip", zip, sample=samples),
		expand("output/qc/fastqc_trimmed/{sample}_2_fastqc.zip", zip, sample=samples),
		expand("output/logs/hisat2/{sample}.log", sample=samples),
		expand("output/logs/trimmomatic/{sample}.log", sample=samples)
	output:
		"output/qc/multiqc_trimmed_reads.html"
	params:
		"--config envs/multiqc_config.yaml"  # Optional: extra parameters for multiqc
	log:
		"output/logs/multiqc_trimmed_reads.log"
	wrapper:
		"0.35.1/bio/multiqc"




#### subsample and bigwig ####

# subsample bam files
rule subsample:
    input:
        "output/mapped_sorted/{sample}.sorted.bam"
    output:
        temp("output/mapped_sorted/{sample}.subsample.bam")
    log:
        "output/logs/subsample/{sample}.log"
    conda:
        "envs/sambamba.yaml"
    threads:
        8
    shell:
        """
        nreads=$(samtools view -c {input})
        rate=$(echo "scale=5;100000/$nreads" | bc)
        sambamba view -f bam -t 5 --subsampling-seed=42 -s $rate {input} | samtools sort -m 4G -@ 8 -T - > {output} 2> {log}
        """

# merge subsamples into one bam file
rule merge_subsamples:
    input:
        expand("output/mapped_sorted/{sample}.subsample.bam", sample=samples)
    output:
        "output/mapped_sorted/merged_subsample.bam"
    params:
        "" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "0.35.1/bio/samtools/merge"


# index merged bam file
rule index_merged_subsample:
    input:
        "output/mapped_sorted/merged_subsample.bam"
    output:
        "output/mapped_sorted/merged_subsample.bai"
    params:
        "" # optional params string
    wrapper:
        "0.35.1/bio/samtools/index"


# generate a bigwig file of merged bam
rule bigwig_subsample:
    input:
        bam="output/mapped_sorted/merged_subsample.bam",
        bai="output/mapped_sorted/merged_subsample.bai"
    output:
        "output/mapped_sorted/merged_subsample.bw"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output}"

