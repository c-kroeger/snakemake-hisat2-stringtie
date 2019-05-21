# Snakefile for HISAT2 Alignment and StringTie Assembly to performe genome-guided de novo assembly of transcripts.
# The pipeline additionally includes quality control of sequencing and alignment, class code assignement with gffcompare, and generation of count_tables.  

import pandas as pd

#Load configuration file
configfile: "config.yaml"


# Load sample_table from configuration file
sample_table = pd.read_csv(config["samples"], delimiter="\t").set_index("sample", drop=False)
samples = list(sample_table.index)

samples_PE = sample_table[sample_table["fq1"].notna() & sample_table["fq2"].notna()].index
samples_PE = list(samples_PE)

samples_SE = sample_table[sample_table["fq1"].notna() & sample_table["fq2"].isna()].index
samples_SE = list(samples_SE)

# input function for fastq files of PE and SE samples

def get_fastq_files(wildcards):
	return sample_table.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

#def get_fastq(wildcards):
 #   if (wildcards.sample in samples_PE) == True:
  #      return sample_table.loc[(wildcards.sample), ["fq1", "fq2"]]  
   # return sample_table.loc[(wildcards.sample), ["fq1"]]



rule all:
	input: 
		expand("ballgown/{sample}/{sample}.gtf", sample=samples),
		"GFFcompare.annotated.gtf",
		"qc/multiqc.html",
		"gene_count_matrix.csv",
		expand("mapped/{sample}.bam", sample=samples)


# Running the HISAT2 wrapper: Aligns reads and pipes output into samtools to convert into BAM file.
rule hisat2:
	input:
		reads = get_fastq_files
	output:
		"mapped/{sample}.bam"
	log:                                
		"logs/hisat2/{sample}.log"
	params:                             # idx is required, extra is optional
		idx=config["reference"]["index"],
		extra="--dta --new-summary" # --dta improves de novo assembly with Stringtie, --new-summary required for multiqc of hisat statistics
    	# threads: 1 # although 1 is chosen, still syntax error for this line, why?                        
	conda:
		"envs/hisat2.yaml"
	script:
		"scripts/hisat2_wrapper0.34.0.py"


# Running samtool sort wrapper: sorts BAM files
rule samtools_sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    params:
        "-m 4G"
    threads: 8
    wrapper:
        "0.34.0/bio/samtools/sort"


# Initial assembly of transcripts with StringTie for each sample 
rule stringtie_initial:
	input: 
		sbam="mapped/{sample}.sorted.bam",
		anno=config["reference"]["annotation"]
	output: 
		"sample_gtf/{sample}.gtf"
	log:
		"logs/stringtie/{sample}.log"
	threads: 8 
	shell: 
		# -l adds a label to assembled transcripts. Is it necessary? 
		"stringtie -v -p {threads} -G {input.anno} -o {output} {input.sbam} 2> {log}"


# Merge assembled transcripts from all samples. This produces a single, consistent set of transcripts and should overcome low coverage (leading to fragmented transcripts in some samples.
rule stringtie_merge:
	input: 
		gtf=expand("sample_gtf/{sample}.gtf", sample=samples), # or prepare a mergelist.txt of files.
		anno=config["reference"]["annotation"]
	output:
		"stringtie_merged.gtf"
	log:
		"logs/strintie_merge.log"
	threads:
		""
	shell:
		"stringtie -v --merge -p {threads} -G {input.anno} -o {output} {input.gtf} 2> {log}"


# Comparison of reference annotation with all transcripts assembled by stringtie --merge. Thereby usefull class codes are assigned describing the relation betwenn reference transcripts and assembled transcripts.
rule gffcompare_transcripts:
	input:
		st_transcripts="stringtie_merged.gtf",
		anno=config["reference"]["annotation"]
	output:
		"GFFcompare.annotated.gtf"
	shell:
		 "gffcompare -G -r {input.anno} -o GFFcompare {input.st_transcripts}"


# New stringtie assembly and quantification restricted to the transcript set produced by stringtie --merge. Additional input for ballgown analysis is prepared.
rule stringtie_ballgown:
	input: 
		anno="stringtie_merged.gtf", 	# important to use set of merged transcripts as annotation file!!
		sbam="mapped/{sample}.sorted.bam"
	output:
		"ballgown/{sample}/{sample}.gtf"
	threads:
		""
	shell: "stringtie -e -B -p {threads} -G {input.anno} -o {output} {input.sbam}"


# generation of gene and transcript count matrices from stringtie quanitification using the python 2.7 script coming with stringtie. 
rule count_tables:
	input:
		expand("ballgown/{sample}/{sample}.gtf", sample=samples)
	output:
		"gene_count_matrix.csv",
		"transcript_count_matrix.csv"
	conda:
		"envs/count_tables.yaml"
	shell:
		"python scripts/prepDE.py {input}"	#default read length is 75



##### quality control

# wrapper for fastqc 
rule fastqc:
    input:
        "chrX_data/samples/{sample_read}.fastq.gz"
    output:
        html="qc/fastqc/{sample_read}_fastqc.html",
        zip="qc/fastqc/{sample_read}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{sample_read}.log"
    wrapper:
        "0.34.0/bio/fastqc"



# multiqc summarising the fastqc files and hisat alignment statistics
rule multiqc:
	input:
		expand("qc/fastqc/{sample}_1_fastqc.zip", zip, sample=samples),
		expand("qc/fastqc/{sample}_2_fastqc.zip", zip, sample=samples_PE),
		expand("logs/hisat2/{sample}.log", sample=samples)
	output:
		"qc/multiqc.html"
	params:
		""  # Optional: extra parameters for multiqc
	log:
		"logs/multiqc.log"
	conda:
		"envs/multiqc.yaml"
	shell: 
		"multiqc {params} --force -n {output} {input} 2> {log}" 
	


