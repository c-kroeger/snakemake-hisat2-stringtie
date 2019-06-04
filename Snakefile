# Snakefile for HISAT2 Alignment and StringTie Assembly to performe genome-guided de novo assembly of transcripts.
# The pipeline additionally includes quality control of sequencing and alignment, class code assignement with gffcompare, and generation of count_tables.  


import pandas as pd

#### Load configfile and sample table ####

#Load configuration file
configfile: "config.yaml"

# Load sample_table from configuration file
sample_table = pd.read_csv(config["samples"], delimiter="\t").set_index("sample", drop=False)
samples = list(sample_table.index)



#### Defining input function ####

#def get_fastq(wildcards):
#    return (sample_table.loc[(wildcards.sample), ["fq1", "fq2"]])


def get_fastq_r1(wildcards):
	return sample_table.loc[(wildcards.sample), "fq1"]

def get_fastq_r2(wildcards):
	return sample_table.loc[(wildcards.sample), "fq2"]



#### target rule ####

rule all:
	input: 
		expand("ballgown/{sample}/{sample}.gtf", sample=samples),
		"GFFcompare.annotated.gtf",
		"qc/multiqc.html",
		"gene_count_matrix.csv",
		#expand("mapped/{sample}.bam", sample=samples),
		#expand("mapped_sorted/{sample}.sorted.bam", sample=samples),
		expand("trimmed/{sample}_1.fastq.gz", sample=samples),
		expand("trimmed/{sample}_2.fastq.gz", sample=samples)


#### snakemake rules ####

# Running Trimmomatic wrapper for PE reads:
rule trimmomatic_pe:
	input:
		r1 = get_fastq_r1,
		r2 = get_fastq_r2
	output:
		r1="trimmed/{sample}_1.fastq.gz",
		r2="trimmed/{sample}_2.fastq.gz",
		#reads whithout mate after trimming
		r1_unpaired="trimmed/{sample}_1.unpaired.fastq.gz",
		r2_unpaired="trimmed/{sample}_2.unpaired.fastq.gz"
	log:
		"logs/trimmomatic/{sample}.log"
	params:
		#list of trimmers:
		trimmer=["ILLUMINACLIP:trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10","LEADING:3","TRAILING:3","SLIDINGWINDOW:4:15", "MINLEN:20"],
		#optional parameters
		extra="",
		compression_level="-9"
	threads:
		2
	wrapper:
		"0.35.0/bio/trimmomatic/pe"


# Running the HISAT2 wrapper: Aligns reads and pipes output into samtools to convert into BAM file.
rule hisat2:
	input:
		reads = ["trimmed/{sample}_1.fastq.gz", "trimmed/{sample}_2.fastq.gz"]
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
		"mapped_sorted/{sample}.sorted.bam"
	params:
		"-m 4G"
	threads: 8
	wrapper:
		"0.34.0/bio/samtools/sort"


# Initial assembly of transcripts with StringTie for each sample 
rule stringtie_initial:
	input: 
		sbam="mapped_sorted/{sample}.sorted.bam",
		anno=config["reference"]["annotation"]
	output: 
		"sample_gtf/{sample}.gtf"
	log:
		"logs/stringtie/{sample}.log"
	threads: 8 
	conda:
		"envs/stringtie.yaml"
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
		"logs/stringtie_merge.log"
	threads:
		""
	conda:
		"envs/stringtie.yaml"
	shell:
		"stringtie -v --merge -p {threads} -G {input.anno} -o {output} {input.gtf} 2> {log}"


# Comparison of reference annotation with all transcripts assembled by stringtie --merge. Thereby usefull class codes are assigned describing the relation betwenn reference transcripts and assembled transcripts.
rule gffcompare_transcripts:
	input:
		st_transcripts="stringtie_merged.gtf",
		anno=config["reference"]["annotation"]
	output:
		"GFFcompare.annotated.gtf"
	conda:
		"envs/gffcompare.yaml"
	shell:
		 "gffcompare -G -r {input.anno} -o GFFcompare {input.st_transcripts}"


# New stringtie assembly and quantification restricted to the transcript set produced by stringtie --merge. Additional input for ballgown analysis is prepared.
rule stringtie_ballgown:
	input: 
		anno="stringtie_merged.gtf", 	# important to use set of merged transcripts as annotation file!!
		sbam="mapped_sorted/{sample}.sorted.bam"
	output:
		"ballgown/{sample}/{sample}.gtf"
	threads:
		""
	conda:
		"envs/stringtie.yaml"
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



#### quality control ####

# wrapper for fastqc 
rule fastqc:
    input:
        "trimmed/{sample_read}.fastq.gz"
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
		expand("qc/fastqc/{sample}_2_fastqc.zip", zip, sample=samples),
		expand("logs/hisat2/{sample}.log", sample=samples),
		expand("logs/trimmomatic/{sample}.log", sample=samples)
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
	


