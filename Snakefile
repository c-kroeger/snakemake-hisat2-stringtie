# Snakefile for HISAT2 Alignment and StringTie Assembly to performe genome-guided de novo assembly of transcripts.
# The pipeline additionally includes quality control of sequencing and alignment, class code assignement with gffcompare, and generation of count_tables.  

configfile: "config.yaml"

rule all:
	input: 
		expand("ballgown/{sample}/{sample}.gtf", sample=config["samples"]),
		"GFFcompare.annotated.gtf",
		expand("qc/fastqc/{sample}_fastqc.html", sample=config["samples"]),
		"qc/multiqc.html",
		"qc/hisat.html",
		"gene_count_matrix.csv"


# Running the HISAT2 wrapper: Aligns reads and pipes output into samtools to convert into BAM file.
rule hisat2:
	input:
		reads = lambda wildcards: config["samples"][wildcards.sample]
	output:
		"mapped/{sample}.bam"
	log:                                
		"logs/hisat2/{sample}.log"
	params:                             # idx is required, extra is optional
		idx="chrX_data/indexes/chrX_tran",
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
		anno="chrX_data/genes/chrX.gtf"
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
		gtf=expand("sample_gtf/{sample}.gtf", sample=config["samples"]), # or prepare a mergelist.txt of files.
		anno="chrX_data/genes/chrX.gtf"
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
		anno="chrX_data/genes/chrX.gtf"
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
		expand("ballgown/{sample}/{sample}.gtf", sample=config["samples"])
	output:
		"gene_count_matrix.csv",
		"transcript_count_matrix.csv"
	conda:
		"envs/test.yaml"
	shell:
	#	"python .snakemake/conda/0cc9cc6b/bin/prepDE.py {input}"	
		"python scripts/prepDE.py {input}"	#default read length is 75
	# is it possible to access path of rule environment by an variable? 



##### quality control

# modified merge rule from Jonas
rule fastq_merge:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
       temp("chrX_data/fastq/{sample}.fastq.gz")
       #"chrX_data/fastq/{sample}.fastq.gz"
    run:
        if len(input) == 1:
            shell("ln -s {input} {output}")
        else:
            shell("cat {input} > {output}")


# wrapper for fastqc 
rule fastqc:
    input:
        "chrX_data/fastq/{sample}.fastq.gz"
    output:
        html="qc/fastqc/{sample}_fastqc.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.34.0/bio/fastqc"


# wrapper for multiqc summarising the fastqc files.
rule multiqc_of_fastqc:
	input:
		expand("qc/fastqc/{sample}_fastqc.zip", zip, sample=config["samples"])
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
	

# wrapper for multiqc summarising the hisat2 alignment.
rule multiqc_of_hisat:
	input:
        	expand("logs/hisat2/{sample}.log", sample=config["samples"])
	output:
		"qc/hisat.html"
	params:
		""  # Optional: extra parameters for multiqc
	log:
		"logs/multiqc_hisat.log"
	conda:
		"envs/multiqc.yaml"
	shell: 
		"multiqc {params} --force -n {output} {input} 2> {log}" 


