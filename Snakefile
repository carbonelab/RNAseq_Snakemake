# Snakefile to analyze RRBS PE data
# 

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/fastqc/raw/multiqc_report.html"),
        expand("data/trimming/out_paired_{dir}_{sample}.fastq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/unpaired/out_unpaired_{dir}_{sample}.fastq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/fastqc/trim/out_paired_{dir}_{sample}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/fastqc/trim/multiqc_report.html"),
        expand("data/star/{sample}.star.ReadsPerGene.out.tab", sample = SAMPLES),
        "data/counts_table.txt",
        "data/ide/filtered_counts.txt"

rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip",
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule multiqc_raw:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"])
    output:
        "data/fastqc/raw/multiqc_report.html"
    params:
        workingDir = "data/fastqc/raw/*"
    shell:
        "multiqc -o data/fastqc/raw/ {params.workingDir}"

rule trimmomatic:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        "data/trimming/out_paired_R1_{sample}.fastq.gz",
        "data/trimming/unpaired/out_unpaired_R1_{sample}.fastq.gz",
        "data/trimming/out_paired_R2_{sample}.fastq.gz",
        "data/trimming/unpaired/out_unpaired_R2_{sample}.fastq.gz"
    params:
        paired_R1 = "data/trimming/out_paired_R1_{sample}.fastq.gz",
        unpaired_R1 = "data/trimming/unpaired/out_unpaired_R1_{sample}.fastq.gz",
        paired_R2 = "data/trimming/out_paired_R2_{sample}.fastq.gz",
        unpaired_R2 = "data/trimming/unpaired/out_unpaired_R2_{sample}.fastq.gz",
        leading = "LEADING:3",
        trailing = "TRAILING:3",
        sldWindow = "SLIDINGWINDOW:4:15",
        minLen = "MINLEN:36"
    shell:
        "trimmomatic PE -phred33 {input.fwd} {input.rev} {params.paired_R1} {params.unpaired_R1} {params.paired_R2} {params.unpaired_R2} {params.leading} {params.trailing} {params.sldWindow} {params.minLen}"

rule fastqc_trim:
    input:
        fwd = "data/trimming/out_paired_R1_{sample}.fastq.gz",
        rev = "data/trimming/out_paired_R2_{sample}.fastq.gz"
    output:
        fwd = "data/fastqc/trim/out_paired_R1_{sample}_fastqc.zip",
        rev = "data/fastqc/trim/out_paired_R2_{sample}_fastqc.zip"
    params:
        outdir = "data/fastqc/trim"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule multiqc_trimmed:
    input:
        expand("data/fastqc/trim/out_paired_{dir}_{sample}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"])
    output:
        "data/fastqc/trim/multiqc_report.html"
    params:
        workingDir = "data/fastqc/trim/*"
    shell:
        "multiqc -o data/fastqc/trim/ {params.workingDir}"

rule star:
    input:
        fwd = "data/trimming/out_paired_R1_{sample}.fastq.gz",
        rev = "data/trimming/out_paired_R2_{sample}.fastq.gz"
    output:
        "data/star/{sample}..star.ReadsPerGene.out.tab"
    params:
        outdir = "data/star/{sample}.star.",
        starDir = config["star"],
        gtf = config["gtf_file"]
    shell:
        "STAR --runThreadN 8 --genomeDir {params.starDir} --readFilesIn {input.fwd} {input.rev} --readFilesCommand zcat --sjdbGTFfile {params.gtf} --outFileNamePrefix {params.outdir} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode GeneCounts --twopassMode Basic --outReadsUnmapped Fastx"


rule compile_counts:
    input:
        "data/star/{sample}.star.ReadsPerGene.out.tab"
    output:
        "data/counts_table.txt"
    params:
        outdir = "data/counts_table.txt",
        indir = "data/star",
        inColumn = "4",
    shell:
        "Rscript scripts/compile_readouts.R -s {params.inColumn} -d {params.inDir} -o {params.outDir}"

rule ide:
    input:
        "data/counts_table.txt"
    output:
        "data/ide/filtered_counts.txt"
    params:
        gene_info_file = config["gene_info_file"]
    shell:
        "Rscript scripts/run_IDE.R {params.gene_info_file}"