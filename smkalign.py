import os

# Load the configuration file
configfile: "config.yaml"

# Extract variables from the config file
reference = config["reference"]
fastq_dir = config["fastq_dir"]

# Detect all .fastq.gz files in the specified directory
input_files = glob_wildcards(f"{fastq_dir}/{{sample}}.fastq.gz").sample

rule all:
    input:
        expand("results/{sample}.bam.bai", sample=input_files),
        expand("results/{sample}.fasta", sample=input_files),
        expand("results/{sample}_depth.txt", sample=input_files),
        expand("results/{sample}.vcf",sample=input_files)
        

rule cutadapt:
    input:
        fastq = f"{fastq_dir}/{{sample}}.fastq.gz"
    output:
        cleaned_fastq = "results/cleaned_{sample}.fastq.gz"
    params:
        adapter = "CTGTCTCTTATACACATCT"
    shell:
        """
        cutadapt -q 20 -a {params.adapter} --minimum-length 50 -o {output.cleaned_fastq} {input.fastq}
        """

rule bwa_mem2:
    input:
        reference = reference,
        cleaned_fastq = "results/cleaned_{sample}.fastq.gz"
    output:
        sam = "results/{sample}.sam"
    shell:
        """
        bwa-mem2 mem -t 16 {input.reference} {input.cleaned_fastq} > {output.sam}
        """

rule samtools_sort:
    input:
        sam = "results/{sample}.sam"
    output:
        bam = "results/{sample}.bam"
    shell:
        """
        samtools view -bS {input.sam} | samtools sort -o {output.bam}
        """

rule samtools_index:
    input:
        bam = "results/{sample}.bam"
    output:
        index = "results/{sample}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

rule samtools_consensus:
    input:
        reference = reference,
        bam = "results/{sample}.bam"
    output:
        fasta = "results/{sample}.fasta"
    shell:
        """
        samtools consensus -f fasta {input.bam} -d 5 -o {output.fasta}
        """

rule samtools_depth:
    input:
        bam = "results/{sample}.bam"
    output:
        depth_txt = "results/{sample}_depth.txt"
    shell:
        """
        samtools depth {input.bam} > {output.depth_txt}
        """

rule vcf:
    input:
        reference = reference,
        bam = "results/{sample}.bam"
    output:
        vcf = "results/{sample}.vcf"
    shell:
        """
        bcftools mpileup -Ou -f {input.reference} {input.bam} | bcftools call --ploidy 1 -mv -Ov | bcftools filter -i 'QUAL >= 50.0' -o {output.vcf}
        """
