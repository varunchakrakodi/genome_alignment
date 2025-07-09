import os
import glob # Don't forget to import glob

##Take details from config yaml file
configfile: "config.yaml"
reference = config["reference"]
fastq_dir = config["fastq_dir"]
input_files = glob.glob(f"{fastq_dir}/*.fastq.gz") # Use glob.glob for more direct file listing
input_samples = [os.path.basename(f).replace(".fastq.gz", "") for f in input_files] # Extract sample names

rule all:
    input:
        expand("results/{sample}.bam.bai", sample=input_samples),
        expand("results/{sample}.fasta", sample=input_samples),
        expand("results/{sample}_depth.txt", sample=input_samples),
        expand("results/{sample}.vcf.gz.csi", sample=input_samples) # Changed to reflect indexed VCF as final
        
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
    threads: 16 # Explicitly define threads
    shell:
        """
        bwa-mem2 mem -t {threads} {input.reference} {input.cleaned_fastq} > {output.sam}
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
        bai = "results/{sample}.bam.bai"
    shell:
        """
        samtools index {input.bam} -o {output.bai}
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
        
rule vcf_index:
    input:
        vcf= "results/{sample}.vcf"
    output:
        gz = "results/{sample}.vcf.gz", # Corrected path
        csi = "results/{sample}.vcf.gz.csi" # Corrected path
    shell:
        """
        bgzip {input.vcf}
        bcftools index {output.gz} -o {output.csi}
        """

rule pysam:
    input:
        bam = "results/{sample}.bam",
        bam_index = "results/{sample}.bam.bai"
    output:
        bed = "results/{sample}.bed"
    shell:
        """
        python3 /home/neurovirology/Apps/varun/nocov.py -d 5 {input.bam} {output.bed}
        """

rule consensus:
    input:
        reference = reference,
        bed = "results/{sample}.bed",
        gz = "results/{sample}.vcf.gz" # This expects the gzipped VCF
    output:
        fasta = "results/{sample}.fasta"
    shell:
        """
        bcftools consensus -f {input.reference} -H 1 -m {input.bed} {input.gz} -o {output.fasta}
        """
