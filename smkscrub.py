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
        expand("scrubbed/{sample}.fastq.gz", sample=input_files)

rule cutadapt:
    input:
        fastq = f"{fastq_dir}/{{sample}}.fastq.gz"
    output:
        cleaned_fastq = "scrubbed/cleaned_{sample}.fastq.gz"
    params:
        adapter = "CTGTCTCTTATACACATCT"
    shell:
        """
        cutadapt -q 20 -a {params.adapter} --minimum-length 50 -o {output.cleaned_fastq} {input.fastq}
        """

rule map:
    input:
        reference = reference,
        fastq = "scrubbed/cleaned_{sample}.fastq.gz"
    output:
        bam = "scrubbed/mapped_{sample}.bam"
    shell:
        """
        minimap2 --split-prefix=tmp$$ -ax sr {input.reference} {input.fastq} | samtools view -bh | samtools sort -o {output.bam}
        """

rule non_host:
    input:
        bam = "scrubbed/mapped_{sample}.bam"
    output:
        fastq = "scrubbed/{sample}.fastq.gz"
    shell:
        """
        samtools fastq -F 3584 -f 77 {input.bam}  | gzip -c > scrubbed/R1_{wildcards.sample}.fastq.gz
        samtools fastq -F 3584 -f 141 {input.bam}  | gzip -c > scrubbed/R2_{wildcards.sample}.fastq.gz
        """
