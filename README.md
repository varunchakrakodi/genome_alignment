There are two scripts provided here.

**1. smkalign.py**

This is a Snakemake file that can be used to process the .fastq.gz files and align against a reference sequence. For this to work, concat R1 and R2 reads for each sample and put them in a folder. The method uses bwa-mem2 based alignment. Hence, the reference sequence needs to be indexed. Also requires a config.yaml file (Template is provided).

**Command:**
It is a good idea to run DAG (Directed Acyclic Graph) before executing the snakemake program to ensure everything is in place. (This step is optional). Also gives a visualisation of the expected process.

snakemake --snakefile smkalign.py --configfile config.yaml --dag | dot -Tpng > dag.png

For executing snakemake

snakemake --snakefile smkalign.py --configfile config.yaml --cores n

Depending on the system configuration you are using, filesystems might cache the status of files. Further, there can be a delay between the completion of a job and the visibility of its output files. This can cause Snakemake to prematurely fail, thereby wrongly assuming the files were not generated.
Adding a latency wait gives the filesystem time to update and make the new files visible. In such cases Use --latency-wait sec

$PATH = "/path/to/working_dir"

$PATH/results containing all the outputs will be created in the working directory for verification and further processing if needed.

Unfortunately, samtools consensus calling creates a header as >reference in each fasta file. This leads to issues in downstream processing if you run cat *.fasta > Samples.fasta

Hence it is recommended to change headers in all .fasta files using the following command structure before further processing.

>sed -i 's/^>.*/>file/' file.fasta

**Dependencies:**
1. Snakemake (https://snakemake.readthedocs.io/en/stable/)
2. Cutadapt ( https://github.com/marcelm/cutadapt )
3. samtools (https://github.com/samtools/samtools/releases/)
4. bwa-mem2 (https://github.com/bwa-mem2/bwa-mem2)
5. If running DAG, requires Graphviz (https://github.com/graphp/graphviz)

.bai file can be viewed using Tablet Alignment viewer (https://ics.hutton.ac.uk/tablet/)

**2. alignstats.sh**

Sometimes it is useful to have a deeper look into the Alignment statistics. The bash script looks into .bam and file_depth.txt files and generates alignment_info.csv file containing data on Number of Reads, Aligned	Coverage (%), Minimum Depth, Maximum Depth and Average Depth

**Command:**
./alignstats.sh $PATH/results

Best works for Illumina Short-reads

