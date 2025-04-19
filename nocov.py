import pysam
import argparse

###Function to call. Keep default depth as 1, but use -d 5 in command.
def find_no_coverage_regions(bam_file, output_bed, min_depth=1):
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")  # Open BAM file for reading.
        header = samfile.header

        with open(output_bed, "w") as bed_file:
            for chrom in header.references:
                chrom_length = header.get_reference_length(chrom)
                coverage = samfile.count_coverage(chrom, start=0, end=chrom_length)

                start = 0
                no_coverage_start = None

                for pos in range(chrom_length):
                    depth = sum(coverage[i][pos] for i in range(len(coverage))) # sum the coverage for all reads at this position.
                    if depth < min_depth:
                        if no_coverage_start is None:
                            no_coverage_start = pos
                    else:
                        if no_coverage_start is not None:
                            bed_file.write(f"{chrom}\t{no_coverage_start}\t{pos}\n")
                            no_coverage_start = None
                #handle the case where the no coverage region extends to the end of the chromosome.
                if no_coverage_start is not None:
                    bed_file.write(f"{chrom}\t{no_coverage_start}\t{chrom_length}\n")

    except FileNotFoundError:
        print(f"Error: BAM file '{bam_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if 'samfile' in locals() and samfile:
            samfile.close()

def main():
    parser = argparse.ArgumentParser(description="Identify regions with no coverage in a BAM file.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("output_bed", help="Path to the output BED file.")
    parser.add_argument("-d", "--min_depth", type=int, default=1, help="Minimum coverage depth (default: 1).")

    args = parser.parse_args()

    find_no_coverage_regions(args.bam_file, args.output_bed, args.min_depth)

if __name__ == "__main__":
    main()
