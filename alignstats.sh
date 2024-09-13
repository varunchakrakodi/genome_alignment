#!/bin/bash

# Output file name
output_csv="alignment_info.csv"

# Header for the CSV file
echo "Filename,Number of Reads Aligned,Coverage (%),Minimum Depth,Maximum Depth,Average Depth" > "$output_csv"

# Iterate over each file.bam
for bam_file in *.bam; do
    # Filename without extension
    filename=$(basename -- "$bam_file")
    filename="${filename%.*}"
    
    # Number of reads aligned
    num_aligned=$(samtools view -c -F 4 "$bam_file")
    
    # Depth file corresponding to the BAM file
    depth_file="${filename}_depth.txt"
    
    # Calculate coverage using awk
    coverage=$(awk '$3 > 0 {covered++} END {print covered/NR * 100}' "$depth_file")
    
    # Calculate minimum, maximum, and average depth using awk
    depth_stats=$(awk '{sum+=$3; if (NR==1 || $3<min) min=$3; if (NR==1 || $3>max) max=$3} END {printf "%.2f,%d,%d", sum/NR, min, max}' "$depth_file")
    IFS=',' read avg_depth min_depth max_depth <<< "$depth_stats"
    
    # Write to CSV
    echo "$filename,$num_aligned,$coverage,$min_depth,$max_depth,$avg_depth" >> "$output_csv"
done

echo "CSV file created: $output_csv"
