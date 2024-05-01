#!/bin/bash

REF_GEN=$1
R1=$2
R2=$3
PREFIX=$4


# Align data to the human reference genome
bwa mem -t 8 $REF_GEN $R1 $R2 > $PREFIX"_aligned.sam"

# Convert SAM to BAM format
samtools view -S -b $PREFIX"_aligned.sam" > $PREFIX"_aligned.bam"

# Sort
samtools sort $PREFIX"_aligned.bam" -o $PREFIX"_aligned_sorted.bam"

# Index
samtools index $PREFIX"_aligned_sorted.bam"

# Subset to the Region of Interest
samtools view -b $PREFIX"_aligned_sorted.bam" chrX:20000000-40000000 > $PREFIX"_subset.bam"
samtools sort $PREFIX"_subset.bam" -o $PREFIX"_subset_sorted.bam"
samtools index $PREFIX"_subset_sorted.bam"

# Generate a read-depth file
samtools depth $PREFIX"_subset_sorted.bam" > $PREFIX"_read_depth.txt"

