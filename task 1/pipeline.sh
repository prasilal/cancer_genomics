#!/bin/bash

REF_GEN=$1
R1=$2
R2=$3


# Align data to the human reference genome
bwa mem -t 8 $REF_GEN $R1 $R2 > aligned.sam

# Convert SAM to BAM format
samtools view -S -b aligned.sam > aligned.bam

# Sort
samtools sort aligned.bam -o aligned_sorted.bam

# Index
samtools index aligned_sorted.bam

# Variant calling
bcftools mpileup -Ou -f $REF_GEN aligned_sorted.bam | bcftools call -mv -Ob -o variants.bcf

# Covert BCF to VCF
bcftools view -Ov -o variants.vcf variants.bcf

#Annotate using VEP
./vep_install/ensembl-vep/vep -i variants.vcf -o result.vcf --species "human" --fasta $REF_GEN --vcf --database
