# Variant Calling

### Tools Used
- BWA 
- Samtools
- Bcftools

### Data Preprocessing
- First I downloaded the data and decompressed all the compressed files using `gunzip`

**Index reference genome**
```
bwa index chr7.fa
```
**Align the reads to ref_genome**
- Then I used BWA to align the raw sequencing reads from both tumor and normal samples to the human reference genome
```
bwa mem -t 8 chr7.fa R1.fastq R2.fastq > aligned.sam
```
**Convert SAM to BAM**
```
samtools view -S -b aligned.sam > aligned.bam
```
**Sort and index**
- After that I sorted and indexed the alignments using Samtools
```
samtools sort aligned.bam -o sorted.bam
samtools index sorted.bam
```
**Variant calling**
- Then I performed the variant calling using bcftools
```
bcftools mpileup -Ou -f chr7.fa sorted.bam | bcftools call -mv -Ob -o variants.bcf
```
**Convert BCF to VCF**
- Then I converted the output from BCF to VCF to use it as an input to VEP
```
bcftools view -Ov -o variants.vcf variants.bcf
```
**Annotate using VEP**
- Lastly I used VEP to annotate the variants
```
 vep -i variants.vcf -o result.vcf --species "human" --fasta chr7.fa --vcf --database
 ```
- The whole pipeline script can be seen in [pipeline.sh](pipeline.sh)

- After analyzing the output of VEP I found one pathogenic variant at position `2915243`. This variant is annotated as *stop gained*, which is pathogenic because it adds a stop codon where it is not supposed to be, which can lead to a loss of function of the protein



