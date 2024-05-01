# Cancer Genomics Data Analysis

### Tools Used
- BWA 
- Samtools

### Data Preprocessing
- First I downloaded the data and decompressed all the compressed files using `gunzip`

### Indexing of the Human Reference Genome
```
bwa index hg19.fa
```

### Alignment to the Human Reference Genome
- Then I used BWA to align the raw sequencing reads from both tumor and normal samples to the human reference genome
```
bwa mem -t 8 hg19.fa tu.r1.fa tu.r2.fa > tu_aligned.sam
bwa mem -t 8 hg19.fa wt.r1.fa wt.r2.fa > wt_aligned.sam
```

### Coverting SAM to BAM
```
samtools view -S -b tu_aligned.sam > tu_aligned.bam
samtools view -S -b wt_aligned.sam > wt_aligned.bam
```

### Sorting and Indexing Alignments
- After that I sorted and indexed the alignments using Samtools
- tumor:
```
samtools sort tu_aligned.bam -o tu_aligned_sorted.bam
samtools index tu_aligned_sorted.bam
```
- wild-type:
```
samtools sort wt_aligned.bam -o wt_aligned_sorted.bam
samtools index wt_aligned_sorted.bam
```

### Step 3: Subsetting to Region of Interest 
- Then I extracted the region of interest (chrX:20000000-40000000) and subsequently sorted and indexed the subsetted files
- tumor:
```
samtools view -b tu_aligned_sorted.bam chrX:20000000-40000000 > tu_subset.bam
samtools sort tu_subset.bam -o tu_subset_sorted.bam
samtools index tu_subset_sorted.bam
```
- wild-type:
```
samtools view -b wt_aligned_sorted.bam chrX:20000000-40000000 > wt_subset.bam
samtools sort wt_subset.bam -o wt_subset_sorted.bam
samtools index wt_subset_sorted.bam
```

### Step 4: Generating Read-Depth Plot
- Calculated read depth at each position in the genome using Samtools.
- Generated a read-depth plot to visualize genomic alterations.


