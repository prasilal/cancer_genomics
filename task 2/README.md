# Cancer Genomics Data Analysis Pipeline

## Tools Used
- BWA 
- Samtools

## Pipeline Steps

### Step 1: Alignment to the Human Reference Genome
- Used BWA to align the raw sequencing reads from both tumor and normal samples to the human reference genome.

### Step 2: Sorting and Indexing Alignments
- Sorted the alignments using Samtools.
- Indexed the sorted BAM files for quick retrieval.

### Step 3: Subsetting to Region of Interest (Optional)
- Subsetting the aligned reads to a specific region of interest in the genome.

### Step 4: Generating Read-Depth Plot
- Calculated read depth at each position in the genome using Samtools.
- Generated a read-depth plot to visualize genomic alterations.


