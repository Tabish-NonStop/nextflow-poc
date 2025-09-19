## Welcome to this PoC demonstrating (one of the many) usecases of [Nextflow](https://www.nextflow.io/) in bioinformatics


# Docs
This PoC is a mini bioinformatics pipeline with three steps:

- QC (Quality Check)  
Takes your FASTQ file (sequencing reads).  
Runs fastqc to check read quality.  
Produces a simple qc_report.html.  
  
- ALIGN (Read Mapping)  
Takes your reads and the reference genome (reference.fa).  
Runs bwa mem to map reads onto the reference.  
Produces an aligned.bam file (where each read aligns on the genome).  

- VARIANT_CALL (Variant Detection)  
Takes the alignment (aligned.bam) and reference.  
Uses samtools + bcftools to look for differences (mutations/variants) between your reads and the reference.  
Produces a variants.vcf file.  
  

In super simple words:  
- Step 1: “Are my DNA reads good?”  
- Step 2: “Where do my reads fit in the reference genome?”  
- Step 3: “What genetic changes do my reads have compared to the reference?”  


# How to run this project?

- 1. Install Nextflow (Linux is preferred)
```
curl -s https://get.sdkman.io | bash
```