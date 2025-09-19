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
### This PoC requires [Linux](https://www.linux.org/) ([wsl](https://ubuntu.com/desktop/wsl) has been used)

### Set-Up
- Step 1. Install Java
1. Install Java (version >= 17)
```
sdk install java 17.0.10-tem
```
2. Confirm installation
```
java -version
```

- Step 2: [Install Nextflow](https://www.nextflow.io/docs/latest/install.html)
1. Download Nextflow
```
curl -s https://get.sdkman.io | bash
```
2. Make Nextflow executable
```
chmod +x nextflow
```
3. Move Nextflow into an executable path. For example:
```
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
```
4. Confirm Nextflow is installed correctly
```
nextflow info
```

### Running the project
- Step 1: Install the required dependencies/libraries/tools
```
sudo apt update
sudo apt install -y fastqc bwa samtools bcftools
```

- Step 2: Index the hoax chromosome file (reference.fa)
```
bwa index data/reference.fa
```

- Step 3: Run the pipeline
```
nextflow run main.nf
```

## Expected outputs
Expected outputs in results/:  
\t qc_report.html → FastQC report.  
\t aligned.bam → aligned reads.  
\t variants.vcf → called variants.  