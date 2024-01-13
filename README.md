# aTEA: active TE Annotation by Long Read RNA-seq (version 0.1)

## Introduction

aTEA (active TE annotation) is an analysis pipeline to identify active TE transcripts (including both autonomously expressed TE and TE-gene chimeric transcripts) from long-read RNA-seq data (e.g., PacBio Iso-seq or ONT cDNA-RNA sequencing). 

## Prerequisite
### Software
Linux system

Python 3.8.17 or latest

Perl 5.32

Isoseq (https://github.com/pacificbiosciences/isoseq/)

BamTools (v2.5.1) (https://github.com/pezmaster31/bamtools)

minimap2 (v2.24) (https://github.com/lh3/minimap2)

TALON (v5.0) (https://github.com/mortazavilab/TALON)

CD-HIT (v4.8.1) (https://github.com/weizhongli/cdhit)

gffread (v0.12.7) (https://github.com/gpertea/gffread)

### Genomic data

Reference genome: `GRCz11.fa` (https://ftp.ensembl.org/pub/release-103/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz)

Gene annotation: `Danio_rerio.GRCz11.103.chr.gtf` (https://ftp.ensembl.org/pub/release-103/gtf/danio_rerio/Danio_rerio.GRCz11.103.chr.gtf)

TE annotation: `GRCz11.TE.fa.out` (https://figshare.com/articles/dataset/zebrafish_TE_annotation_RepeatMasker_/24991395)

BED12: `GRCz11.bed12` which can be obtained by using this command `paftools.js gff2bed  Danio_rerio.GRCz11.103.chr.gtf > GRCz11.bed12`

## Install and Run
This tool is not required to install on your local computer. To run it, 

1. Download the package `aTEA-0.1.tar.gz` to your local server
2. Unpack it using the command `tar -zxvf aTEA-0.1.tar.gz`
3. Put all required genomic data into the same folder
4. Make sure the required tools to be added into your PATH
5. Run `bash aTEA.sh`

## Demo
We provide a demo data for testing, which is a CCS bam file containing 20,000 CCS reads. You can expect to get transcript identification results and TE classification results from this demo analysis. In addition, several other relevant data for analysis are also included. This demo run should be finished in 10 mins on a standard Linux computer with 10 cores.

### Input
1. CCS bam file: `demo.ccs.bam`
2. Zebrafish reference genome: `GRCz11.fa`
3. Zebrafish genome annotation: `Danio_rerio.GRCz11.103.chr.gtf`
4. Zebrafish TE annotation: `GRCz11.TE.fa.out`
5. Isoseq primer sequences: `primers.fasta`
6. BED12 file: `GRCz11.bed12`
7. TALON config file: `config_file_demo`

### Output
1. Transcript identification: `final_talon.gtf` (transcript identification results in gtf format) and `final_talon.transcripts.fa` (sequences file)
2. Transcript quantification: `*.combined.filterlow` (reads count, TPM for each transcript)
3. Transcript classification: `Demo.TE.annotation.list2` (TE-alone transcript with TE annotation); `Demo.TE-Gene.annotation.list2` (TE-gene transcript with TE annotation) and `Demo.Gene.annotation.list2` (TE-free gene transcripts)

## For Your Data
1. Generate CCS reads from subreads bam: `ccs your.subreads.bam your.ccs.bam`
2. Prepare genomic data for your own species, including genome sequence, gene annotation, TE annotation, etc.
3. Replace the name of demo file with your data: In aTEA.sh, replace the **demo** with **your**. If you have multiple samples/replicates to run, simply list all name in the loop, e.g., `sample1 sample2 sample3 ..`.
4. Prepare config file for TALON: replace the path for the labeled SAM file in the file **config_file_demo**, e.g., `/your/home/workdir/labeled/your_sample_labeled.sam`.
5. Run `bash aTEA.sh` or submit the job to a computing cluster.


