# aTEA: active TE Annotation by Long Read RNA-seq (version 0.1)

## Introduction

aTEA (active TE annotation) is an analysis pipeline to identify active TE transcripts (including both autonomously expressed TE and TE-gene chimeric transcripts) from long-read RNA-seq data (e.g., PacBio Iso-seq or ONT cDNA-RNA sequencing). 

## Prerequisite

Linux system

python 3.8.17 or latest

Perl 5.32

Isoseq (https://github.com/pacificbiosciences/isoseq/)

BamTools (v2.5.1) (https://github.com/pezmaster31/bamtools)

minimap2 (v2.24) (https://github.com/lh3/minimap2)

TALON (v5.0) (https://github.com/mortazavilab/TALON)

CD-HIT (v4.8.1) (https://github.com/weizhongli/cdhit)

gffread (v0.12.7) (https://github.com/gpertea/gffread)

## Install and Run

1. Download the package to your local server
2. Unpack it using the command `tar -zxvf aTEA-0.1.tar.gz`
3. Run `bash aTEA.sh`

## Demo
We provide a demo data for testing, which is a CCS bam file containing 20,000 CCS reads. You can expect to get transcript identification results and TE classification results from this demo analysis. In addition, all relevant genomic data for analysis are also included. This demo run should be finished in 10 mins on a standard Linux computer with 10 cores.

### input
1. CCS bam file: `demo.ccs.bam`
2. zebrafish reference genome: `GRCz11.fa`
3. zebrafish genome annotation: `Danio_rerio.GRCz11.103.chr.gtf`
4. zebrafish TE annotation: `GRCz11.TE.fa.out`
5. Isoseq primer sequences: `primers.fasta`
6. BED12 file: `GRCz11.bed12`
7. TALON config file: `config_file_demo`

### output
1. transcript identification: `final_talon.gtf` (transcript identification results in gtf format) and `final_talon.transcripts.fa` (sequences file)
2. transcript quantification: `*.combined.filterlow` (reads count, TPM for each transcript)
3. transcript classification: `Demo.TE.annotation.list2` (TE-alone transcript with TE annotation); `Demo.TE-Gene.annotation.list2` (TE-gene transcript with TE annotation) and `Demo.Gene.annotation.list2` (TE-free gene transcripts)

## For Your Data
1. generate CCS reads from subreads bam: `ccs your.subreads.bam your.ccs.bam`
2. replace the name of demo file with your data: In aTEA.sh, replace the **demo** with **your**. If you have multiple samples/replicates to run, simply list all name in the loop, e.g., `sample1 sample2 sample3 ..`.
3. prepare config file for TALON: replace the path for the labeled SAM file in the file **config_file_demo**, e.g., `/your/home/workdir/labeled/your_sample_labeled.sam`.
4. run `bash aTEA.sh` or submit the job to a computing cluster.


