# aTEA: active TE Annotation by Long Read RNA-seq (version 0.1)

## Introduction

aTEA (active TE annotation) is an analysis pipeline to identify active TE transcripts (including both autonomously expressed TE and TE-gene chimeric transcripts) from long-read RNA-seq data (e.g., PacBio Iso-seq or ONT cDNA-RNA sequencing). 

## Prerequisite

Linux system

python 3.8.17 or latest

Perl 5.32

Isoseq https://github.com/pacificbiosciences/isoseq/

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

### input
1. CCS bam file: `demo.ccs.bam`
2. zebrafish reference genome: `GRCz11.fa`
3. zebrafish genome annotation: `Danio_rerio.GRCz11.103.chr.gtf`
4. zebrafish TE annotation: `GRCz11.TE.fa.out`
5. Isoseq primer sequences: `primers.fasta`
6. BED12 file: `GRCz11.bed12`
7. TALON config file: `config_file_demo`

### output
1. transcript identification: `final_talon.gtf` and `final_talon.transcripts.fa`
2. transcript quantification: `*.combined.filterlow`
3. transcript classification: `Demo.TE.annotation.list2`; `Demo.TE-Gene.annotation.list2` and `Demo.Gene.annotation.list2`

## For Your Data



