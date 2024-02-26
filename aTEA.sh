
###########################################################
###  active TE Annotation pipeline                      ###
###  Bo Li                                              ### 
###  01-11-2024                                         ###
###  email: libowyj@med.umich.edu                       ### 
###########################################################

#############################################################
## This pipeline is used for TE transcript identification from Long read RNA-seq data

## install the following tools before runing it.

#To install isoseq, using conda install -c bioconda isoseq
#To install bamtools, using conda install -c bioconda bamtools
#To install minimap2, using conda install -c bioconda minimap2
#To install TALON, see https://github.com/mortazavilab/TALON
#To install CD-HIT, see https://github.com/weizhongli/cdhit
#To install gffread, see https://github.com/gpertea/gffread

## other data required
# 1. genomic data: GRCz11.fa (zebrafish reference genome) (downloaded from https://ftp.ensembl.org/pub/release-103/fasta/danio_rerio/dna/)
# 2. genome annotation Danio_rerio.GRCz11.103.chr.gtf (downloaded from https://ftp.ensembl.org/pub/release-103/gtf/danio_rerio/) 
# 3. bed12 files: using the following command line to genereate a bed12 file
#    paftools.js gff2bed  Danio_rerio.GRCz11.103.chr.gtf > GRCz11.bed12
# 4. Iso-seq primer sequences: primers.fasta
# 5. TALON configure file: config_file_final





#############################################################
                      ## Main steps ##

## run ccs to extract ccs reads from subreads

#for id in demo 
#do
#	ccs cell0_1.subreads.bam $id.ccs.bam --all --log-level INFO --report-json cell0_1.report.json --log-file cell0_1.ccs.log --report-file cell0_1.report.txt --metrics-json cell0_1.zmw_metrics.json.gz -j 28
#done

## separate all ccs reads into full length and non-full length
export PATH="/home/libowyj/miniconda3/envs/pacbio/bin/:/home/libowyj/miniconda3/envs/bamtools/bin/:$PATH"

for bam in demo
do
	lima $bam.ccs.bam primers.fasta $bam.fl.bam --isoseq --peek-guess --store-unbarcoded
    isoseq refine $bam.fl.NEB_5p--NEB_Clontech_3p.bam primers.fasta $bam.flnc.bam --require-polya
    bamtools convert -in $bam.flnc.bam -format fasta -out $bam.flnc.fa
    bamtools convert -in  $bam.fl.unbarcoded.bam -format fasta -out $bam.lq.fa
    trim_isoseq_polyA -i $bam.lq.fa -G > $bam.lq.clean.fa 2>$bam.lq.lima.trimA.log 
done

## prepare reversed complement sequences for lq data because these sequences are not oriented withion adapter and polyA detection.

perl reverse_complement.pl

## run minimap2 with both full-length and non-full-length data

minimap2 -t 10 --MD -ax splice:hq -uf --secondary=no --junc-bed GRCz11.bed12 GRCz11.fa demo.flnc.fa  > demo.ccs.hq.sam
minimap2 -t 10 --MD -ax splice:hq -uf --secondary=no --junc-bed GRCz11.bed12 GRCz11.fa demo.lq.clean.fa demo.lq.clean.rc.fa  > demo.ccs.lq.sam

## run talon for transcript identification

mkdir -p labeled

for sam in demo
do
	talon_label_reads --f $sam.ccs.hq.sam --g GRCz11.fa --t 20 --deleteTmp --o labeled/$sam.hq
	talon_label_reads --f $sam.ccs.lq.sam --g GRCz11.fa --t 20 --deleteTmp --o labeled/$sam.lq
done

talon_initialize_database --f Danio_rerio.GRCz11.103.chr.gtf --g GRCz11 --a GRCz11_annot --o talon
talon --f config_file_demo --db talon.db --build GRCz11 --o zebrafish -t 16
talon_filter_transcripts --db talon.db -a GRCz11_annot  --maxFracA 0.5 --minCount 1 --allowGenomic --minDatasets=1 --o filtered_transcripts_minC1minD1.csv
talon_create_GTF --whitelist filtered_transcripts_minC1minD1.csv --db talon.db -a GRCz11_annot --build GRCz11 --o filtered_transcripts
talon_abundance --whitelist filtered_transcripts_minC1minD1.csv --db talon.db --build GRCz11 -a GRCz11_annot --o zebrafish_all

## filter the talon results 

perl filter_talon.pl --talon_tsv zebrafish_all_talon_abundance_filtered.tsv --known_trans_reads 1 --novel_trans_reads 2 --novel_gene_reads 2

## classify all transcripts into TE-alone, TE-gene and gene

bash TE_classification.sh

## prepare results

mkdir -p minimap2_dir
mkdir -p TALON_dir
mkdir -p TE_classification_dir
mkdir -p main_results
mv *.list2 final_talon.* *.TE.sum main_results/
mv *.fl.* *.flnc.* *.lq.* *.hq.* minimap2_dir/
mv talon_final_id.list *.csv *talon*.tsv* *_QC.log filtered_transcripts_talon.gtf talon.db labeled talon_tmp TALON_dir/
mv *.TE_overlap_with_transcript.gtf Demo_transcript_output Demo_transcript_output.format *_transcript_output_exon.tsv *_transcript_output_TE.tsv TE_classification_dir/








