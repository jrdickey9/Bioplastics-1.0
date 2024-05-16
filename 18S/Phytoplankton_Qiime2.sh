#!/bin/bash

#file path to dir with folder of seqs
cd /Users/lab/Desktop/Jonathan_Dickey_LabMac/Bioplastics_Shurin/18S

###Qiime2 lines begin here
##Import
#qiime tools import \
#--type 'SampleData[PairedEndSequencesWithQuality]' \
#--input-path seqs_18S \
#--input-format CasavaOneEightSingleLanePerSampleDirFmt \
#--output-path Phyto-demux-paired-endR1.qza

##Visualize
#qiime demux summarize \
#--i-data Phyto-demux-paired-endR1.qza \
#--o-visualization Phyto-demux-paired-endR1.qzv

#James et al. 2022 in Nature Communications does not quality filter

##denoise
#qiime dada2 denoise-paired \
#  --i-demultiplexed-seqs Phyto-demux-paired-endR1.qza \
#  --p-trunc-len-f 140 \
#  --p-trunc-len-r 138 \
#  --p-max-ee-f 2 \
#  --p-max-ee-r 2 \
#  --p-trunc-q 2 \
#  --p-pooling-method "independent" \
#  --p-chimera-method "consensus" \
#  --o-representative-sequences Phyto-asv-sequencesR1.qza \
#  --o-table Phyto-asv-tableR1.qza \
#  --o-denoising-stats Phyto-denoising-statsR1.qza

##View denoising
#qiime metadata tabulate \
#  --m-input-file Phyto-denoising-statsR1.qza \
#  --o-visualization Phyto-denoising-statsR1.qzv

##View asv sequences
#qiime feature-table tabulate-seqs \
#  --i-data Phyto-asv-sequencesR1.qza \
#  --o-visualization Phyto-asv-sequencesR1.qzv

## View feature table
#qiime feature-table summarize \
#  --i-table Phyto-asv-tableR1.qza \
#  --m-sample-metadata-file 18S_Phyto_metadata.tsv \
#  --o-visualization Phyto-asv-tableR1.qzv

## View feature table
#qiime feature-table summarize \
#  --i-table Phyto-asv-tableR1.qza \
#  --o-visualization Phyto-asv-tableR1.qzv

#export asv table
#qiime tools export \
#  --input-path Phyto-asv-tableR1.qza \
#  --output-path Phyto_Dada_table

#Seq alignment for fasttree
#qiime phylogeny align-to-tree-mafft-fasttree \
#  --i-sequences Phyto-asv-sequencesR1.qza \
#  --o-alignment aligned-Phyto-asv-sequences.qza \
#  --o-masked-alignment masked-aligned-Phyto-asv-sequences.qza \
#  --o-tree unrooted-tree-Phyto-asv-sequences.qza \
#  --o-rooted-tree rooted-tree-Phyto-asv-sequences.qza

#export phylogeny
#qiime tools export \
# --input-path rooted-tree-Phyto-asv-sequences.qza \
# --output-path rooted-tree-Phyto-asv-sequences

######################################################
####################Taxonomy w PR2####################
######################################################

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path pr2_version_5.0.0_SSU_mothur.fasta \
  --output-path pr2_v5.0.0.qza

#import taxonomy
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path pr2_version_5.0.0_SSU_mothur.tax \
  --output-path pr2_5.0.0_tax.qza

qiime feature-classifier extract-reads \
  --i-sequences pr2_v5.0.0.qza \
  --p-f-primer TTGTACACACCGCCC \
  --p-r-primer CCTTCYGCAGGTTCACCTAC \
  --o-reads pr2_v5.0.0_v9_extracts.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads pr2_v5.0.0_v9_extracts.qza \
  --i-reference-taxonomy pr2_5.0.0_tax.qza \
  --o-classifier pr2_v5.0.0_v9_classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier pr2_v5.0.0_v9_classifier.qza \
  --i-reads Phyto-asv-sequencesR1.qza \
  --o-classification merged_asv_tax_sklearn.qza

qiime tools export \
  --input-path merged_asv_tax_sklearn.qza \
  --output-path asv_tax_dir

