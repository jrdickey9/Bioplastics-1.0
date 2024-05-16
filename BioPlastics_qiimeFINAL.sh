#!/bin/bash

cd /Users/lab/Desktop/Jonathan_Dickey_LabMac/Bioplastics_Shurin/16S

#IMPORT DATA; completed on 04/23/24 
#qiime tools import \
#  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#  --input-path Seqs/ \
#  --output-path BioP-sequencesFINAL.qza

#VISUALIZE DATA IMPORT; completed on 4/16/24 
#qiime demux summarize \
# --i-data BioP-sequencesFINAL.qza \
#  --o-visualization BioP-sequencesFINAL.qzv

#FILTER BY Q SCORE IF USING DEBLUR (we aren't though, so note how is commented (#) out). 
#qiime quality-filter q-score \
# --i-demux BioP-sequencesR1.qza \
# --o-filtered-sequences BioP-filtered.qza \
# --o-filter-stats BioP-filter-stats.qza

#VISUALIZE FILTER STATS
#qiime metadata tabulate \
#  --m-input-file BioP-filter-stats.qza \
#  --o-visualization BioP-filter-stats.qzv

#PAIRED END DATA (Forward and Reverse) and DADA2; completed on 4/16/24 
#qiime dada2 denoise-paired \
#  --i-demultiplexed-seqs BioP-sequencesFINAL.qza \
#  --p-trunc-len-r 130 \
#  --p-trunc-len-f 138 \
#  --p-max-ee-f 2 \
#  --p-max-ee-r 3 \
#  --p-pooling-method "pseudo" \
#  --p-chimera-method "consensus" \
#  --o-representative-sequences BioP-paired-asv-sequencesFINAL.qza \
#  --o-table BioP-paired-asv-tableFINAL.qza \
#  --o-denoising-stats BioP-paired-denoising-statsFINAL.qza

##VISUALIZE PAIRED END DATA; completed on 4/16/24 
#qiime metadata tabulate \
#  --m-input-file BioP-paired-denoising-statsFINAL.qza \
#  --o-visualization BioP-paired-denoising-statsFINAL.qzv

qiime feature-table tabulate-seqs \
  --i-data BioP-paired-asv-sequencesFINAL.qza \
  --o-visualization BioP-paired-asv-sequencesFINAL.qzv

qiime feature-table summarize \
  --i-table BioP-paired-asv-tableFINAL.qza \
  --m-sample-metadata-file 16S_BioP_metadata.tsv \
  --o-visualization BioP-paired-asv-tableFINAL.qzv

#SILVA 138; TRAIN -- make sure its for 515-806
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads silva-138-99-seqs-515-806.qza \
#  --i-reference-taxonomy silva-138-99-tax-515-806.qza \
#  --o-classifier silva-138-99-515F-806R-nb-classifierFINAL.qza

#CLASSIFY AGAINST SILVA, PAIRED
#qiime feature-classifier classify-sklearn \
#  --i-classifier silva-138-99-515F-806R-nb-classifierFINAL.qza \
#  --i-reads BioP-paired-asv-sequencesFINAL.qza \
#  --output-dir BioP-paired-taxonomyFINAL

#VISUALIZE TAXONOMY PAIRED
#qiime metadata tabulate \
#  --m-input-file BioP-paired-taxonomyFINAL/classification.qza \
#  --o-visualization BioP-paired-taxonomyFINAL.qzv

#BUILD PHYLOGENY PAIRED 
#qiime phylogeny align-to-tree-mafft-fasttree \
# --i-sequences BioP-paired-asv-sequencesFINAL.qza \
# --o-alignment BioP-aligned-paired-seqsFINAL.qza \
# --o-masked-alignment BioP-masked-aligned-paired-seqsFINAL.qza \
# --o-tree BioP-unrooted-tree-paired-seqsFINAL.qza \
# --o-rooted-tree BioP-rooted-tree-paired-seqsFINAL.qza

#EXPORT PAIRED
#qiime tools export \
#  --input-path BioP-paired-asv-tableFINAL.qza \
#  --output-path BioP_paired_dadatableFINAL

qiime tools export \
  --input-path BioP-paired-taxonomyFINAL/classification.qza \
  --output-path BioP_paired_taxaFINAL
