#!/bin/bash

cd ../myresults/hw
pwd

# Building a feauter table, which gives the length of each sample sequences
# This will help me determine what the sample depth should be used when choosing reads to be used in diversity calculations
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /data/project_data/16S/pyc_manifest

# Creating a seqencing table gives sequence length statistics as well as the exac sequence of each feature ID
# Sequences can be BLASTed on the NCBI database to look up if a sequences matches that of a known taxon
qiime feature-table tabulate-seqs \
  --i-data rep-seqs1.qza \
  --o-visualization rep-seqs1.qzv

# denoising-stats1.qsv gives statistics about the filtering process such as what percent of the each paired-end passed the filter, what percentage got merged, etc 
qiime metadata tabulate \
  --m-input-file denoising-stats1.qza \
  --o-visualization denoising-stats1.qzv

# The command below will build a phylogenetic tree of the identified taxa by using the sequencing table created above (rep-seqs1.qza) 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs1.qza \
  --o-alignment aligned-rep-seqs1.qza \
  --o-masked-alignment masked-aligned-rep-seqs1.qza \
  --o-tree unrooted-tree1.qza \
  --o-rooted-tree rooted-tree1.qza

# Now using our rooted phyogenetic tree, we can calculate a who suite of diversity metrics
# These results will go into the newly-created directory core-metrics-results_HW
# The sampling depth was choosen based on the shortest read, which was 31,623 based on results from table.qsv
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree1.qza \
  --i-table table.qza \
  --p-sampling-depth 31000 \
  --m-metadata-file /data/project_data/16S/pyc_manifest \
  --output-dir core-metrics-results_HW

# Creating a visualization file to see how Faith Phylogenetic Diversity (a measure of alpha-diversity) varies across different metadata variables, such as animal health and site status 
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_HW/faith_pd_vector.qza \
  --m-metadata-file /data/project_data/16S/pyc_manifest \
  --o-visualization core-metrics-results_HW/faith-pd-group-significance.qzv

# Creating a visualization file to see how species evenness (also a measure of alpha-diversity) varies across different metadata variables
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results_HW/evenness_vector.qza \
  --m-metadata-file /data/project_data/16S/pyc_manifest \
  --o-visualization core-metrics-results_HW/evenness-group-significance.qzv

# Creating a rarefaction plot in order to visualize the effect of our trimming on various diversity metrics by metadata factors 
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree1.qza \
  --p-max-depth 100000 \
  --m-metadata-file /data/project_data/16S/pyc_manifest \
  --o-visualization alpha-rarefaction

