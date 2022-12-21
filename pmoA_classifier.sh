#TAXONOMIC ASSIGNMENT
#training a classifier
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path pmoA_tax/100.fasta \
  --output-path taxonomy/100_tax/100.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path pmoA_tax/pmoA_cultivated_taxonomy.txt \
  --output-path taxonomy/100_tax/pmoA_cultivated_taxonomy.qza

qiime feature-classifier extract-reads \
  --i-sequences taxonomy/100_tax/100.qza \
  --p-f-primer GGNGACTGGGACTTCTGG \
  --p-r-primer CCGGMGCAACGTCYTTACC \
  --o-reads taxonomy/100_tax/100_pmoA_ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads taxonomy/100_tax/100_pmoA_ref-seqs.qza \
  --i-reference-taxonomy taxonomy/pmoA_cultivated_taxonomy.qza \
  --o-classifier taxonomy/100_tax/100_pmoA_ref-seqs_classifier.qza

#assigning taxonomy to rep-seqs
qiime feature-classifier classify-sklearn \
  --i-classifier  taxonomy/100_tax/100_pmoA_ref-seqs_classifier.qza \
  --i-reads output/output/rep-seqs_pmoA.qza  \
  --p-confidence 0.8 \
  --o-classification output/output/taxonomy_pmoA.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file output/output/taxonomy_pmoA_12453.qza \
  --o-visualization output/output/taxonomy_pmoA_12453.qzv

#view
qiime tools view taxonomy.qzv
  
#optional
#creates interactive barplots
#for a first glance at sample compositions
qiime taxa barplot \
  --i-table output/output/otu-table_denoise_pmoA_min3k_norare.qza \
  --i-taxonomy output/output/taxonomy_pmoA.qza \
  --m-metadata-file metadata/metadata.txt \
  --o-visualization output/output/pmoa_tax_barplot.qzv
  
 