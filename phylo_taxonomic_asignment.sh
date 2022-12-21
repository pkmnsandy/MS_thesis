
#phylogenetics
#sequence alignment, mafft ist die software
qiime alignment mafft \
  --i-sequences output/output/rep-seqs_pmoA.qza \
  --o-alignment output/output/rep-seqs_pmoA_aligned.qza \
  --verbose
  
# Hinweis am Ende des Befehls: 
#"Strategy:
 #FFT-NS-2 (Fast but rough)
 #Progressive method (guide trees were built 2 times.)

#If unsure which option to use, try 'mafft --auto input > output'.
#For more information, see 'mafft --help', 'mafft --man' and the mafft page.

#The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
#It tends to insert more gaps into gap-rich regions than previous versions.
#To disable this change, add the --leavegappyregion option."

#mask
qiime alignment mask \
  --i-alignment output/output/rep-seqs_pmoA_aligned.qza \
  --o-masked-alignment output/output/rep-seqs_pmoA_aligned_masked.qza \
  --verbose

#build unrooted tree
qiime phylogeny fasttree \
  --i-alignment output/output/rep-seqs_pmoA_aligned_masked.qza \
  --o-tree output/output/rep-seqs_pmoA_aligned_masked_unrooted.qza \
  --verbose
  
  
#rooted tree, using longest root
qiime phylogeny midpoint-root \
  --i-tree output/output/rep-seqs_pmoA_aligned_masked_unrooted.qza \
  --o-rooted-tree output/output/rep-seqs_pmoA_aligned_masked_unrooted_rooted.qza \
  --verbose

#export rooted tree
qiime tools export \
  --input-path output/output/rep-seqs_pmoA_aligned_masked_unrooted_rooted.qza \
  --output-path output/output/exported_data
  
#tax-assignment - takes long to run! (bei mir nur so 10 Minuten?)
#using QIIME built in naive Bayesian classifier - has tendency to over-classify reads
#see: https://usda-ars-gbru.github.io/Microbiome-workshop/tutorials/qiime2/ (accessed: 07 Aug 2019)
#first step: train classifier OR use pre-trained set
#here:  classifier model trained on the Silva 99% database trimmed to the V4 region

#all levels
#TAXONOMIC ASSIGNMENT
#training a classifier
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path taxonomy/pmoA_12500.fasta \
  --output-path taxonomy/tax_12453/fastseq_12500.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path taxonomy/taxonomy_12500.txt \
  --output-path taxonomy/tax_12453/taxonomy_12500.qza

qiime feature-classifier extract-reads \
  --i-sequences taxonomy/tax_12453/fastseq_12500.qza \
  --p-f-primer GGNGACTGGGACTTCTGG \
  --p-r-primer CCGGMGCAACGTCYTTACC \
  --o-reads taxonomy/tax_12453/pmoA_ref-seqs_12500.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads taxonomy/tax_12453/pmoA_ref-seqs_12500.qza \
  --i-reference-taxonomy taxonomy/tax_12453/taxonomy_12500.qza \
  --o-classifier taxonomy/tax_12453/pmoA_ref-seqs_12500_classifier.qza

#testing the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy/tax_12453/pmoA_ref-seqs_12500_classifier.qza \
  --i-reads output/output/filtered_rep-seqs.qza \
  --p-confidence 0.7 \
  --o-classification output/output/taxonomy_pmoA_12500_taxonomy_0.7.qza \
  --verbose


qiime metadata tabulate \
  --m-input-file output/output/taxonomy_pmoA_12500.qza \
  --o-visualization output/output/taxonomy_pmoA_12500.qzv
  
#optional
#creates interactive barplots
#for a first glance at sample compositions
qiime taxa barplot \
  --i-table output/output/otu-table_denoise_pmoA_470_min3k.qza \
  --i-taxonomy output/output/taxonomy_pmoA_12500_taxonomy.qza \
  --m-metadata-file metadata/metadata.txt \
  --o-visualization output/output/pmoa_tax_barplot.qzv




#classifier_update
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path taxonomy/pmoA_seq.fasta \
  --output-path taxonomy/tax_12453/fastseq_12500.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path taxonomy/pmoA_tax_file.txt \
  --output-path taxonomy/tax_12453/taxonomy_10k_new.qza

#site-specific classifier
qiime feature-classifier extract-reads \
  --i-sequences taxonomy/tax_12453/fastseq_12500.qza \
  --p-f-primer GGNGACTGGGACTTCTGG \
  --p-r-primer CCGGMGCAACGTCYTTACC \
  --o-reads taxonomy/tax_12453/pmoA_ref-seqs_10k-cleaned-repremoved_new_site.qza

#making the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  taxonomy/tax_12453/pmoA_ref-seqs_10k-cleaned-repremoved_new_site.qza \
  --i-reference-taxonomy taxonomy/tax_12453/taxonomy_10k_new.qza \
  --o-classifier taxonomy/tax_12453/pmoA_ref-seqs_10k-cleaned-repremoved-classifier_new.qza

#running the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy/tax_12453/pmoA_taxonomy_classifier.qza \
  --i-reads output/output/rep-seqs_pmoA.qza \
  --p-confidence 0.7 \
  --o-classification output/output/taxonomy_pmoA_10k_taxonomy_new_conf-70.qza \
  --verbose

qiime metadata tabulate \
  --m-input-file output/output/taxonomy_pmoA_10k_taxonomy_new_conf-70.qza \
  --o-visualization output/output/taxonomy_pmoA_10k_taxonomy_new_conf-70.qzv

#building a classifier using RESCRIPT

#import taxonomy and fasta file
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path taxonomy/pmoA_seq.fasta \
  --output-path taxonomy/tax_12453/fastseq_12500.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path taxonomy/pmoA_tax_file.txt \
  --output-path taxonomy/tax_12453/taxonomy_10k_new.qza

  #filter taxonomy file to match sequence file
  qiime rescript filter-taxa \
    --i-taxonomy taxonomy/tax_12453/taxonomy_10k_new.qza \
    --m-ids-to-keep-file taxonomy/tax_12453/fastseq_12500.qza \
    --o-filtered-taxonomy taxonomy/tax_12453/taxonomy_10k_new_keep.qza
 
qiime rescript evaluate-taxonomy \
    --i-taxonomies taxonomy/tax_12453/taxonomy_10k_new_keep.qza \
    --o-taxonomy-stats taxonomy/tax_12453/taxonomy_10k_new_keep_eval.qzv

qiime metadata tabulate \
    --m-input-file taxonomy/tax_12453/taxonomy_10k_new_keep.qza \
    --o-visualization taxonomy/tax_12453/taxonomy_10k_new_keep.qzv
  
qiime rescript evaluate-seqs \
    --i-sequences taxonomy/tax_12453/fastseq_12500.qza \
    --p-kmer-lengths 32 16 8 \
    --o-visualization taxonomy/tax_12453/fastseq_12500_eval.qzv

#Build and evaluate our classifier
qiime rescript evaluate-fit-classifier \
    --i-sequences taxonomy/tax_12453/fastseq_12500.qza \
    --i-taxonomy taxonomy/tax_12453/taxonomy_10k_new_keep.qza \
    --p-n-jobs 2 \
    --o-classifier taxonomy/tax_12453/pmoA_taxonomy_classifier.qza \
    --o-evaluation taxonomy/tax_12453/pmoA_taxonomy_evaluation.qzv \
    --o-observed-taxonomy taxonomy/tax_12453/pmoA_taxonomy_predicted_taxonomy.qza
  
qiime rescript evaluate-taxonomy \
  --i-taxonomies taxonomy/tax_12453/taxonomy_10k_new_keep.qza taxonomy/tax_12453/pmoA_taxonomy_predicted_taxonomy.qza\
  --p-labels ref-taxonomy predicted-taxonomy \
  --o-taxonomy-stats taxonomy/tax_12453/pmoA_ref-taxonomy-evaluation.qzv

#create phyloseq objects
mkdir -p Phyloseq

#exporting ASV table
qiime tools export --input-path output/output/otu-table_denoise_pmoA_min3k.qza --output-path Phyloseq/

#exporting taxonomy
qiime tools export --input-path output/output/taxonomy_pmoA_10k_taxonomy_new_conf-70.qza --output-path Phyloseq/

#exporting unrooted tree
qiime tools export --input-path output/output/rep-seqs_pmoA_aligned_masked_unrooted.qza --output-path Phyloseq/

#reformatting taxonomy file
sed 's/Feature ID/#OTUID/' Phyloseq/taxonomy.tsv | sed 's/Taxon/taxonomy/' | sed 's/Confidence/confidence/' > Phyloseq/biom-taxonomy.tsv

#merging taxonomy with biom file
biom add-metadata \
    -i Phyloseq/feature-table.biom \
    -o Phyloseq/taxa_table.biom \
    --observation-metadata-fp Phyloseq/biom-taxonomy.tsv \
    --sc-separated taxonomy

qiime tools export table.qza --output-dir exported
biom convert -i Phyloseq/feature-table.biom -o feature-table.tsv --to-tsv