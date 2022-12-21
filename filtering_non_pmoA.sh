#filtering out non-pmoA features from denoise artifact
qiime feature-table filter-features \
--i-table output/denoise.qza \
--m-metadata-file filter.txt \
--p-exclude-ids \
--o-filtered-table output/output/otu-table_denoise_pmoA.qza

qiime feature-table summarize \
  --i-table output/output/otu-table_denoise_pmoA.qza \
  --o-visualization output/output/otu-table_denoise_pmoA.qzv \
  --m-sample-metadata-file metadata/metadata.txt

#filtering out non-pmoA sequences from rep-seqs artifact
qiime feature-table filter-seqs \
--i-data output/rep-seqs.qza \
--m-metadata-file filter.txt \
--p-exclude-ids \
--o-filtered-data output/output/rep-seqs_pmoA.qza

qiime feature-table tabulate-seqs \
  --i-data output/output/rep-seqs_pmoA.qza \
  --o-visualization output/output/rep-seqs_pmoA.qzv


#filtering out sequences with less than 470 bp length
qiime feature-table filter-seqs \
    --i-data output/output/rep-seqs_pmoA.qza \
    --m-metadata-file output/output/rep-seqs_pmoA.qza \
    --p-where 'length(sequence) > 470' \
    --o-filtered-data output/output/rep-seqs_pmoA_400.qza 

qiime feature-table tabulate-seqs \
  --i-data output/output/rep-seqs_pmoA_400.qza \
  --o-visualization output/output/rep-seqs_pmoA_400.qzv

#filtering out features with less than 470 bp length
qiime feature-table filter-features \
--i-table output/output/otu-table_denoise_pmoA.qza \
--m-metadata-file features-to-keep.txt \
--o-filtered-table output/output/otu-table_denoise_pmoA_470.qza

qiime feature-table summarize \
  --i-table output/output/otu-table_denoise_pmoA.qza \
  --o-visualization output/output/otu-table_denoise_pmoA_470.qzv \
  --m-sample-metadata-file metadata/metadata.txt

#filtering out samples with less than 3000 reads
  qiime feature-table filter-samples \
  --i-table output/output/otu-table_denoise_pmoA_470.qza \
  --p-min-frequency 3000 \
  --o-filtered-table output/output/otu-table_denoise_pmoA_470_min3k.qza

  qiime feature-table summarize \
  --i-table output/output/otu-table_denoise_pmoA_470_min3k.qza \
  --o-visualization output/output/otu-table_denoise_pmoA_470_min3k.qzv \
  --m-sample-metadata-file metadata/metadata.txt

#filtering out rare ASVs
qiime feature-table filter-features \
   --i-table output/output/otu-table_denoise_pmoA_470_min3k.qza \
   --p-min-frequency 5 \
   --p-min-samples 2 \
   --o-filtered-table output/output/otu-table_denoise_pmoA_470_min3k_norare.qza
   

#generating ASV table
qiime feature-table summarize \
  --i-table output/output/otu-table_denoise_pmoA_470_min3k.qza \
  --o-visualization output/output/otu-table_denoise_pmoA_470_min3k_ASV_table.qzv \
  --m-sample-metadata-file metadata/metadata.txt


#filtering rep-seqs.qza for taxonomic assignment

qiime feature-table filter-seqs \
--i-data output/output/rep-seqs_pmoA_400.qza \
--m-metadata-file output/output/features-to-keep.txt \
--o-filtered-data output/output/filtered_rep-seqs.qza

qiime feature-table tabulate-seqs \
  --i-data output/output/filtered_rep-seqs.qza \
  --o-visualization output/output/filtered_rep-seqs.qzv