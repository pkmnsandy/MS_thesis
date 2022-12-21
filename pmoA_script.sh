#! /bin/bash



#################################################
#################################################
#####                                        ####
##### Amplicon analysis using QIIME2 & DADA2 ####
#####         Several libraries              ####
#####                                        ####
#################################################
#################################################


#change to working directory
cd /Volumes/Pkmnsandy_HDD/amplicon_seq

##preprocessing - reorganise R1 and R2 reads
##depending on what R1 starts with (primerFwd or primerREV)
#without adapter trimming or joining!

#find sequences that start with FWD primer in R1 files
#two output files:
#R1_FWD - true forward reads
#R2_FWD - true reverse reads

#find sequences that start with REV Primer in R1 files
#R1_REV - actually reverse reads
#R2_FWD - actually forward reads

# die ungezippten Dateien manuell in neuen Ordner kopieren, in dem die Auswertung erfolgen soll. Dann die working directory dahin ändern:


# Benennung der Datei prüfen!!!

for i in *R1*.fastq
do
  #echo $i
   var1=$(echo "$i" | (sed "s/\(^[[:digit:]]*_.*_R\)[1-2]\(_.*fastq$\)/\12\2/" ))
   var2=$(echo "$i" | (grep -o -a  "[[:upper:]]" | head -1))
   var3=$(echo "$i" | (grep -o -a "L00[1-2]"))
   
  # echo $var1
  # echo $var2
  # echo $var3
           #fwd primer (189f), rev primer (mb661r), normal direction
    cutadapt -g GGNGACTGGGACTTCTGG -A CCGGMGCAACGTCYTTACC --match-read-wildcards --action=none -o R1_${var2%}_FWD_${var3%}.fastq -p R2_${var2%}_FWD_${var3%}.fastq -e 0.25 --discard-untrimmed $i $var1  2>&1 | tee -a FWD.log
            #rev primer (189r), fwd primer (mb661f), normal direction
    cutadapt -g CCGGMGCAACGTCYTTACC -A GGNGACTGGGACTTCTGG --match-read-wildcards --action=none -o R1_${var2%}_REV_${var3%}.fastq -p R2_${var2%}_REV_${var3%}.fastq -e 0.25 --discard-untrimmed $i $var1  2>&1 | tee -a REV.log

done


#reorganise reads from above
#combine reads starting with fwd primer --> forward.fastq
#combine reads starting with rev primer --> reverse.fastq
#nomenclature needed for qiime import!

# im folgenden Befehl findet er jeweils L002 nicht, aber das liegt daran, dass im paired-end mode gearbeitet wird und es nur l001 Dateien gibt?
# Wiebke hatte hier am Ende der Schleife die Dateien in den libraries gezipped --> Problem: dann können die Sequenzen darin nicht gelesen werden und es entstehen leere Dateien... 

for i in *L001_R1*.fastq
do
  #echo $i
   #var1=$(echo "$i" | (sed "s/\(^[[:digit:]]*_.*_R\)[1-2]\(_.*fastq$\)/\12\2/" ))
   var2=$(echo "$i" | (grep -o -a  "[[:upper:]]" | head -1))
  # echo $var1
  # echo $var2
  echo "processing library $var2"
  
  mkdir lib$var2
  
  cat R1_${var2%}_FWD_L001.fastq R2_${var2%}_REV_L001.fastq  > lib$var2/forward.fastq
  
  cat R2_${var2%}_FWD_L001.fastq  R1_${var2%}_REV_L001.fastq  > lib$var2/reverse.fastq
  
  #mit zip entstehen leere Dateien...
  #gzip lib$var2/*.fastq
  
  done
  
  
  #### neue Version Skript Wiebke
  
  #making library list & Sub-SampleSheets

for i in *L001_R1*.fastq 
do
   var2=$(echo "$i" | (grep -o -a  "[[:upper:]]" | head -1))
   echo $var2 >> libs.txt
done

# hier muss Name von meinem Sample Sheet eingetragen werden?


map=metadata.txt


#while read i; do
    
    sed -n 1,2p metadata4.csv > Sample${i}.csv
    awk -v i="$i" -F , '{ if ($3 == i) { print } }' metadata4.csv >> Sample${i}.csv
    cat Sample${i}.csv | sed "s/,/\t/g"  > Sample${i}.tsv
    rm Sample${i}.csv
    
 done < libs.txt

# ich habe im obigen Befehl zusätzlich angegeben, dass meine Spalten Tabulatoren getrennt sind und die Zeile mit awk geändert:

while read i; do          
sed -n 1,2p $map > Sample${i}.csv     
awk -v i="$i" -F $'\t'  '{ if ($3 == i) { print } }' $map >> Sample${i}.csv 
cat Sample${i}.csv | sed "s/,/\t/g"  > Sample${i}.tsv
rm Sample${i}.csv
done < libs.txt 


#demultiplex with cutadapt
#at the beginning of each library, there will be an error message
#you can ignore this, it is because the first line of libs.txt
#is a header, which will of course not be found as a barcode ;)

 while read i; do
    
    
    echo "processing library $i \n"
    
    while read line; do
     varID=$(echo $line | grep -o -a "^S[[:digit:]]*"  )
    
     varBC=$(echo $line | grep -o -a "[ACGT]*" | head -1)
   
    

    echo "processing $varID \n"
    
    cutadapt -g $varBC -O 5 -o lib$i/${varID%}_R1.fastq -p lib$i/${varID%}_R2.fastq -e 0 --discard-untrimmed lib$i/forward.fastq lib$i/reverse.fastq  2>&1 | tee -a lib$i/demux.log 
    
    done < Sample$i.tsv
    
    
done < libs.txt

 

mkdir demux-all

#move demux files from lib directories to demux-all directory
while read i; do
    
    #echo $i
    mv lib${i}/S*.fastq demux-all 
        
done < libs.txt



cd demux-all

#rename files
#hier kommt ein hinweis/fehler, den man ignorieren kann
rename 's/(^S[0-9]*_)(R[1-2])/\1\1L001_\2_001/g' *.fastq

#die folgenden zwei Befehle habe ich ausprobiert, damit der Import zu QIIME2 funktionert. Tut er
#aber nicht. -n bewirkt, dass die Dateien erstmal nicht umbenannt werden, sondern nur die neuen 
#Namen in der Konsole angezeigt werden. Wenn die stimmen, kann -n entfernt werden
#rename -n 's/R1_001/R1_001_forward/' *fastq.gz
#rename -n 's/R2_001/R2_001_reverse/' *fastq.gz

#compress fastq files
gzip *.fastq

#### create manifest file, Version 2018.11, für neuere Version s.u. ####
#header line
#echo "sample-id,filename,direction" > MANIFEST.tsv

#write info for forward reads into MANIFEST
#for i in *_R1_*.fastq.gz
#do
    #echo $i 
    #var2=$(echo "$i" | (grep -o -a  "S[0-9]*" | head -1))
    #make corresponding R2-file name
    #var3=$(echo "$i" | (sed "s/\(^S[0-9]*_.*_R\)[1-2]\(_.*fastq.gz$\)/\12\2/" ))
    #echo $var3
    #echo $var2
    
    #echo $var2,$i,"forward" >> MANIFEST.tsv
    #echo $var2,$var3,"reverse" >> MANIFEST.tsv
    
#done

#cat MANIFEST.tsv | sed 's/\,/\t/g' > MANIFEST

# ODER:

#### create manifest file, Version 2019.7! ####
#header line
echo "sample-id,forward-absolute-filepath,reverse-absolute-filepath" > MANIFEST.tsv

#write info for forward reads into MANIFEST
for i in *_R1_*.fastq.gz
do
    
    var2=$(echo "$i" | (grep -o -a  "S[0-9]*" | head -1))
   
    var3=$(echo "$i" | (sed "s/\(^S[0-9]*_.*_R\)[1-2]\(_*.fastq.gz$\)/\12\2/" ))
    echo $var2,"/Volumes/Pkmnsandy_HDD/amp/demux-all/$i","/Volumes/Pkmnsandy_HDD/amp/demux-all/$var3" >> MANIFEST.tsv
done
#manually edit the manifest file to correct the path and the delimiter


#cat MANIFEST.tsv | sed $'s\_/\t/g' > MANIFEST.csv | sed $'s\_/\t/g'  > map_file.tsv - I manually converted csv to txt
#optional:
#rm MANIFEST.tsv

#einfach tsv in ubuntu aus manifest löschen

#write info for reverse reads in MANIFEST
#for i in *_R2_*.fastq.gz
#do
    #echo $i 
 #   var2=$(echo "$i" | (grep -o -a  "S[0-9]*" | head -1))
    #echo $var2
    
  #  echo $var2,$i,"reverse" >> MANIFEST
#done 

#change to directory above
cd ../
 
#activate qiime
source activate /home/shmoo/Illumina/bin/miniconda/miniconda3/envs/qiime2-2019.7
source activate /home/kato/miniconda3/envs/qiime2-2020.6


# the following steps can also be performed using the GUI (Q2-Studio)
# (except combining the libraries, further down)!
# Q2-Studio is installed and can be launched and exited with the commands below
#
# but I have not tried this out
# although it should be rather easy to use
# feel free to try it out and explore
#
# launch Q2-Studio
cd q2studio-2019.7.0
npm start

# leave Q2-Studio
cd ../


# from here, these are command line steps in QIIME2  
#import reads into qiime2
#if importable types need to be changed, check possibilities using:
#qiime tools import --show-importable-types

    
#import as whole artifact



# metadata.yml Datei fehlt? Online steht jedoch, dass die nicht selbst kreiert werden muss
#bei alter Version von qiime (2018.11) wäre dieser Befehl richtig:
#qiime tools import \input-format PairedEndFastqManifestPhred33V2
 #  --type 'SampleData[PairedEndSequencesWithQuality]' \
 #  --input-path demux-all \
 #  --output-path reads-demux.qza 
    source activate /home/shmoo/Illumina/bin/miniconda/miniconda3/envs/qiime2-2019.7
    
#es gibt neuere versionen! aufpassen mit welchem classifier man arbeitet, ansonsten funktioniert taxonomy nicht
#source activate /home/shmoo/miniconda3/envs/qiime2-2020.11


#bei neuer Version (2019.7):
# In this variant of the fastq manifest format, there must be forward and reverse read fastq.gz / fastq files for each sample ID. This format assumes that the PHRED offset used for the positional quality scores in all of the fastq.gz / fastq files is 33
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path demux-all/map_file.txt \
  --output-path qzav/reads-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2    
  

#####archaea, hat geklappt, wichtig hier das im manuscript 'absolute filepath' angegeben ist##############
#qiime tools import \
 # --type 'SampleData[SequencesWithQuality]' \
  #--input-path demux-all2/MANIFEST \
 # --output-path reads-demux.qza \
  #--input-format SingleEndFastqManifestPhred33V2
  

#check
qiime tools peek reads-demux.qza


#primer trimming
#sequences are 189f and mb661r
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences qzav/reads-demux.qza \
    --p-front-r CCGGMGCAACGTCYTTACC \
    --p-front-f GGNGACTGGGACTTCTGG \
    --p-error-rate 0.2\
    --o-trimmed-sequences qzav/trim-seqs.qza \
    --verbose
    
##################archaea##############################    

#primer trimming, nicht paired end
#sequences are 515f für archaea
#qiime cutadapt trim-single \
  #  --i-demultiplexed-sequences reads-demux.qza \
 #   --p-front GTGCCAGCMGCCGCGGTAA \
  #  --p-error-rate 0.15\
  # --o-trimmed-sequences trim-seqs.qza \
   # --verbose


    
    
#create visualisation
qiime demux summarize \
  --i-data qzav/trim-seqs.qza \
  --o-visualization qzav/trim-seqs.qzv

  
# export für ENA, demultiplexed trimmed
cd rawdata/archaea #path anpassen!
mkdir enasource activate /home/shmoo/Illumina/bin/miniconda/miniconda3/envs/qiime2-2019.7
qiime tools export --input-path trim-seqs.qza --output-path ena/ 





  
  
#visualise
#check for decreasing quality
qiime tools view trim-seqs.qzv


#denoising using dada2 - takes long to run! (ich habe es um 16:00 angemacht, am nächsten Morgen um 9:00 war der Befehl durchgelaufen)
#trun-len --> decide based on qualities from above
#amd change if required
#Daten Romy: ich habe den forward read bei 230 getrimmt und den reverse read bei 240
#run1 - denoise
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qzav/trim-seqs.qza \
    --p-trunc-len-f 240 \
    --p-trunc-len-r 240 \
    --p-n-threads 0 \
    --o-table qzav/denoise_1.qza \
    --o-representative-sequences qzav/rep-seqs_1.qza \
    --o-denoising-stats qzav/denoising-stats_1.qza \
    --verbose
    
#run2 - denoise
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qzav/trim-seqs.qza \
    --p-trunc-len-f 252 \
    --p-trunc-len-r 250 \
    --p-n-threads 0 \
    --o-table qzav/denoise_2.qza \
    --o-representative-sequences qzav/rep-seqs_2.qza \
    --o-denoising-stats qzav/denoising-stats_2.qza \
    --verbose
  
#run3 - denoise
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qzav/trim-seqs.qza \
    --p-trunc-len-f 260 \
    --p-trunc-len-r 250 \
    --p-n-threads 0 \
    --o-table qzav/denoise_3.qza \
    --o-representative-sequences qzav/rep-seqs_3.qza \
    --o-denoising-stats qzav/denoising-stats_3.qza \
    --verbose

#run4 - denoise
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qzav/trim-seqs.qza \
    --p-trunc-len-f 255 \
    --p-trunc-len-r 250 \
    --p-n-threads 0 \
    --o-table qzav/denoise_4.qza \
    --o-representative-sequences qzav/rep-seqs_4.qza \
    --o-denoising-stats qzav/denoising-stats_4.qza \
    --verbose
    
#qiime dada2 denoise-paired \
#    --i-demultiplexed-seqs trim-seqs.qza \
#    --p-trunc-len-f 186 \
#    --p-trunc-len-r 230 \
#    --p-n-threads 0 \
 #   --o-table denoise2.qza \
 #   --o-representative-sequences rep-seqs2.qza \
 #   --o-denoising-stats denoising-stats2.qza \
 #   --verbose
  
  
 ##############archaea#########single
 #qiime dada2 denoise-single \
  #  --i-demultiplexed-seqs trim-seqs.qza \
 ##   --p-trunc-len 200 \
  #  --p-n-threads 0 \
  #  --o-table denoise.qza \
  #  --o-representative-sequences rep-seqs.qza \
  #  --o-denoising-stats denoising-stats.qza \
  #  --verbose
  
##################create visualisations#####################################################
#first: overall denoising output ('OTU-table')
#incl. metadata
#change name of metadata-file (= SampleSheet)    
#if sample sheet is .csv change to .tsv
#change comma to tab in samplesheet
cat pmoA_wald-map.tsv | sed "s/,/\t/g"  > map_file.tsv


#run1
qiime feature-table summarize \
  --i-table qzav/denoise_1.qza \
  --o-visualization qzav/denoise_1.qzv \
  --m-sample-metadata-file metadata.txt

#run2
qiime feature-table summarize \
  --i-table qzav/denoise_2.qza \
  --o-visualization qzav/denoise_2.qzv \
  --m-sample-metadata-file metadata.txt

#run3
qiime feature-table summarize \
  --i-table qzav/denoise_3.qza \
  --o-visualization qzav/denoise_3.qzv \
  --m-sample-metadata-file metadata.txt

#run4
qiime feature-table summarize \
  --i-table qzav/denoise_4.qza \
  --o-visualization qzav/denoise_4.qzv \
  --m-sample-metadata-file metadata.txt

#visualise
qiime tools view denoise.qzv
  
  
#second: representative sequences
#run1
qiime feature-table tabulate-seqs \
  --i-data qzav/rep-seqs_1.qza \
  --o-visualization qzav/rep-seqs_1.qzv

#run2
qiime feature-table tabulate-seqs \
  --i-data qzav/rep-seqs_2.qza \
  --o-visualization qzav/rep-seqs_2.qzv

#run3
qiime feature-table tabulate-seqs \
  --i-data qzav/rep-seqs_3.qza \
  --o-visualization qzav/rep-seqs_3.qzv

#run4
qiime feature-table tabulate-seqs \
  --i-data qzav/rep-seqs_4.qza \
  --o-visualization qzav/rep-seqs_4.qzv
  

#visualise
qiime tools view rep-seqs.qzv

#third: denoising stats
#run1
qiime metadata tabulate \
  --m-input-file qzav/denoising-stats_1.qza \
  --o-visualization qzav/denoising-stats_1.qzv 

#run2
qiime metadata tabulate \
  --m-input-file qzav/denoising-stats_2.qza \
  --o-visualization qzav/denoising-stats_2.qzv

#run3
qiime metadata tabulate \
  --m-input-file qzav/denoising-stats_3.qza \
  --o-visualization qzav/denoising-stats_3.qzv

#run4
qiime metadata tabulate \
  --m-input-file qzav/denoising-stats_4.qza \
  --o-visualization qzav/denoising-stats_4.qzv 

#visualise
qiime tools view denoising-stats.qzv
 ############################################################# 

#if you have several runs, you should be able to combine them now, using the following commands
#NOTE: these are at the moment commented out
#you need to remove the # -signs 
#qiime feature-table merge \
#  --i-tables table-1.qza \
#  --i-tables table-2.qza \
#  --o-merged-table table.qza

#qiime feature-table merge-seqs \
#  --i-data rep-seqs-1.qza \
#  --i-data rep-seqs-2.qza \
#  --o-merged-data rep-seqs.qza

qiime feature-table filter-seqs \
    --i-data rep-seqs.qza \
    --m-metadata-file rep-seqs.qza \
    --p-where 'length(sequence) > 310' \
    --o-filtered-data filter_length-seqs.qza 
    
    #second: filtered representative sequences  
qiime feature-table tabulate-seqs \
  --i-data filter_length-seqs.qza \
  --o-visualization filter_length-seqs.qzv

#visualise
qiime tools view filter_length-seqs.qzv


#phylogenetics
#sequence alignment, mafft ist die software
qiime alignment mafft \
  --i-sequences output/rep-seqs_pmoA_filtered.qza \
  --o-alignment output/aligned-rep-seqs_pmoA.qza \
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
  --i-alignment output/aligned-rep-seqs_pmoA.qza \
  --o-masked-alignment output/masked_aligned-rep-seqs_pmoA.qza \
  --verbose

#build unrooted tree
qiime phylogeny fasttree \
  --i-alignment output/masked_aligned-rep-seqs_pmoA.qza \
  --o-tree output/unrooted-tree_pmoA.qza \
  --verbose
  
  
#rooted tree, using longest root
qiime phylogeny midpoint-root \
  --i-tree output/unrooted-tree_pmoA.qza \
  --o-rooted-tree output/rooted-tree_pmoA.qza \
  --verbose

#export rooted tree
qiime tools export \
  --input-path output/rooted-tree_pmoA.qza \
  --output-path exported-tree
  
#tax-assignment - takes long to run! (bei mir nur so 10 Minuten?)
#using QIIME built in naive Bayesian classifier - has tendency to over-classify reads
#see: https://usda-ars-gbru.github.io/Microbiome-workshop/tutorials/qiime2/ (accessed: 07 Aug 2019)
#first step: train classifier OR use pre-trained set
#here:  classifier model trained on the Silva 99% database trimmed to the V4 region

#TAXONOMIC ASSIGNMENT
qiime feature-classifier classify-sklearn \
  --i-classifier  taxonomy/pmoA_classifier_2.qza \
  --i-reads output/rep-seqs_pmoA_filtered.qza  \
  --p-confidence 0.8 \
  --o-classification output/taxonomy-pmoA_filtered_2.qza \
  --verbose
#working above

  
qiime feature-classifier classify-sklearn \
  --i-classifier  taxonomy/pmoa4rdp_qiime.tax \
  --i-reads rep-seqs.qza  \
  --p-reads-per-batch 2000 \
  --o-classification taxonomy_silva_138.qza 

 
#----------aktuelle phylogeny mit 138 silva
qiime feature-classifier classify-sklearn \
  --i-classifier  taxonomy/pmoa4rdp_qiime.tax \
  --i-reads rep-seqs.qza  \
  --p-reads-per-batch 2000 \
  --o-classification taxonomy_silva.qza 
  
 #-------------- phylogeny Max
 #Run the classifier - Consensus
#qiime feature-classifier classify-sklearn \
  --i-classifier taxonomy/pmoa4rdp_qiime.tax \
  --i-reads output/rep-seqs_pmoA_filtered.qza \
  --o-classification taxonomy138.qza \
  --p-n-jobs -2 \
  --p-confidence 0.75 \
  --verbose
  
qiime metadata tabulate \
  --m-input-file taxonomy_silva_138.qza \
  --o-visualization taxonomy.qzv
#view
qiime tools view taxonomy.qzv   


#training a classifier - Sandy
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path pmoA_tax/100.fasta \
  --output-path taxonomy/100.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path pmoA_tax/pmoA_cultivated_taxonomy.txt \
  --output-path taxonomy/pmoA_cultivated_taxonomy.qza

qiime feature-classifier extract-reads \
  --i-sequences taxonomy/100.qza \
  --p-f-primer GGNGACTGGGACTTCTGG \
  --p-r-primer CCGGMGCAACGTCYTTACC \
  --o-reads taxonomy/ref-seqs_pmoA_2.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads taxonomy/ref-seqs_pmoA_2.qza \
  --i-reference-taxonomy taxonomy/pmoA_cultivated_taxonomy.qza \
  --o-classifier taxonomy/pmoA_classifier_2.qza


#Filter out rare ASVs: features that appear less then 20 times and not in least 5 samples are thrown out
qiime feature-table filter-features \
   --i-table output/denoise.qza \
   --p-min-frequency 5 \
   --p-min-samples 2 \
   --o-filtered-table output/outputtable-filter.qza
   
qiime feature-table summarize \
  --i-table table-filter.qza \
  --o-visualization table-filter.qzv \
  --m-sample-metadata-file metadata.txt
  
qiime tools view table-filter.qzv

#Filter out contaminant and unclassified ASVs
qiime taxa filter-table \
   --i-table table-filter.qza \
   --i-taxonomy taxonomy_silva_138.qza \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-table table-filter2.qza
   
qiime feature-table summarize \
  --i-table table-filter2.qza \
  --o-visualization table-filter2.qzv \
  --m-sample-metadata-file map_file.tsv
 
 
 
 
 
  
#create visualisation
qiime metadata tabulate \
  --m-input-file taxonomy_silva_138.qza \
  --o-visualization silva/taxonomy.qzv

#view
qiime tools view taxonomy.qzv
  
#optional
#creates interactive barplots
#for a first glance at sample compositions
qiime taxa barplot \
  --i-table output/denoise_pmoA_filtered.qza \
  --i-taxonomy output/taxonomy-pmoA_filtered_2.qza \
  --m-metadata-file metadata.txt \
  --o-visualization output/pmoa_tax_barplot_2.qzv
  
 
#view 
qiime tools view silva/taxa-bar-plots2.qzv





  
################alpha und betadiversity###############################################
qiime diversity alpha-rarefaction \
  --i-table table-filter.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 3000 \
  --m-metadata-file metadata.txt \
  --o-visualization alpha-rarefaction.qzv
  
  
  qiime diversity alpha-rarefaction \
  --i-table temp2.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 6000 \
  --m-metadata-file metadata.txt \
  --o-visualization alpha-rarefaction2.qzv


#archaea
 # qiime diversity alpha-rarefaction \
 # --i-table denoise.qza \
 # --i-phylogeny rooted-tree.qza \
 # --p-max-depth 1200 \
#  --m-metadata-file map_file.tsv \
 # --o-visualization alpha-rarefaction.qzv
  
  
qiime tools view alpha-rarefaction.qzv
qiime tools view alpha-rarefaction2.qzv
  
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table denoise.qza \
  --p-sampling-depth 6800 \
  --m-metadata-file map_file.tsv \
  --output-dir core-metrics-results
  
  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-filter2.qza \
  --p-sampling-depth 6500 \
  --m-metadata-file map_file.tsv \
  --output-dir core-metrics-results4
  
  
  
  
qiime emperor plot \
  --i-pcoa core-metrics-results4/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file map_file.tsv \
  --o-visualization core-metrics-results4/unweighted-unifrac-emperor-elevation.qzv
  
qiime emperor plot \
  --i-pcoa core-metrics-results4/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file map_file.tsv \
  --o-visualization core-metrics-results4/weighted-unifrac-emperor-elevation.qzv
  
  #achtung 2,3 und 4 version existiert (3 und 4 sind die gefilterten, einmal nach Max (10 reads in min. 3 Proben), 4 mit 3 reads in mind 2 Proben)
qiime emperor plot \
  --i-pcoa core-metrics-results4/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file map_file.tsv \
  --o-visualization core-metrics-results4/unweighted-unifrac-emperor-elevation.qzv  
  
qiime tools view core-metrics-results4/weighted-unifrac-emperor-elevation.qzv


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results4/shannon_vector.qza \
  --m-metadata-file map_file.tsv \
  --o-visualization core-metrics-results4/shannon_group-significance.qzv
  
qiime tools view core-metrics-results4/shannon_group-significance.qzv



#export to .biom --> can be imported to R for analysis
#optional: use QIIME2 for further analysis (alpha, beta diversity, etc.)
#see Tutorials, examples given in the explanatory doc
 
 
 
 #first, export your data as a .biom
qiime tools export \
  --input-path output/output/otu-table_denoise_pmoA_min3k.qza \
  --output-path output/output/exported_data
  
qiime tools export \
  --input-path output/output/rep-seqs_pmoA.qza \
  --output-path output/output/exported_data

#convert to tsv  

  
biom convert \
  -i output/output/exported_data/feature-table.biom \
  -o output/output/exported_data/otu_table.txt \
  --to-tsv --header-key taxonomy 
  

#change header
var="OTUID"
sed -i "2s/^#OTU ID/$var/" "output/output/exported_data/otu_table.txt"
  
  
#then export taxonomy info
qiime tools export \
  --input-path output/output/taxonomy_pmoA_12453.qza \
  --output-path output/output/exported_data/taxonomy_pmoA.txt
 
# Manually change "feature ID" to "OTUID" die Raute ist hier wichtig!!!
#Change the first line of taxonomy.tsv (i.e. the header)
var="#OTUID \ttaxonomy \tconfidence" \
sed -i "1s/.*/$var/" output/output/exported_data/taxonomy.txt
 

 
 
qiime tools export \
  --input-path output/rooted-tree_pmoA.qza \
  --output-path exported-data


conda deactivate

############taxonomy info


cp exported-data4/taxonomy.tsv biom-taxonomy.tsv



biom add-metadata -i exported-data4/feature-table.biom -o table-with-taxonomy3.biom \
--observation-metadata-fp taxonomy3.tsv \
--sc-separated taxonomy
  

 
    
biom convert \
  -i table-with-taxonomy3.biom \
  -o ASV-final.txt --to-tsv --header-key taxonomy
  
  
biom convert -i table-with-taxonomy.biom -o ASV.txt --to-tsv --header-key taxonomy
  
#export to .biom --> can be imported to R for analysis
#optional: use QIIME2 for further analysis (alpha, beta diversity, etc.)
#see Tutorials, examples given in the explanatory doc

############export mit taxonomy####################### level ist gewünschte taxonmische Auflösung 2=phyla, 3=class, 6 =genera, 7 =strain
 qiime taxa collapse \
  --i-table table-filter2.qza \
  --i-taxonomy taxonomy_silva_138.qza \
  --p-level 6 \
  --output-dir taxtable6/
  
 cd taxtable6
 
 qiime tools export \
  --input-path collapsed_table.qza \
  --output-path exported-data
  
 biom convert \
  -i exported-data/feature-table.biom \
  -o exported-data/otu_table.txt \
  --to-tsv 
 
 
 
 
 #first, export your data as a .biom
qiime tools export \
  --input-path denoise.qza \
  --output-path exported-data

#convert to tsv  
biom convert \
  -i exported-data/feature-table.biom \
  -o exported-data/otu_table.txt \
  --to-tsv --header-key taxonomy

  biom convert \
  -i exported-data/feature-table.biom \
  -o exported-data/otu_table.txt \
  --to-tsv --header-key taxonomy
  
#change header
var="OTUID"
sed -i "2s/^#OTU ID/$var/" exported-data/otu_table.txt
  
  
#then export taxonomy info
qiime tools export \
  --input-path output/taxonomy-pmoA_filtered_2.qza \
  --output-path exported-data
 
# Manually change "feature ID" to "OTUID" 
#Change the first line of taxonomy.tsv (i.e. the header)
var="OTUID \t taxonomy \t confidence"
sed -i "1s/.*/$var/" exported-data/taxonomy.tsv 
 
qiime tools export \
  --input-path output/rooted-tree_pmoA.qza \
  --output-path exported-data


conda deactivate

#move to R









