#PBS -S /bin/bash
#PBS -N Taxonomy_Boxplots
#PBS -q bik_q
#PBS -l nodes=1:ppn=5
#PBS -l mem=50gb
#PBS -l walltime=12:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUT_TAX=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2_BLAST_Taxonomy.qza
OUTPUT_TAX=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2_BLAST_Taxonomy.qzv
TABLE=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureTable_DADA2.qza
MAPPING=/scratch/tjp99569/fepe/qiime-files/16S-FEPE-samples-mapping-file-08-20-2020.txt
BOXPLOT=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2_BLAST_Taxonomy_BoxPlots.qzv

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Visualize Taxonomy
qiime metadata tabulate \
  --m-input-file $INPUT_TAX \
  --o-visualization $OUTPUT_TAX \

#Generate boxplots
qiime taxa barplot \
  --i-table $TABLE \
  --i-taxonomy $INPUT_TAX \
  --m-metadata-file $MAPPING \
  --o-visualization $BOXPLOT \
  --verbose
