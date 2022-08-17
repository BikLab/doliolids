#PBS -S /bin/bash
#PBS -N Assign_Taxonomy
#PBS -q bik_q
#PBS -l nodes=1:ppn=5
#PBS -l mem=50gb
#PBS -l walltime=72:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUT=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2.qza
REFSEQ=/scratch/tjp99569/fepe/database/SILVA138-nr99_sequences_16S.qza
REFTAX=/scratch/tjp99569/fepe/database/SILVA138-nr99_taxonomy_16S.qza
OUTPUT=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2_BLAST_Taxonomy.qza

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Blast+ parameters
qiime feature-classifier classify-consensus-blast \
--i-query $INPUT \
--i-reference-reads $REFSEQ \
--i-reference-taxonomy $REFTAX \
--p-perc-identity 0.8 \
--o-classification $OUTPUT \
--verbose 

