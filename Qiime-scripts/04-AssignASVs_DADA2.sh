#PBS -S /bin/bash
#PBS -N AssignASVs_DADA2
#PBS -q bik_q
#PBS -l nodes=1:ppn=5
#PBS -l mem=50gb
#PBS -l walltime=12:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUT=/home/tjp99569/fepe/dada2-results/01-FEPE_Samples_ArtifactFile_trimmed.qza
TABLE=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureTable_DADA2.qza
SEQS=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2.qza
STATS=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_Denoise_Stats_DADA2.qza

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc


#Dada2 parameters
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $INPUT \
  --p-trim-left-f 0 \
  --p-trunc-len-f 237 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 253 \
  --p-n-threads 5 \
  --o-table $TABLE \
  --o-representative-sequences $SEQS \
  --o-denoising-stats $STATS \
  --verbose
