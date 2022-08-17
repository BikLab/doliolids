#PBS -S /bin/bash
#PBS -N AssignASVs_DADA2
#PBS -q bik_q
#PBS -l nodes=1:ppn=1
#PBS -l mem=50gb
#PBS -l walltime=12:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUTABLE=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureTable_DADA2.qza
OUTPUTABLE=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureTable_DADA2.qzv
MAPPING=/scratch/tjp99569/fepe/qiime-files/16S-FEPE-samples-mapping-file-08-20-2020.txt
INPUTSEQ=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2.qza
OUTPUTSEQ=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureSeq_DADA2.qzv
INPUTSTATS=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_Denoise_Stats_DADA2.qza
OUTPUTSTATS=/home/tjp99569/fepe/dada2-results/

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc


#Visualize Feature Table Summary
qiime feature-table summarize \
  --i-table $INPUTABLE \
  --o-visualization $OUTPUTABLE \
  --m-sample-metadata-file $MAPPING \

#Visualize Rep-Seqs
qiime feature-table tabulate-seqs \
  --i-data $INPUTSEQ \
  --o-visualization $OUTPUTSEQ \

#Export Stats
  qiime tools export \
  --input-path $INPUTSTATS \
  --output-path $OUTPUTSTATS

