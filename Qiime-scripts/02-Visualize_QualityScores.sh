#PBS -S /bin/bash
#PBS -N Visualize_Quality
#PBS -q bik_q
#PBS -l nodes=1:ppn=1
#PBS -l mem=50gb
#PBS -l walltime=12:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUT=/home/tjp99569/fepe/dada2-results/01-FEPE_Samples_ArtifactFile_trimmed.qza
OUTPUT=/home/tjp99569/fepe/dada2-results/01-FEPE_Samples_ArtifactFile_trimmed.qzv

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Visualize Quality Scores
qiime demux summarize \
  --i-data $INPUT \
  --o-visualization $OUTPUT

