#PBS -S /bin/bash
#PBS -N Removing_Primers
#PBS -q bik_q
#PBS -l nodes=1:ppn=5
#PBS -l mem=50gb
#PBS -l walltime=12:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUT=/home/tjp99569/fepe/dada2-results/01-FEPE_Samples_ArtifactFile.qza
OUTPUT=/home/tjp99569/fepe/dada2-results/01-FEPE_Samples_ArtifactFile_trimmed.qza

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#cutadapt parameters
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $INPUT \
  --p-cores 5 \
  --p-adapter-f TTACCGCGGCKGCTGRCACACAATTACCATA \
  --p-adapter-r ATTAGAWACCCBNGTAGTCCGGCTGGCTGACT \
  --p-error-rate 0.11 \
  --p-discard-untrimmed True \
  --o-trimmed-sequences $OUTPUT
