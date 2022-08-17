#PBS -S /bin/bash
#PBS -N Export_FeatureTable
#PBS -q bik_q
#PBS -l nodes=1:ppn=1
#PBS -l mem=50gb
#PBS -l walltime=12:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUTABLE=/home/tjp99569/fepe/dada2-results/FEPE_trimmed_FeatureTable_DADA2.qza
OUTPUTABLE=/home/tjp99569/fepe/dada2-results/

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Export Feature table
  qiime tools export \
  --input-path $INPUTABLE \
  --output-path $OUTPUTABLE

