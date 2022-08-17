#PBS -S /bin/bash
#PBS -N ImportFastqfepe
#PBS -q bik_q
#PBS -l nodes=1:ppn=1
#PBS -l mem=50gb
#PBS -l walltime=12:00:00
#PBS -M tjp99569@uga.edu
#PBS -m abe
#PBS -j oe

#Path Variables
INPUT=/scratch/tjp99569/fepe/qiime-files/16S-FEPE-fastq-file-manifest-08-20-2020.txt
OUTPUT=/home/tjp99569/fepe/dada2-results/01-FEPE_Samples_ArtifactFile.qza

#Activate QIIME2/2020.6
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

module load QIIME2/2020.6

echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

#Script
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $INPUT \
  --output-path $OUTPUT \
  --input-format PairedEndFastqManifestPhred33V2
