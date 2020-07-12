#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=20G   # Request 20GB RAM
#$ -j y

echo "***********************************************************************"
echo "                 Polygenic scoring - height   "
echo "***********************************************************************"

PRSice_linux \
 --A1 Allele1 \
 --A2 Allele2 \
 --beta \
 --pvalue p \
 --stat b \
 --snp MarkerName \
 --base ~[file path]/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt \
 --fastscore \
 --all-score \
 --no-regress \
 --print-snp \
 --bar-levels 0.0001,0.01,0.05,0.1,0.5,1 \
 --target ~[file path and file name] \
 --out ~[file path and file name]
