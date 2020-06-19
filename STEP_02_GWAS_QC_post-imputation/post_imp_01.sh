#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -j y

# Unpack your imputed genotype files 

printf "vcfparse=[location]vcfparse.pl" > Config.conf

source ./Config.conf

password="[password given by Michigan imputation server]"

cd [file location]

for chr in $(seq 1 22);
  do
    7za x chr_"$chr".zip -p$password
  done

for chr in $(seq 1 22);
  do
    gunzip chr"$chr".info.gz
    gunzip chr"$chr".dose.vcf.gz
  done

perl $vcfparse -d [file path] -o [file path] -g

exit
