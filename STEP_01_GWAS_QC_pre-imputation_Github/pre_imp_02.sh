#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=10G   # Request 10GB RAM
#$ -j y

echo "*************************************************************************"
echo
echo " Use checkVCF to ensure that the VCF files are valid "
echo
echo "*************************************************************************"
echo

cd [data file path]

cp [file path]/human_g1k_v37.fasta .

cp [file path]/human_g1k_v37.fasta.fai .

module unload anaconda3/5.2.0

module load python/2.7.15

for chr in -chr1 -chr2 -chr3 -chr4 -chr5 -chr6 -chr7 -chr8 -chr9 -chr10 -chr11 -chr12 -chr13 -chr14 -chr15 -chr16 -chr17 -chr18 -chr19 -chr20 -chr21 -chr22

do

vcf="filename${chr}.vcf.gz"

python [location]checkVCF.py -r human_g1k_v37.fasta -o out $vcf

dir="checkVCF${chr}/"

rm -r $dir

mkdir $dir

mv out.check* $dir

done

module unload python/2.7.15

module load anaconda3/5.2.0

exit

