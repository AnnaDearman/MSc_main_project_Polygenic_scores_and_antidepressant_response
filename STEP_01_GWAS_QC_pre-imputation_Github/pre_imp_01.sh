#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=60G   # Request 60GB RAM
#$ -j y

echo "*************************************************************************"
echo
echo " Create and run config file "
echo
echo "*************************************************************************"
echo

printf "root=[data_file_path]
plink=[plink_binaries_filepath]" > Config2.conf

source ./Config2.conf

echo "*************************************************************************"
echo
echo " Create a frequency file "
echo
echo "*************************************************************************"
echo

$plink --freq --bfile $root.het_cleaned --out $root.freq

echo "*************************************************************************"
echo
echo " Execute check-bim script "
echo
echo "*************************************************************************"
echo

perl HRC-1000G-check-bim-NoReadKey.pl -b $root.het_cleaned.bim -f $root.freq.frq -r /data/scratch/bt19624/gendep/1000GP_Phase3_combined.legend -g -p EUR

sh Run-plink.sh

mv *_cleaned* [new file path]

echo "*************************************************************************"
echo
echo " Compress vcf files "
echo
echo "*************************************************************************"
echo

cd [new file path]

for chr in -chr1 -chr2 -chr3 -chr4 -chr5 -chr6 -chr7 -chr8 -chr9 -chr10 -chr11 -chr12 -chr13 -chr14 -chr15 -chr16 -chr17 -chr18 -chr19 -chr20 -chr21 -chr22

do

vcf="filename${chr}.vcf"

bgzip $vcf

done

mv *cleaned-updated*.log updated/
mv *cleaned-updated*.bim updated/
mv *cleaned-updated*.fam updated/
mv *cleaned-updated*.bed updated/

exit
