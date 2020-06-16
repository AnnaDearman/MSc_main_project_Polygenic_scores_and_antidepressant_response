#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=50G   # Request 50GB RAM
#$ -j y

# Run IC to produce diagnostic plots

printf "ic=[location of software]ic.pl" > Config2.conf

source ./Config2.conf

cd [file path]

for chr in $(seq 1 22);
do
  gunzip chr"$chr".dose.vcf.cut.gz
done

perl $ic -d [file path] -r [file path]1000GP_Phase3_combined.legend -g -p EUR [file path]

exit
