#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -j y

echo "*************************************************************************"
echo
echo " Plot histogram of differences in MAF between GENDEP and 1000G "
echo
echo "*************************************************************************"
echo

R --file=PlotFrequencies.R

mv Diff_MAF_hist.pdf [new location]

exit
