#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -j y

source ./Config.conf

echo "*************************************************************************"
echo
echo " Run EIGENSOFTplus to generate new plots "
echo
echo "*************************************************************************"
echo

# Although Eigensoft has already been run, I was unable to generate
# population stratification outliers and a screeplot, so I have made
# this script for that sole purpose

cd [location_of_data]

# N.B. You need to look at the screeplot before settling on "numoutlierevec"

R --vanilla --slave --args stem=[filename] numoutevec=100 numoutlieriter=5 nsnpldregress=0 numgamma=100 gamplot="YES" heatplot="YES" numoutlierevec=3 outliersigmathresh=6 kmeans="NO" snpweightout="YES" ESOFTdir=[Location_of_software] < [Location_of_software]