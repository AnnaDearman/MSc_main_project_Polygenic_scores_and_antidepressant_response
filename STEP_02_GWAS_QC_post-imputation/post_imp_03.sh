#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=6G   # Request 6GB RAM
#$ -j y

cd [file path]
 
for chnum in $(seq 1 22);
  do

  echo " **************************************************************************  " 
  echo "                       Converting vcf file to plink"
  echo " **************************************************************************  " 

  plink --vcf [file path]chr"$chnum".dose.vcf --make-bed --out [project]_"$chnum"
  
  echo " **************************************************************************  " 
  echo "                     Attempting merge to identify failures"
  echo " **************************************************************************  " 

  plink --bfile [project]_"$chnum" --bmerge [project]_"$chnum" --merge-mode 6
  
  echo " **************************************************************************  " 
  echo "                           Excluding failed SNPs"
  echo " **************************************************************************  " 
  
  plink --bfile [project]_"$chnum" --exclude plink.diff --make-bed --out [project]no_fails_"$chnum"
  
  echo " **************************************************************************  " 
  echo "                           Identifying duplicates" 
  echo " **************************************************************************  " 

  plink --bfile [project]no_fails_"$chnum" --list-duplicate-vars
  
  echo " **************************************************************************  " 
  echo "                           Excluding duplicate SNPs"
  echo " **************************************************************************  " 

  plink --bfile [project]no_fails_"$chnum" --exclude plink.dupvar --make-bed --out [project]no_fails_no_dups_"$chnum"
  
  echo " **************************************************************************  "   
  echo "                     Filter by impuation quality (r2=0.4)"
  echo " **************************************************************************  " 

  plink --bfile [project]no_fails_no_dups_"$chnum" --qual-scores [file path]chr"$chnum".info 7 1 1 --qual-threshold 0.4 --make-bed --out [file path]/chr"$chnum"

done
