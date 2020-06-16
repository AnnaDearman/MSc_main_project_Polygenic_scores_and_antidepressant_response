#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=20G   # Request 20GB RAM
#$ -j y

cd [file path]

echo " **************************************************************************  " 
echo "                  Making a merge list for chromosomes 2-22"
echo " **************************************************************************  " 

rm merge.list

for i in $(seq 2 22);

do
    printf "chr"$i"\n" >> merge.list 

done

echo " **************************************************************************  " 
echo "               Merging chromosomes 2-22 with chromosome 1 data"
echo " **************************************************************************  " 

plink --bfile chr1 --merge-list merge.list --make-bed --out [file path]/project_1kg_imputed_0.4_pre_naming_conversion

echo " ************************************************************************** "
echo "              Converting bim file to annovar input file format"
echo " ************************************************************************** "

gawk '{if(match($0,/:([0-9]+):/,m))print $1 "\t" m[1] "\t" m[1] "\t" $6 "\t" $5 "\t" $2}' project_1kg_imputed_0.4_pre_naming_conversion.bim > variants.txt

echo " ************************************************************************** "
echo "                       Getting rsIDs using Annovar"
echo " ************************************************************************** "

# Download SNP database

cd software/annovar/

# Convert the IDs in our data files using Annovar

perl annotate_variation.pl --filter --dbtype snp138 --build hg19 [file path]variants.txt --out [file path]variants_rsids.txt humandb/

cd ../../

echo " **************************************************************************  "
echo "  Getting a list of the rsIDs and removing any duplicates"
echo " **************************************************************************  "

cat variants_rsids.txt.hg19_snp138_dropped | awk '{print $8 "\t" $2}' | sort -u -k2 | sort -k1 > SNPs.txt

echo " **************************************************************************** "
echo "                   Converting SNP IDs to rsIDs"
echo " ************************************************************************** "

plink \
  --bfile project_1kg_imputed_0.4_pre_naming_conversion \
  --update-name SNPs.txt \
  --make-bed \
  --out project_1kg_imputed_0.4

echo " **************************************************************************  " 
echo "                Getting a list of SNPs with rsIDs" 
echo " **************************************************************************  " 

fgrep rs project_1kg_imputed_0.4.bim | awk '{print $2}' | awk '!seen[$1]++'  >  project_1kg_imputed_0.4.rsids

echo " **************************************************************************  " 
echo "                  Limiting data to variants with rs IDs"
echo " **************************************************************************  " 

plink \
 --bfile project_1kg_imputed_0.4 \
 --extract project_1kg_imputed_0.4.rsids \
 --make-bed \
 --out project_1kg_imputed_0.4_snps

echo " **************************************************************************  " 
echo "                           Filtering common SNPs"
echo " **************************************************************************  " 

plink \
 --bfile project_1kg_imputed_0.4_snps \
 --maf 0.01 \
 --make-bed \
 --out project_1kg_imputed_0.4_snps_maf0.01.common

echo " **************************************************************************  " 
echo "      Run Iterative Missingness script to remove SNPs, then samples,"
echo "      at increasingly high missingness cutoffs "
echo " **************************************************************************  " 

printf "root=[file path]/project_1kg_imputed_0.4_snps_maf0.01
plink=[binaries file path]plink" > Config.conf

sh [file path]/Iterative_Missingness.sh 90 99 1

source ./Config.conf

# N.B. output files are $root.filtered

echo " **************************************************************************  " 
echo "                     Getting sex info to put back in to the plink files"
echo " **************************************************************************  " 

awk '{print $1, $2, $5}' [file path]/project_kept_snp99_ind95_b37.het_cleaned.fam > project.update.sex

echo " **************************************************************************  " 
echo "                           Updating sex"
echo " **************************************************************************  " 

plink \
 --bfile $root.filtered \
 --make-bed \
 --update-sex project.update.sex \
 --out  project_1kg_imputed_0.4_snps_maf0.01_geno0.01_mind0.01_sex

exit 
