#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=5G   # Request 5GB RAM
#$ -j y

source ./Config.conf

echo "**************************************************************"
echo
echo " Get 1KG genotype and phenotype data "
echo
echo "**************************************************************"
echo

# Obtain 1KG Phase 1 data from PLINK2 website

cd [location_of_data_large_capacity_drive]

wget https://www.dropbox.com/s/k9ptc4kep9hmvz5/1kg_phase1_all.tar.gz
tar -xvzf 1kg_phase1_all.tar.gz

# Obtain 1KG Population info from 1KG

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped

grep -f <(awk '{print $2}' 1kg_phase1_all.fam) <(awk 'NR > 1 {print 0, $2, $7}' 20130606_g1k.ped) > 1KG_Phenos.txt

# Changes 1KG phenotypes (populations) to numbers starting at 3
# When merged with my data, my phenotpyes should all be 2.
# Case-control data can have phenotypes 1 and 2, which fits in with the below system

sed -i -e 's/ASW/3/g' -e 's/CEU/4/g' -e 's/CHB/5/g' -e 's/CHS/6/g' -e 's/CLM/7/g' -e 's/FIN/8/g' -e 's/GBR/10/g' -e 's/IBS/11/g' -e 's/JPT/12/g' -e 's/LWK/13/g' -e 's/MXL/14/g' -e 's/PUR/15/g' -e 's/TSI/16/g' -e 's/YRI/17/g' 1KG_Phenos.txt

echo "***************************************************************"
echo
echo " Restrict 1KG data and root data to common rsID SNPs "
echo
echo "***************************************************************"
echo

# Limit files to SNPs with rs IDs

fgrep rs $root.IBD_cleaned.bim > $root.IBD_cleaned.rsids.txt

# Get rs ID variant names

awk '{print $2}' $root.IBD_cleaned.rsids.txt > $root.IBD_cleaned.rsid_names.txt

# Extract rs IDs from root

$plink \
--bfile $root.IBD_cleaned \
--extract $root.IBD_cleaned.rsid_names.txt \
--chr 1-22 \
--make-bed \
--out $root.IBD_cleaned.rsids.autosomal

# Extract rs IDs from 1KG (and add phenotypes)

$plink \
--bfile 1kg_phase1_all \
--extract $root.IBD_cleaned.rsid_names.txt \
--pheno 1KG_Phenos.txt \
--make-bed \
--out 1kg_phase1_all.rsids.autosomal

# Obtain SNPs present in both files

awk '{print $2}' 1kg_phase1_all.rsids.autosomal.bim > 1kg_phase1_all.rsids_names.txt

# Extract 1KG SNPs from root

$plink \
--bfile $root.IBD_cleaned.rsids.autosomal \
--extract 1kg_phase1_all.rsids_names.txt \
--make-bed \
--out $root.IBD_cleaned.intersection

echo "*****************************************************************"
echo
echo " Test-merge 1KG and root data sets and output mismatching calls "
echo
echo "*****************************************************************"
echo

$plink \
--bfile $root.IBD_cleaned.intersection \
--bmerge 1kg_phase1_all.rsids.autosomal \
--merge-mode 6 \
--out $root.1KG.IBD_cleaned_failures

echo "*****************************************************************"
echo
echo " Check for SNPs with multiple positions "
echo
echo "*****************************************************************"
echo

# Log file would list problematic SNPs with quotes around them

fgrep \'rs $root.1KG.IBD_cleaned_failures.log |\

# Presumably column 7 in the log file contains the rsID flanked with ' and .

awk '{print $7}' |\
sed -e "s/'//g" -e "s/.//g" > $root.1KG.IBD_cleaned_failures.multiple.positions.txt

cat $root.1KG.IBD_cleaned_failures.missnp $root.1KG.IBD_cleaned_failures.multiple.positions.txt > $root.1KG.IBD_cleaned_failures.multiple.positions.missnp 

echo "*****************************************************************"
echo
echo " Exclude mismatched SNPs and those with multiple positions "
echo
echo "*****************************************************************"
echo

$plink \
--bfile $root.IBD_cleaned.intersection \
--exclude $root.1KG.IBD_cleaned_failures.multiple.positions.missnp \
--make-bed \
--out $root.IBD_cleaned.intersection_for_merge

echo "*****************************************************************"
echo
echo " Merge root and 1KG data sets "
echo
echo "*****************************************************************"
echo

$plink \
--bfile $root.IBD_cleaned.intersection_for_merge \
--bmerge 1kg_phase1_all.rsids.autosomal \
--out $root.1kg.pop_strat

echo "*****************************************************************"
echo
echo " Filter merged data for missing variants, rare variants and HWE fails "
echo
echo "*****************************************************************"
echo

$plink \
--bfile $root.1kg.pop_strat \
--geno 0.01 \
--maf 0.01 \
--hwe 0.0001 \
--make-bed \
--out $root.1kg.pop_strat.for_prune

echo "*****************************************************************"
echo
echo " LD pruning of merged data "
echo
echo "*****************************************************************"
echo

# Make prune.in and prune.out files

$plink \
--bfile $root.1kg.pop_strat.for_prune \
--indep-pairwise 1500 150 0.2 \
--out $root.1kg.pop_strat.prune

# Keep prune.in SNPs only

$plink \
--bfile $root.1kg.pop_strat.for_prune \
--extract $root.1kg.pop_strat.prune.prune.in \
--make-bed \
--out $root.1kg.LD_pop_strat

echo "*****************************************************************"
echo
echo " Run convertf to make EIGENSTRAT format file  "
echo
echo "*****************************************************************"
echo

convertf -p <(printf "genotypename: $root.1kg.LD_pop_strat.bed
             snpname: $root.1kg.LD_pop_strat.bim
             indivname: $root.1kg.LD_pop_strat.fam
             outputformat: EIGENSTRAT
             genotypeoutname: $root.1kg.LD_pop_strat.eigenstratgeno
             snpoutname: $root.1kg.LD_pop_strat.snp
             indivoutname: $root.1kg.LD_pop_strat.ind")

echo "*****************************************************************"
echo
echo " Make list of unique phenotype codes for projection "
echo
echo "*****************************************************************"
echo

awk '{print $3}' 1KG_Phenos.txt | sort | uniq > $root.1kg.LD_poplist.txt

echo "*****************************************************************"
echo
echo " Run Smartpca, projecting phenotype on 1KG samples only "
echo
echo "*****************************************************************"
echo

cd [location_of_scripts]

module load gnuplot

smartpca.perl \
-i $root.1kg.LD_pop_strat.eigenstratgeno \
-a $root.1kg.LD_pop_strat.snp \
-b $root.1kg.LD_pop_strat.ind \
-o $root.1kg.LD_pop_strat.pca \
-p $root.1kg.LD_pop_strat.plot \
-e $root.1kg.LD_pop_strat.eigenvalues \
-l $root.1kg.LD_pop_strat.log \
-w $root.1kg.LD_poplist.txt \
-m 0

module unload gnuplot

echo "*****************************************************************"
echo
echo " Modify merged .pca.evec file for R "
echo
echo "*****************************************************************"
echo

awk 'NR > 1 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12"CHANGE"}' $root.1kg.LD_pop_strat.pca.evec > $root.1kg.LD_pop_strat.pca.evec_RENAMED

# Change "Project" below to your project name

sed -i -e 's/CaseCHANGE/Project/g' -e 's/13CHANGE/LWK/g' -e 's/14CHANGE/MXL/g' -e 's/15CHANGE/PUR/g' -e 's/16CHANGE/TSI/g' -e 's/17CHANGE/YRI/g' -e 's/3CHANGE/ASW/g' -e 's/4CHANGE/CEU/g' -e 's/5CHANGE/CHB/g' -e 's/6CHANGE/CHS/g' -e 's/7CHANGE/CLM/g' -e 's/8CHANGE/FIN/g' -e 's/10CHANGE/GBR/g' -e 's/11CHANGE/IBS/g' -e 's/12CHANGE/JPT/g'  $root.1kg.LD_pop_strat.pca.evec_RENAMED

echo "*****************************************************************"
echo
echo " Plot PCs after merging with 1KG data "
echo
echo "*****************************************************************"
echo

R --file=PC_Plot_1KG.R --args $root

echo "*****************************************************************"
echo
echo " Heterozygosity Test "
echo
echo "*****************************************************************"
echo

# Test for unusual patterns of genome-wide heterogeneity in LD-pruned data

$plink \
--bfile $root.LD_IBD \
--ibc \
--out $root.het

# Exclude samples identified as outliers (user-defined number of SDs)
# I added the ability to specify the number of SDs

R --file=IdHets.R --args $root.het sdcut=3

$plink \
--bfile $root.LD_IBD \
--remove $root.het.LD_het_outliers_sample_exclude \
--make-bed \
--out $root.LD_het_cleaned

$plink \
--bfile $root.IBD_cleaned \
--remove $root.het.LD_het_outliers_sample_exclude \
--make-bed \
--out $root.het_cleaned

echo "*****************************************************************"
echo
echo " Zip old files "
echo
echo "*****************************************************************"
echo

cd [location_of_scripts]

mv *.pdf [location_of_data_to_keep]

cd [location_of_data]

tar -czvf [location_and_desired_filename].tar.gz *.IBD_outliers.* *.LD_IBD.* *.IBD_cleaned* *pop_strat.* *smartpca.log* *LD_poplist* *.het.*

rm *.IBD_outliers.* *.LD_IBD.* *.IBD_cleaned* *pop_strat.* *smartpca.log* *LD_poplist* *.het.*

mv [location_and_filename_of_new_archive] [location_of_data]

cd [location_of_data]

mv *het_cleaned* [location_of_data_to_keep]

exit
