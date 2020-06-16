#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -j y

echo "**********************************************************************"
echo 
echo " Defining variables and making config file "
echo
echo "**********************************************************************"
echo

printf "root=[project_filepath]
plink=[plink_binaries_filepath]" > Config.conf

source ./Config.conf

echo "**********************************************************************"
echo 
echo " Creating new PLINK files with only common (more than 1%) variants "
echo
echo "**********************************************************************"
echo

$plink \
  --bfile $root \
  --maf 0.01 \
  --make-bed \
  --out $root.common

echo "**********************************************************************"
echo
echo " Run Iterative Missingness script to remove SNPs, then samples, "
echo " at increasingly high cut-offs "
echo
echo "**********************************************************************"
echo

# N.B. This script uses --geno and --mind multiple times

sh ./Iterative_Missingness.sh 90 99 1 # % removal to start at, % removal to end at, step size (%)

source ./Config.conf

echo "**********************************************************************"
echo
echo " Review call rates to ensure all missing SNPs and individuals "
echo " have been dropped "
echo
echo "**********************************************************************"
echo

$plink \
  --bfile $root.filtered \
  --missing \
  --out $root.filtered_missing

echo "**********************************************************************"
echo
echo " Examine the lowest call rates for variants and individuals "
echo
echo "**********************************************************************"
echo

sort -k 5 -gr $root.filtered_missing.lmiss | head # Missing loci
sort -k 6 -gr $root.filtered_missing.imiss | head # Missing individuals

echo "**********************************************************************"
echo
echo " Assess SNPs for deviation from Hardy-Weinberg Equilibrium "
echo
echo "**********************************************************************"
echo

# Calculate p values for HWE test

$plink \
  --bfile $root.filtered \
  --hardy \
  --out $root.hw_p_values

# Remove deviant SNPs past a given threshold

$plink \
  --bfile $root.filtered \
  --hwe 0.00001 \
  --make-bed \
  --out $root.hw_dropped

echo "**********************************************************************"
echo
echo " Prune data for linkage disequilibrium "
echo
echo "**********************************************************************"
echo

# Look for pairs of loci that are uncorrelated
# indep-pairwise arguments are window size, shift size, r2 cutoff

$plink \
  --bfile $root.hw_dropped \
  --indep-pairwise 1500 150 0.2 \
  --out $root.LD_one

# Make a new bed file using only uncorrelated SNPs

$plink \
  --bfile $root.hw_dropped \
  --extract $root.LD_one.prune.in \
  --make-bed \
  --out $root.LD_two

# Exclude high-LD and non-autosomal regions from the pruned file

awk -f highLDregions4bim_b37.awk $root.LD_two.bim > highLDexcludes
awk '($1 < 1) || ($1 > 22) {print $2}' $root.LD_two.bim > autosomeexcludes
cat highLDexcludes autosomeexcludes > highLD_and_autosomal_excludes

$plink \
  --bfile $root.LD_two \
  --exclude highLD_and_autosomal_excludes \
  --make-bed \
  --out $root.LD_three

echo "*********************************************************************"
echo
echo " Pairwise identity-by-descent (IBD) check "
echo
echo "*********************************************************************"
echo

$plink \
  --bfile $root.LD_three \
  --genome \
  --make-bed \
  --out $root.IBD

# If applicable, remove samples with pi-hat (% IBD) above threshold 0.1875

awk '$10 >= 0.1875 {print $1, $2}' $root.IBD.genome > $root.IBD_outliers.txt

$plink \
  --bfile $root.IBD \
  --remove $root.IBD_outliers.txt \
  --make-bed \
  --out $root.no_close_relatives

# Calculate average IBD per individual using R, output outliers (>[user-defined] SDs above mean)

R --file=IndividualIBD.R --args $root 4

# Exclude outliers from both LD-stripped and all SNP binary files

$plink \
  --bfile $root.LD_three \
  --remove $root.IBD_INDIV_outliers.txt \
  --make-bed \
  --out $root.LD_IBD

$plink \
  --bfile $root.hw_dropped \
  --remove $root.IBD_INDIV_outliers.txt \
  --make-bed \
  --out $root.IBD_cleaned

echo "*********************************************************************"
echo
echo " Population stratification by PCA using Eigensoft "
echo
echo "*********************************************************************"
echo

# Convert files to EIGENSOFT format using CONVERTF

convertf -p <(printf "genotypename: "$root".LD_IBD.bed
snpname: "$root".LD_IBD.bim
indivname: "$root".LD_IBD.fam
outputformat: EIGENSTRAT
genotypeoutname: "$root".pop_strat.eigenstratgeno
snpoutname: "$root".pop_strat.snp
indivoutname: "$root".pop_strat.ind")

# Run SmartPCA

module load gnuplot

smartpca.perl \
-i $root.pop_strat.eigenstratgeno \
-a $root.pop_strat.snp \
-b $root.pop_strat.ind \
-o $root.pop_strat.pca \
-p $root.pop_strat.plot \
-e $root.pop_strat.eval \
-l $root.pop_strat_smartpca.log \
-m 0 \
-t 100 \
-k 100 \
-s 6

module unload gnuplot

# Remove leading tab from eigenvector file, and split ID into two columns, for R

sed -i -e 's/^[ \t]*//' -e 's/:/ /g' $root.pop_strat.pca.evec

# Plot principal components

# This is essentially a duplication of the gnuplot step from earlier

R --file=PlotPCs.R --args $root.pop_strat 1 2

# Re-run to assess which components to include as covariates in the final analysis

# This step has essentiall already been done but we're 
# giving the output files different names now

convertf -p <(printf "genotypename: $root.LD_IBD.bed
snpname: $root.LD_IBD.bim
indivname: $root.LD_IBD.fam
outputformat: EIGENSTRAT
genotypeoutname: $root.PCs_for_covariates.eigenstratgeno
snpoutname: $root.PCs_for_covariates.snp
indivoutname: $root.PCs_for_covariates.ind")

# Run SmartPCA

module load gnuplot

smartpca.perl \
-i $root.PCs_for_covariates.eigenstratgeno \
-a $root.PCs_for_covariates.snp \
-b $root.PCs_for_covariates.ind \
-o $root.PCs_for_covariates.pca \
-p $root.PCs_for_covariates.plot \
-e $root.PCs_for_covariates.eval \
-l $root.PCs_for_covariates_smartpca.log \
-m 0 \
-t 100 \
-k 100 \
-s 6 \

module unload gnuplot

echo "***********************************************************************"
echo
echo " Create .tar.gz archive 'qc_01.tar.gz' "
echo
echo "***********************************************************************"
echo

# Create a compressed archive of files no longer needed, and delete originals

cd [location_of_data]

tar -czvf [location_and_desired_filename].tar.gz *.IBD_INDIV.* *.IBD_INDIV_outliers.* *.no_close_relatives.* *.IBD_outliers.* *.IBD.* *.LD_three.* *.LD_two.* *.LD_one.* *.hw_dropped.* *.hw_p_values.* *.filtered_missing.* *.filtered.* *.common* 
cd [location_of_data]

rm *.IBD_INDIV.* *.IBD_INDIV_outliers.* *.no_close_relatives.* *.IBD_outliers.* *.IBD.* *.LD_three.* *.LD_two.* *.LD_one.* *.hw_dropped.* *.hw_p_values.* *.filtered_missing.* *.filtered.* *.common*

cd [location_of_scripts]

rm highLDexcludes autosomeexcludes highLD_and_autosomal_excludes

mv [location_and_filename_of_new_archive] [location_of_data]

mv *.pdf [location_of_data_to_keep]

exit
