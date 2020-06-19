#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=20G   # Request 20GB RAM
#$ -j y

echo "***********************************************************************"
echo "                 Polygenic scoring - anorexia   "
echo "***********************************************************************"

cd ~[file path]

# Removes header lines - always check summary stats files first

sed '1,70d' pgcAN2.2019-07.vcf.tsv > pgcAN2.2019-07_1.vcf.tsv

PRSice_linux \
 --A1 REF \
 --A2 ALT \
 --pvalue PVAL \
 --stat BETA \
 --snp ID \
 --base ~[file path]/pgcAN2.2019-07_1.vcf.tsv \
 --fastscore \
 --info-base IMPINFO,0.9 \
 --all-score \
 --no-regress \
 --print-snp \
 --bar-levels 0.0001,0.01,0.05,0.1,0.5,1 \
 --target ~[file path]/[project]_1kg_imp_QCd_popstratout_removed \
 --out ~[file path]/[project]_anorexia_PGS

echo "***********************************************************************"
echo "                 Polygenic scoring - reaction time   "
echo "***********************************************************************"

PRSice_linux \
 --A1 Effect_allele \
 --A2 Other_allele \
 --pvalue P \
 --stat Beta \
 --snp MarkerName \
 --base ~[file path]/Davies2018_UKB_RT_summary_results_29052018.txt \
 --fastscore \
 --all-score \
 --no-regress \
 --print-snp \
 --bar-levels 0.0001,0.01,0.05,0.1,0.5,1 \
 --target ~[file path]/[project]_1kg_imp_QCd_popstratout_removed \
 --out ~[file path]/[project]_ReactTime_PGS

echo "***********************************************************************"
echo "                 Polygenic scoring - neuroticism   "
echo "***********************************************************************"

PRSice_linux \
 --A1 effect_allele \
 --A2 other_allele \
 --pvalue p_value \
 --stat beta \
 --snp variant_id \
 --base ~[file path]/29255261-GCST005232-EFO_0007660-build37.f.tsv \
 --fastscore \
 --all-score \
 --no-regress \
 --print-snp \
 --bar-levels 0.0001,0.01,0.05,0.1,0.5,1 \
 --target ~[file path]/[project]_1kg_imp_QCd_popstratout_removed \
 --out ~[file path]/[project]_neuroticism_PGS

echo "***********************************************************************"
echo "                 Polygenic scoring - schizophrenia    "
echo "***********************************************************************"

PRSice_linux \
 --A1 A1 \
 --A2 A2 \
 --pvalue P \
 --stat OR \
 --snp SNP\
 --base ~[file path]/baseSCZ.assoc \
 --fastscore \
 --all-score \
 --no-regress \
 --print-snp \
 --bar-levels 0.0001,0.01,0.05,0.1,0.5,1 \
 --target ~[file path]/[project]_1kg_imp_QCd_popstratout_removed \
 --out ~[file path]/[project]_schizophrenia_PGS

exit