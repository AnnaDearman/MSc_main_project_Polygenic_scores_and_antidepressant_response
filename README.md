### MSc_main_project_Polygenic_scores_and_antidepressant_response
A respository of files used in my six-month Master's project: "Polygenic scores for psychological traits and disorders and response and adverse drug reactions to antidepressants"

This was my six-month Master's project in 2020; my first big project in bioinformatics. At the time of writing, this is still a work-in-progress.
Overall, I am proud of how autonomously I was able to use and troubleshoot UNIX command-line tools, and bash and R scripts, given my limited experience and previous biology-only background! 
The command-line tools I used were all written by other people. Full credit goes to their respective authors.

A few notes to anyone who wants to follow along for their own project:

- This is not yet as user-friendly and generic as it could be. It is currently posted to demonstrate the kind of work I have done, not to help other researchers who are new to this type of analysis, though I would like it to be so in future, and would be willing to answer any questions people may have about these scripts! e-mail annalivings@yahoo.co.uk.
- There are some commands that compress and delete files. I had to compress my files to save space - others may like to omit these steps.
- I have redacted some file paths - look through the scripts and fill in your own file paths.
- I advise anyone following this to thoroughly research the tools I've used. Go to https://www.cog-genomics.org/plink/2.0, https://github.com/JoniColeman/gwas_scripts, https://sites.google.com/site/mikeweale/software/eigensoftplus, https://www.well.ox.ac.uk/~wrayner/tools/, https://github.com/zhanxw/checkVCF, https://www.well.ox.ac.uk/~wrayner/tools/Post-Imputation.html, https://doc-openbio.readthedocs.io/projects/annovar/en/latest/, https://www.prsice.info/, etc.
- There may be some files and programs that I use that are not included in my repository - if so, apologies: just try searching for the script name online!
- My GWAS data did not include phenotypes or controls - your data, and therefore scripting needs, may be very different! Try https://github.com/JoniColeman/gwas_scripts
- My GWAS data had already been sex-check-cleaned by somebody else.
- Check whether the "module load" and "module unload" statements are applicable to you.

*STEP_01_GWAS_QC_pre-imputation

The order in which these scripts should be run is:

qc_01.sh, qc_02.sh, qc_03.sh, pre_imp_01.sh, pre_imp_02.sh, pre_imp_03.sh. The other scripts are invoked by running these scripts.

The commands in qc_01.sh and qc_02.sh are slightly adapted, but are basically copied, from the respository https://github.com/JoniColeman/gwas_scripts (many thanks!). Some steps were removed because my data did not contain controls and phenotypes. 
In addition to these, there are other scripts and files copied or adapted from https://github.com/JoniColeman/gwas_scripts.
qc_03.sh, which runs EigensoftPlus, was added retrospectively to remove population stratification outliers, because I had been unable to do so when using Eigensoft earlier in the process. This seems a little messy.

The commands in pre_imp_01.sh and pre_imp_02.sh are based on the advice given on the Michigan Imputation Server website: https://imputationserver.readthedocs.io/en/latest/prepare-your-data/
I have included HRC-1000G-check-bim-NoReadKey.pl simply because it contains a tiny edit by me - my first ever taste of perl scripting! Search for "Added by ARD" and you will find it. Full credit for the remaining 99.999% of the script goes to https://www.well.ox.ac.uk/~wrayner/tools/ !
checkVCF.py can be copied from from https://github.com/zhanxw/checkVCF - I have not included it here.

I have not yet produced a proper allele frequency plot. The plot produced in pre_imp_03.sh is not sufficient.

*STEP_02_GWAS_QC_post-imputation

The order in which these scripts should be run is:

post_imp_01.sh, post_imp_02.sh, post_imp_03.sh, post_imp_04.sh
These scripts unpack your imputed data, run diagnostic plots with IC, run some PLINK commands, and replace the SNP IDs using Annovar. I adapted one of my supervisor's previous scripts to run the PLINK commands, but I figured out how to use the other tools myself.

*STEP_03_Power_calculations

I read up on many different psychiatric traits and their GWASs to decide which to generate polygenic scores for, to use as explanatory variables in my analysis of antidepressant response. I also performed a power calculation for each of them, using the avengeme package, polygenescore() function.
In order to understand the effects of different values for the different inputs into the power calculation for my polygenic scores, I made a few scripts that produced plots: Binary_training_trait_power_calc_with_plots.R, Continuous_training_trait_power_calc_with_plots.R and Power_calc_as_function_of_n_ind_training_trait_and_prevalence_with_plots.R.
When all the comments and vectors of numbers got a bit too much to look at, I used Quick_power_calc.R to quickly and easily produce a single power calculation for a set of inputs.

*STEP_04_Polygenic_scores

One simple script to run PRSice to generate polygenic scores for four psychiatric disorders and traits, for my participants.
N.B. I did not have to convert the human genome annotation to a different version / "build" in order to run PRSice on this occasion. If so, you need to make non-binary versions of the genotype data in PLINK, convert it using LiftOver, then make binary versions again. 

*STEP_05_Data_analysis

I am still working on this. It will be a very long R script eventually!

Areas for growth:
- The commands I ran are separated into separate scripts somewhat arbitrarily.
- Ideally these scripts could be made more user-friendly by coding them so that a user could submit a few key pieces of information as arguments and run the scripts in a more hands-off way. However, having to go through the code with a fine tooth comb would give the user greater expertise!
