##### 01) INSTALL PACKAGES #####

if ( !require( haven ) ) install.packages( "haven" ) # to read in Stata .dta files
if ( !require( dplyr ) ) install.packages( "dplyr" ) # for manipulating data etc
if ( !require( tidyr ) ) install.packages( "tidyr" ) # for tidy data
if ( !require( ggplot2 ) ) install.packages( "ggplot2" ) # for plotting
if ( !require( gridExtra ) ) install.packages( "gridExtra" ) # for PDFs
if ( !require( ggcorrplot ) ) install.packages( "gridExtra" ) # For correlation plots
if ( !require( lme4 ) ) install.packages( "lme4" ) # linear mixed effects models
if ( !require( lmerTest ) ) install.packages( "lmerTest" ) # linear mixed effects models
if ( !require( broom.mixed ) ) install.packages( "broom.mixed" )
if( !require( predictmeans)) install.packages( "predictmeans" ) # for diagnostic plots
if( !require( LMERConvenienceFunctions )) install.packages( "LMERConvenienceFunctions" ) # for diagnostic plots
if ( !require( MuMIn ) ) install.packages( "MuMIn" ) # to assess models using explained variance
if ( !require( gplots )) install.packages( "gplots" ) # for heatmaps
if ( !require( survival )) install.packages( "survival" ) # for survival analyses
if ( !require( survminer )) install.packages( "survminer" ) # for survival analyses
if ( !require( reshape2 )) install.packages( "reshape2" ) # for survival analyses
if( !require( survMisc)) install.packages( "survMisc" ) # fo get R2 from survival analysis

##### 02) IMPORT LIBRARIES #####

library( haven )
library( dplyr )
library( tidyr )
library( ggplot2 )
library( gridExtra )
library( ggcorrplot )
library( lme4 )
library( lmerTest )
library( broom.mixed )
library( predictmeans )
library( LMERConvenienceFunctions )
library( MuMIn )
library( gplots )
library( survival )
library( survminer )
library( reshape2 )
library( survMisc )

##### 03) SET WORKING DIRECTORY #####

setwd( [redacted] )

##### 04) IMPORT DATA FILES #####

# Import Polygenic risk score files

anorexia <- read.table( "gendep_anorexia_PGS.all.score", 
                        skip = 1,
                        col.names = c( "FID", "bloodsampleid", "A_0_0001", "A_0_01", "A_0_05", "A_0_1", "A_0_5", "A_1" ) )

cognitive <- read.table( "gendep_ReactTime_PGS.all.score", 
                         skip = 1,
                         col.names = c( "FID", "bloodsampleid", "C_0_0001", "C_0_01", "C_0_05", "C_0_1", "C_0_5", "C_1" ) )

neuroticism <- read.table( "gendep_neuroticism_PGS.all.score", 
                           skip = 1,
                           col.names = c( "FID", "bloodsampleid", "N_0_0001", "N_0_01", "N_0_05", "N_0_1", "N_0_5", "N_1" ) )

schizophrenia <- read.table( "gendep_schizophrenia_PGS.all.score", 
                             skip = 1,
                             col.names = c( "FID", "bloodsampleid", "S_0_0001", "S_0_01", "S_0_05", "S_0_1", "S_0_5", "S_1" ) )

height <- read.table( "gendep_height_PGS.all.score",
                      skip = 1,
                      col.names = c( "FID", "bloodsampleid", "H_0_0001", "H_0_01", "H_0_05", "H_0_1", "H_0_5", "H_1" ) )

# Import side effect data (ASEC) and main GENDEP data

ASEC <- read_dta( [redacted] )
GENDEP <- read_dta( [redacted] )
GENDEP_raw <- read_dta( [redacted] )

# Import principal component score covariates
# N.B. I manually removed the first line that says "eigvals"
# prior to import

PCs <- read.table( "[redacted].pca.evec" )

##### 05) MAKE IMPORTANT CHANGES TO DATA FILES #####

# Make certain variables factors, not numerical

GENDEP$sex <- as.factor(GENDEP$sex)
GENDEP$drug <- as.factor(GENDEP$drug)
GENDEP$centreid <- as.factor(GENDEP$centreid)
GENDEP$subjectid <- as.factor(GENDEP$subjectid)
ASEC$subjectid <- as.factor(ASEC$subjectid)
GENDEP$random<- as.factor(GENDEP$random)

# Make a new "dropout" column that relates to "dropping out"
# of the FIRST drug, even if they didn't "drop out" of the study

GENDEP$dropout2 <- GENDEP$completer
GENDEP$dropout2[GENDEP$dropout2==0] <- "d"
GENDEP$dropout2[GENDEP$dropout2==1] <- 0
GENDEP$dropout2[GENDEP$dropout2=="d"] <- 1
GENDEP$dropout2 <- as.numeric(GENDEP$dropout2)

# Tidy up and label principal components data

PCs <- PCs[,2:102] # Ignore first column; it's effectively a duplicate of the second
rownames( PCs ) <- PCs[,1] # Use IDs as row names
colnames( PCs ) <- c( "bloodsampleid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
                      "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", 
                      "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", 
                      "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36", 
                      "PC37", "PC38", "PC39", "PC40", "PC41", "PC42", "PC43", "PC44", "PC45", "PC46", 
                      "PC47", "PC48", "PC49", "PC50", "PC51", "PC52", "PC53", "PC54", "PC55", "PC56", 
                      "PC57", "PC58", "PC59", "PC60", "PC61", "PC62", "PC63", "PC64", "PC65", "PC66", 
                      "PC67", "PC68", "PC69", "PC70", "PC71", "P72C", "PC73", "PC74", "PC75", "PC76", 
                      "PC77", "PC78", "PC79", "PC80", "PC81", "PC82", "PC83", "PC84", "PC85", "PC86", 
                      "PC87", "PC88", "PC89", "PC90", "PC91", "PC92", "PC93", "PC94", "PC95", "PC96", 
                      "PC97", "PC98", "PC99", "PC100" )

# Make scaled principal component columns for PCs 1, 2 and 3

PCs <- PCs %>%
  mutate( zPC1 = scale(PC1))
PCs <- PCs %>%
  mutate( zPC2 = scale(PC2))
PCs <- PCs %>%
  mutate( zPC3 = scale(PC3))

# Merge the GENDEP and ASEC files

GENDEP <- merge( GENDEP, ASEC, by=c("subjectid", "week"),
                 all.x=TRUE, all.y=TRUE )

# Add "totasec" (total side effect count) column

GENDEP <- GENDEP %>%
  mutate( totasec = asec1wk + asec2wk + asec3wk +
            asec4wk + asec5wk + asec6wk + asec7wk +
            asec8wk + asec9wk + asec10wk + asec11wk + 
            asec12wk + asec13wk + asec14wk + asec15wk + 
            asec16wk + asec17wk + asec18wk + asec19wk + asec20wk + asec21wk)

# Add "ztotasec", scaled total ASEC side effect count

GENDEP <- GENDEP %>%
  mutate( ztotasec = scale(totasec) )

# Get total baseline ASEC per subject

bltotasec <- data.frame(
  cbind(as.character(subset(GENDEP,week==0)$subjectid),
        subset(GENDEP,week==0)$totasec)
)
colnames(bltotasec)<-c("subjectid","bltotasec")
bltotasec$subjectid <- as.factor(bltotasec$subjectid)

# Add baseline total ASEC as a column by merging

GENDEP <- merge( GENDEP, bltotasec, by="subjectid",
                 all.x=TRUE, all.y=TRUE )
GENDEP$bltotasec <- as.numeric(GENDEP$bltotasec)

# Remove observations where week <1 or week >12 

GENDEP <- subset( GENDEP, week>0 & week<13 )

# Keep only the participants with PRS:

## Get list of participants in PRS file

participants <- anorexia$bloodsampleid

## Subset GENDEP file based on these IDs

GENDEP <- subset( GENDEP, bloodsampleid %in% participants )

# Make selection vectors for each drug, and a vector 
# called "Drug" that can be used later to colour plots

Esc_boo <- GENDEP$drug==0
Nor_boo <- GENDEP$drug==1
Drug <- NULL
Drug[Esc_boo] <- "Escitalopram"
Drug[Nor_boo] <- "Nortriptyline"

# For each trait whose PRS is being considered,
# center and scale the participants' PRSs (i.e. 
# create Z scores), for each p-value threshold
# and add these as extra columns 

Z_thresh_A <- c( "AZ_0_0001", "AZ_0_01", "AZ_0_05", "AZ_0_1", "AZ_0_5", "AZ_1" )
Z_thresh_C <- c( "CZ_0_0001", "CZ_0_01", "CZ_0_05", "CZ_0_1", "CZ_0_5", "CZ_1" )
Z_thresh_N <- c( "NZ_0_0001", "NZ_0_01", "NZ_0_05", "NZ_0_1", "NZ_0_5", "NZ_1" )
Z_thresh_S <- c( "SZ_0_0001", "SZ_0_01", "SZ_0_05", "SZ_0_1", "SZ_0_5", "SZ_1" )
Z_thresh_H <- c( "HZ_0_0001", "HZ_0_01", "HZ_0_05", "HZ_0_1", "HZ_0_5", "HZ_1" )

names <- colnames( anorexia ) # Grab the column names as they are
for ( i in 3:8 ){ # For every column with PRSs in
  anorexia <- cbind( anorexia, scale( anorexia[,i] ) ) # add the Z scores to the data frame
}
colnames( anorexia ) <- c( names, Z_thresh_A ) # Add the old and new column names

names <- colnames( cognitive )
for ( i in 3:8 ){
  cognitive <- cbind( cognitive, scale( cognitive[,i] ) )
}
colnames( cognitive ) <- c( names, Z_thresh_C )

names <- colnames( neuroticism )
for ( i in 3:8 ){
  neuroticism <- cbind( neuroticism, scale( neuroticism[,i] ) )
}
colnames( neuroticism ) <- c( names, Z_thresh_N )

names <- colnames( schizophrenia )
for ( i in 3:8 ){
  schizophrenia <- cbind( schizophrenia, scale( schizophrenia[,i] ) )
}
colnames( schizophrenia ) <- c( names, Z_thresh_S )

names <- colnames( height )
for ( i in 3:8 ){
  height <- cbind( height, scale( height[,i] ) )
}
colnames( height ) <- c( names, Z_thresh_H )

# Bind the PRS results with the GENDEP file

GENDEP <- merge( GENDEP, anorexia, by="bloodsampleid" )
GENDEP$FID <- NULL
GENDEP <- merge( GENDEP, cognitive, by="bloodsampleid" )
GENDEP$FID <- NULL
GENDEP <- merge( GENDEP, neuroticism, by="bloodsampleid" )
GENDEP$FID <- NULL
GENDEP <- merge( GENDEP, schizophrenia, by="bloodsampleid" )
GENDEP$FID <- NULL
GENDEP <- merge( GENDEP, height, by="bloodsampleid" )
GENDEP$FID <- NULL

# Add the PC covariates to the GENDEP file

GENDEP <- merge( GENDEP, PCs, by="bloodsampleid" )

# Sort

GENDEP <- GENDEP[order(GENDEP$subjectid, as.numeric(GENDEP$week)),]

# Make GENDEP sub-tables for each drug

GENDEP_escit <- subset( GENDEP, drug==0 ) 
GENDEP_nortrip <- subset( GENDEP, drug==1 )

# Make GENDEP sub-tables for randomised and non-randomised

GEND_nonrand <- subset( GENDEP, random==0 )
GEND_rand <- subset( GENDEP, random==1 )

# Break this down further, by drug

GEND_rand_e <- subset(GEND_rand,drug==0)
GEND_rand_n <- subset(GEND_rand,drug==1)
GEND_nonrand_e <- subset(GEND_nonrand,drug==0)
GEND_nonrand_n <- subset(GEND_nonrand,drug==1)

# Make dataset with one obs per individual.
# To help with survival analysis later, choose 
# the "last" observation for each individual (end week)

GEND_end <- subset(GENDEP,endweek==week)

# Break this down further for descriptive statistics

GE_esc <- subset(GEND_end,drug==0)
GE_nor <- subset(GEND_end,drug==1)
GE_M <- subset(GEND_end,sex==1)
GE_F <- subset(GEND_end,sex==2)
GE_nran <- subset(GEND_end, random==0)
GE_ran <- subset(GEND_end, random==1)
GE_drop <- subset(GEND_end, dropout==1)
GE_nondrop <- subset(GEND_end, dropout==0)
GE_nonswitch <- subset(GEND_end, switcher==0)
GE_ns_esc <- subset(GE_nonswitch, drug==0)
GE_ns_nor <- subset(GE_nonswitch, drug==1)
GE_ran_ns <- subset(GE_ran, switcher==0)
GE_ran_ns_e <- subset(GE_ran_ns, drug==0)
GE_ran_ns_n <- subset(GE_ran_ns, drug==1)
GE_ran_e <- subset(GE_ran, drug==0)
GE_ran_n <- subset(GE_ran, drug==1)

##### 06) MAKE IMPORTANT PRS NAMING VECTORS FOR LOOPING AND NAMING #####

PRSs <- c( "AZ_0_0001", "AZ_0_01", "AZ_0_05", "AZ_0_1", "AZ_0_5", "AZ_1",
           "CZ_0_0001", "CZ_0_01", "CZ_0_05", "CZ_0_1", "CZ_0_5", "CZ_1",
           "NZ_0_0001", "NZ_0_01", "NZ_0_05", "NZ_0_1", "NZ_0_5", "NZ_1",
           "SZ_0_0001", "SZ_0_01", "SZ_0_05", "SZ_0_1", "SZ_0_5", "SZ_1",
           "HZ_0_0001", "HZ_0_01", "HZ_0_05", "HZ_0_1", "HZ_0_5", "HZ_1" )

# Vector for PRSs, with "No_PRS" at the end

PRSs_31 <- c( "AZ_0_0001", "AZ_0_01", "AZ_0_05", "AZ_0_1", "AZ_0_5", "AZ_1",
              "CZ_0_0001", "CZ_0_01", "CZ_0_05", "CZ_0_1", "CZ_0_5", "CZ_1",
              "NZ_0_0001", "NZ_0_01", "NZ_0_05", "NZ_0_1", "NZ_0_5", "NZ_1",
              "SZ_0_0001", "SZ_0_01", "SZ_0_05", "SZ_0_1", "SZ_0_5", "SZ_1",
              "HZ_0_0001", "HZ_0_01", "HZ_0_05", "HZ_0_1", "HZ_0_5", "HZ_1",
              "No_PRS")

##### 07) EXPLORE DATA #####

##### 07a) Visualise descriptive statistics #####

pdf( "Descriptive plots.pdf")

# Sex breakdown by drug

esex <- matrix(as.vector(table(GE_esc$sex)),nrow=1,dimnames=
                 list(c("Escitalopram"),c("Male","Female")))
nsex <- matrix(as.vector(table(GE_nor$sex)),nrow=1,dimnames=
                 list(c("Nortriptyline"),c("Male","Female")))
sex_drug <- rbind(esex,nsex)
barplot(sex_drug,names=c("Males","Females"),col=c("lightpink1","yellowgreen"),
        border=NA,main="Sex breakdown by first drug",beside=TRUE,
        args.legend = list(title = "Drug", x = "topright"),
        xlim = c(0, 8), ylim = c(0, 300),
        legend=TRUE)
text(x=c(1.5,2.5,4.5,5.5),y=75,labels=sex_drug)

# Drug breakdown by sex

mdrug <- matrix(as.vector(table(GE_M$drug)),nrow=1,dimnames=
                  list(c("Males"),c("Escitalopram","Nortriptyline")))
fdrug <- matrix(as.vector(table(GE_F$drug)),nrow=1,dimnames=
                  list(c("Females"),c("Escitalopram","Nortriptyline")))
drug_sex <- rbind(mdrug,fdrug)
barplot(drug_sex,names=c("Escitalopram","Nortriptyline"),col=c("purple","orange"),
        border=NA,main="First drug breakdown by sex",beside=TRUE,
        args.legend = list(title = "Sex", x = "topright"),
        xlim = c(0, 8), ylim = c(0, 300),
        legend=TRUE)
text(x=c(1.5,2.5,4.5,5.5),y=75,labels=drug_sex)

# Dropout rate by drug

edrop <- matrix(as.vector(table(GE_esc$dropout)),nrow=1,dimnames=
                  list(c("Escitalopram"),c("No dropout","Dropout")))
ndrop <- matrix(as.vector(table(GE_nor$dropout)),nrow=1,dimnames=
                  list(c("Nortriptyline"),c("No dropout","Dropout")))
drop_drug <- rbind(edrop,ndrop)
barplot(drop_drug,names=c("No dropout","Dropout"),col=c("lightpink1","yellowgreen"),
        border=NA,main="Dropout rate by first drug",beside=TRUE,
        args.legend = list(title = "Drug", x = "topright"),
        xlim = c(0, 8), ylim = c(0, 300),
        legend=TRUE)
text(x=c(1.5,2.5,4.5,5.5),y=50,labels=drop_drug)

# Dropout rate by sex

mdrop <- matrix(as.vector(table(GE_M$dropout)),nrow=1,dimnames=
                  list(c("Males"),c("No dropout","Dropout")))
fdrop <- matrix(as.vector(table(GE_F$dropout)),nrow=1,dimnames=
                  list(c("Females"),c("No dropout","Dropout")))
drop_sex <- rbind(mdrop,fdrop)
barplot(drop_sex,names=c("No dropout","Dropout"),col=c("purple","orange"),
        border=NA,main="Dropout rate by sex",beside=TRUE,
        args.legend = list(title = "Sex", x = "topright"),
        xlim = c(0, 8), ylim = c(0, 350),
        legend=TRUE)
text(x=c(1.5,2.5,4.5,5.5),y=30,labels=drop_sex)

# End week by drug

eend <- matrix(as.vector(table(GE_esc$endweek)),nrow=1,dimnames=
                 list(c("Escitalopram"),c(1:12)))
nend <- matrix(as.vector(table(GE_nor$endweek)),nrow=1,dimnames=
                 list(c("Nortriptyline"),c(1:12)))
end_drug <- rbind(eend,nend)
barplot(end_drug,names=c(1:12),col=c("lightpink1","yellowgreen"),
        border=NA,main="End week (of first drug) by first drug",beside=TRUE,
        args.legend = list(title = "Drug", x = "topright"),
        xlim = c(0, 40), ylim = c(0, 350),
        legend=TRUE)

# End week by sex

mend <- matrix(as.vector(table(GE_M$endweek)),nrow=1,dimnames=
                 list(c("Males"),c(1:12)))
fend <- matrix(as.vector(table(GE_F$endweek)),nrow=1,dimnames=
                 list(c("Females"),c(1:12)))
end_sex <- rbind(mend,fend)
barplot(end_sex,names=c(1:12),col=c("purple","orange"),
        border=NA,main="End week (of first drug) by sex",beside=TRUE,
        args.legend = list(title = "Sex", x = "topright"),
        xlim = c(0, 40), ylim = c(0, 400),
        legend=TRUE)

# End week by randomized (yes/no)

nonrand_end <- matrix(as.vector(table(GE_nran$endweek)),nrow=1,dimnames=
                        list(c("Non_random"),c(1:12)))
rand_end <- matrix(as.vector(table(GE_ran$endweek)),nrow=1,dimnames=
                     list(c("Random"),c(1:12)))
end_rand <- rbind(nonrand_end,rand_end)
barplot(end_rand,names=c(1:12),col=c("red","black"),
        border=NA,main="End week (of first drug) by randomization",beside=TRUE,
        args.legend = list(title = "Randomized?", x = "topright"),
        xlim = c(0, 40), ylim = c(0, 400),
        legend=TRUE)

# Randomized yes/no by sex

mrand <- matrix(as.vector(table(GE_M$random)),nrow=1,dimnames=
                  list(c("Males"),c("Non-random","Random")))
frand <- matrix(as.vector(table(GE_F$random)),nrow=1,dimnames=
                  list(c("Females"),c("Non-random","Random")))
rand_sex <- rbind(mrand,frand)
barplot(rand_sex,names=c("Non-random","Random"),col=c("purple","orange"),
        border=NA,main="Randomized (yes/no) by sex",beside=TRUE,
        args.legend = list(title = "Sex", x = "topright"),
        ylim = c(0, 400),
        legend=TRUE)
text(x=c(1.5,2.5,4.5,5.5),y=30,labels=rand_sex)

# Randomized yes/no by drug

erand <- matrix(as.vector(table(GE_esc$random)),nrow=1,dimnames=
                  list(c("Escitalopram"),c("Non-random","Random")))
nrand <- matrix(as.vector(table(GE_nor$random)),nrow=1,dimnames=
                  list(c("Nortriptyline"),c("Non-random","Random")))
rand_drug <- rbind(erand,nrand)
barplot(rand_drug,names=c("Non-random","Random"),col=c("lightpink1","yellowgreen"),
        border=NA,main="Randomized (yes/no) by first drug",beside=TRUE,
        args.legend = list(title = "Drug", x = "topright"),
        ylim = c(0, 300),
        legend=TRUE)
text(x=c(1.5,2.5,4.5,5.5),y=30,labels=rand_drug)

# End week by drug: randomized "yes"

GE_esc_rand <- subset(GE_esc, random==1)
GE_nor_rand <- subset(GE_nor, random==1)

esc_end_rand <- matrix(as.vector(table(GE_esc_rand$endweek)),nrow=1,dimnames=
                         list(c("Escitalopram"),c(1,3,4,5,6,7,8,9,10,11,12)))
esc_end_rand <- matrix(as.vector(c(esc_end_rand[1],0,esc_end_rand[2:11])),nrow=1,
                       dimnames=list(c("Escitalopram"),c(1:12)))
nor_end_rand <- matrix(as.vector(table(GE_nor_rand$endweek)),nrow=1,dimnames=
                         list(c("Nortriptyline"),c(1:12)))
end_drug_rand <- rbind(esc_end_rand,nor_end_rand)
barplot(end_drug_rand,names=c(1:12),col=c("lightpink1","yellowgreen"),
        border=NA,main="End week by first drug, randomized only",beside=TRUE,
        args.legend = list(title = "Drug", x = "topright"),
        xlim = c(0, 40), ylim = c(0, 250),
        legend=TRUE)

# End week by drug: randomized "no"

GE_esc_nonrand <- subset(GE_esc, random==0)
GE_nor_nonrand <- subset(GE_nor, random==0)

esc_end_nrand <- matrix(as.vector(table(GE_esc_nonrand$endweek)),nrow=1,dimnames=
                          list(c("Escitalopram"),c(1:12)))
nor_end_nrand <- matrix(as.vector(table(GE_nor_nonrand$endweek)),nrow=1,dimnames=
                          list(c("Nortriptyline"),c(1:12)))
end_drug_nrand <- rbind(esc_end_nrand,nor_end_nrand)
barplot(end_drug_nrand,names=c(1:12),col=c("lightpink1","yellowgreen"),
        border=NA,main="End week by first drug, non-randomized only",beside=TRUE,
        args.legend = list(title = "Drug", x = "topright"),
        xlim = c(0, 40), ylim = c(0, 250),
        legend=TRUE)

# Looking at PRSs and drug allocation (random v non-random)

PRS_means_rand <- GEND_end %>%
  group_by(random) %>%
  summarize(
    AZ_0_0001_mean = mean(AZ_0_0001),
    AZ_0_01_mean = mean(AZ_0_01),
    AZ_0_05_mean = mean(AZ_0_05),
    AZ_0_1_mean = mean(AZ_0_1),
    AZ_0_5_mean = mean(AZ_0_5),
    AZ_1_mean = mean(AZ_1),
    CZ_0_0001_mean = mean(CZ_0_0001),
    CZ_0_01_mean = mean(CZ_0_01),
    CZ_0_05_mean = mean(CZ_0_05),
    CZ_0_1_mean = mean(CZ_0_1),
    CZ_0_5_mean = mean(CZ_0_5),
    CZ_1_mean = mean(CZ_1),
    NZ_0_0001_mean = mean(NZ_0_0001),
    NZ_0_01_mean = mean(NZ_0_01),
    NZ_0_05_mean = mean(NZ_0_05),
    NZ_0_1_mean = mean(NZ_0_1),
    NZ_0_5_mean = mean(NZ_0_5),
    NZ_1_mean = mean(NZ_1),
    SZ_0_0001_mean = mean(SZ_0_0001),
    SZ_0_01_mean = mean(SZ_0_01),
    SZ_0_05_mean = mean(SZ_0_05),
    SZ_0_1_mean = mean(SZ_0_1),
    SZ_0_5_mean = mean(SZ_0_5),
    SZ_1_mean = mean(SZ_1),
    HZ_0_0001_mean = mean(HZ_0_0001),
    HZ_0_01_mean = mean(HZ_0_01),
    HZ_0_05_mean = mean(HZ_0_05),
    HZ_0_1_mean = mean(HZ_0_1),
    HZ_0_5_mean = mean(HZ_0_5),
    HZ_1_mean = mean(HZ_1),
    
    AZ_0_0001_sd = sd(AZ_0_0001),
    AZ_0_01_sd = sd(AZ_0_01),
    AZ_0_05_sd = sd(AZ_0_05),
    AZ_0_1_sd = sd(AZ_0_1),
    AZ_0_5_sd = sd(AZ_0_5),
    AZ_1_sd = sd(AZ_1),
    CZ_0_0001_sd = sd(CZ_0_0001),
    CZ_0_01_sd = sd(CZ_0_01),
    CZ_0_05_sd = sd(CZ_0_05),
    CZ_0_1_sd = sd(CZ_0_1),
    CZ_0_5_sd = sd(CZ_0_5),
    CZ_1_sd = sd(CZ_1),
    NZ_0_0001_sd = sd(NZ_0_0001),
    NZ_0_01_sd = sd(NZ_0_01),
    NZ_0_05_sd = sd(NZ_0_05),
    NZ_0_1_sd = sd(NZ_0_1),
    NZ_0_5_sd = sd(NZ_0_5),
    NZ_1_sd = sd(NZ_1),
    SZ_0_0001_sd = sd(SZ_0_0001),
    SZ_0_01_sd = sd(SZ_0_01),
    SZ_0_05_sd = sd(SZ_0_05),
    SZ_0_1_sd = sd(SZ_0_1),
    SZ_0_5_sd = sd(SZ_0_5),
    SZ_1_sd = sd(SZ_1),
    HZ_0_0001_sd = sd(HZ_0_0001),
    HZ_0_01_sd = sd(HZ_0_01),
    HZ_0_05_sd = sd(HZ_0_05),
    HZ_0_1_sd = sd(HZ_0_1),
    HZ_0_5_sd = sd(HZ_0_5),
    HZ_1_sd = sd(HZ_1)
  )

PRS_means_rand <- PRS_means_rand[,2:61]
PRS_means_rand <- melt(PRS_means_rand)
PRS_means_rand <- cbind(PRS_means_rand[1:60,],
                        PRS_means_rand[61:120,2])
PRS_means_rand <- cbind(PRS_means_rand,rep(c("Non-random","Random"),30))
colnames(PRS_means_rand) <- c("PRS", "Group_mean", "SD", "Random")

ggplot(PRS_means_rand, 
       aes(x = Group_mean, y = reorder(PRS, desc(PRS)), 
           col = Random)) +
  geom_point() + 
  geom_errorbar(aes(xmin=Group_mean-SD, 
                    xmax=Group_mean+SD), width=.2,
                position=position_dodge(.9))+
  xlab('Polygenic score') +
  ylab('Group mean')+
  theme(axis.title=element_text(size=7))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="Were they randomized based on their polygenic score?")+
  scale_colour_manual(labels=c("Non-random", "Random"),
                      values=c("red", "black"))

# The furthest apart dots

ggplot(GEND_end, aes(x=HZ_0_1,color=random))+
  geom_density()+
  scale_colour_manual(labels=c("Non-random", "Random"),
                      values=c("red", "black"))+
  labs(title="Were they randomized based on their HZ_0_1 polygenic score?")

# Check whether polygenic score means are significantly
# different for randomized and non-randomized people

t.test(GE_nran$HZ_0_1,GE_ran$HZ_0_1)

# Looking at PRSs and drug allocation (escitalopram v nortriptyline)
# non-random only

PRS_means_drug <- GEND_nonrand %>%
  group_by(drug) %>%
  summarize(
    AZ_0_0001_mean = mean(AZ_0_0001),
    AZ_0_01_mean = mean(AZ_0_01),
    AZ_0_05_mean = mean(AZ_0_05),
    AZ_0_1_mean = mean(AZ_0_1),
    AZ_0_5_mean = mean(AZ_0_5),
    AZ_1_mean = mean(AZ_1),
    CZ_0_0001_mean = mean(CZ_0_0001),
    CZ_0_01_mean = mean(CZ_0_01),
    CZ_0_05_mean = mean(CZ_0_05),
    CZ_0_1_mean = mean(CZ_0_1),
    CZ_0_5_mean = mean(CZ_0_5),
    CZ_1_mean = mean(CZ_1),
    NZ_0_0001_mean = mean(NZ_0_0001),
    NZ_0_01_mean = mean(NZ_0_01),
    NZ_0_05_mean = mean(NZ_0_05),
    NZ_0_1_mean = mean(NZ_0_1),
    NZ_0_5_mean = mean(NZ_0_5),
    NZ_1_mean = mean(NZ_1),
    SZ_0_0001_mean = mean(SZ_0_0001),
    SZ_0_01_mean = mean(SZ_0_01),
    SZ_0_05_mean = mean(SZ_0_05),
    SZ_0_1_mean = mean(SZ_0_1),
    SZ_0_5_mean = mean(SZ_0_5),
    SZ_1_mean = mean(SZ_1),
    HZ_0_0001_mean = mean(HZ_0_0001),
    HZ_0_01_mean = mean(HZ_0_01),
    HZ_0_05_mean = mean(HZ_0_05),
    HZ_0_1_mean = mean(HZ_0_1),
    HZ_0_5_mean = mean(HZ_0_5),
    HZ_1_mean = mean(HZ_1),
    
    AZ_0_0001_sd = sd(AZ_0_0001),
    AZ_0_01_sd = sd(AZ_0_01),
    AZ_0_05_sd = sd(AZ_0_05),
    AZ_0_1_sd = sd(AZ_0_1),
    AZ_0_5_sd = sd(AZ_0_5),
    AZ_1_sd = sd(AZ_1),
    CZ_0_0001_sd = sd(CZ_0_0001),
    CZ_0_01_sd = sd(CZ_0_01),
    CZ_0_05_sd = sd(CZ_0_05),
    CZ_0_1_sd = sd(CZ_0_1),
    CZ_0_5_sd = sd(CZ_0_5),
    CZ_1_sd = sd(CZ_1),
    NZ_0_0001_sd = sd(NZ_0_0001),
    NZ_0_01_sd = sd(NZ_0_01),
    NZ_0_05_sd = sd(NZ_0_05),
    NZ_0_1_sd = sd(NZ_0_1),
    NZ_0_5_sd = sd(NZ_0_5),
    NZ_1_sd = sd(NZ_1),
    SZ_0_0001_sd = sd(SZ_0_0001),
    SZ_0_01_sd = sd(SZ_0_01),
    SZ_0_05_sd = sd(SZ_0_05),
    SZ_0_1_sd = sd(SZ_0_1),
    SZ_0_5_sd = sd(SZ_0_5),
    SZ_1_sd = sd(SZ_1),
    HZ_0_0001_sd = sd(HZ_0_0001),
    HZ_0_01_sd = sd(HZ_0_01),
    HZ_0_05_sd = sd(HZ_0_05),
    HZ_0_1_sd = sd(HZ_0_1),
    HZ_0_5_sd = sd(HZ_0_5),
    HZ_1_sd = sd(HZ_1)
  )

PRS_means_drug <- PRS_means_drug[,2:61]
PRS_means_drug <- melt(PRS_means_drug)
PRS_means_drug <- cbind(PRS_means_drug[1:60,],
                        PRS_means_drug[61:120,2])
PRS_means_drug <- cbind(PRS_means_drug,rep(c("Escitalopram","Nortriptyline"),30))
colnames(PRS_means_drug) <- c("PRS", "Group_mean", "SD", "Drug")

ggplot(PRS_means_drug, 
       aes(x = Group_mean, y = reorder(PRS, desc(PRS)), 
           col = Drug)) +
  geom_point() + 
  geom_errorbar(aes(xmin=Group_mean-SD, 
                    xmax=Group_mean+SD), width=.2,
                position=position_dodge(.9))+
  xlab('Polygenic score') +
  ylab('Group mean')+
  theme(axis.title=element_text(size=7))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title="Were the non-randomized people assigned a drug based on their polygenic score?")+
  scale_colour_manual(labels=c("Escitalopram", "Nortriptyline"),
                      values=c("lightpink1", "yellowgreen"))

ggplot(GEND_nonrand, aes(x=HZ_0_5,color=drug))+
  geom_density()+
  labs(title="Were the non-randomized people assigned a drug based on their HZ_0_5 polygenic score?")+
  scale_colour_manual(labels=c("Escitalopram", "Nortriptyline"),
                      values=c("lightpink1", "yellowgreen"))

# Check whether polygenic score means are significantly
# different for each drug, in the non-randomized subset

t.test(
  subset(GE_esc,random==0)$HZ_0_5,
  subset(GE_nor,random==0)$HZ_0_5
)

dev.off()

##### 07b) Visualize distributions of outcomes over time #####

# Plot ZMADRS change over time for every individual

plot1 <- ggplot(GENDEP, aes(x = week, y = zmadrs)) +
  geom_line(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('ZMADRS') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ZMADRS over time

plot2 <- ggplot(GENDEP, aes(x = factor(week), y = zmadrs)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("ZMADRS") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot total ASEC change over time for each individual

plot3 <- ggplot(GENDEP, aes(x = week, y = ztotasec)) +
  geom_line(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('Total ASEC, z score') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of total side effects over time

plot4 <- ggplot(GENDEP, aes(x = factor(week), y = ztotasec)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("Total ASEC, z score") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot zf1score change over time for every individual

plot5 <- ggplot(GENDEP, aes(x = week, y = zf1score)) +
  geom_line(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('zF1score') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of zf1score over time

plot6 <- ggplot(GENDEP, aes(x = factor(week), y = zf1score)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("zF1score") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot zf2score change over time for every individual

plot7 <- ggplot(GENDEP, aes(x = week, y = zf2score)) +
  geom_line(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('zF2score') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of zf2score over time

plot8 <- ggplot(GENDEP, aes(x = factor(week), y = zf2score)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("zF2score") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot zf3score change over time for every individual

plot9 <- ggplot(GENDEP, aes(x = week, y = zf3score)) +
  geom_line(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('zF3score') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of zf3score over time

plot10 <- ggplot(GENDEP, aes(x = factor(week), y = zf3score)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("zF3score") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot suiscore change over time for every individual

plot11 <- ggplot(GENDEP, aes(x = week, y = suiscore)) +
  geom_line(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('Suicidality score') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of suiscore over time

plot12 <- ggplot(GENDEP, aes(x = factor(week), y = suiscore)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("Suicidality score") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC1 change over time for each individual

plot13 <- ggplot(GENDEP, aes(x = week, y = asec1wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec1wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC1 over time

plot14 <- ggplot(GENDEP, aes(x = factor(week), y = asec1wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec1wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC2 change over time for each individual

plot15 <- ggplot(GENDEP, aes(x = week, y = asec2wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec2wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC2 over time

plot16 <- ggplot(GENDEP, aes(x = factor(week), y = asec2wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec2wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC3 change over time for each individual

plot17 <- ggplot(GENDEP, aes(x = week, y = asec3wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec3wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC3 over time

plot18 <- ggplot(GENDEP, aes(x = factor(week), y = asec3wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec3wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC4 change over time for each individual

plot19 <- ggplot(GENDEP, aes(x = week, y = asec4wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec4wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC4 over time

plot20 <- ggplot(GENDEP, aes(x = factor(week), y = asec4wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec4wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC5 change over time for each individual

plot21 <- ggplot(GENDEP, aes(x = week, y = asec5wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec5wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC5 over time

plot22 <- ggplot(GENDEP, aes(x = factor(week), y = asec5wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec5wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC6 change over time for each individual

plot23 <- ggplot(GENDEP, aes(x = week, y = asec6wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec6wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC6 over time

plot24 <- ggplot(GENDEP, aes(x = factor(week), y = asec6wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec6wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC7 change over time for each individual

plot25 <- ggplot(GENDEP, aes(x = week, y = asec7wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec7wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC7 over time

plot26 <- ggplot(GENDEP, aes(x = factor(week), y = asec7wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec7wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC8 change over time for each individual

plot27 <- ggplot(GENDEP, aes(x = week, y = asec8wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec8wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC8 over time

plot28 <- ggplot(GENDEP, aes(x = factor(week), y = asec8wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec8wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC9 change over time for each individual

plot29 <- ggplot(GENDEP, aes(x = week, y = asec9wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec9wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC9 over time

plot30 <- ggplot(GENDEP, aes(x = factor(week), y = asec9wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec9wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC10 change over time for each individual

plot31 <- ggplot(GENDEP, aes(x = week, y = asec10wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec10wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC10 over time

plot32 <- ggplot(GENDEP, aes(x = factor(week), y = asec10wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec10wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC11 change over time for each individual

plot33 <- ggplot(GENDEP, aes(x = week, y = asec11wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec11wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC11 over time

plot34 <- ggplot(GENDEP, aes(x = factor(week), y = asec11wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec11wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC12 change over time for each individual

plot35 <- ggplot(GENDEP, aes(x = week, y = asec12wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec12wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC12 over time

plot36 <- ggplot(GENDEP, aes(x = factor(week), y = asec12wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec12wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC13 change over time for each individual

plot37 <- ggplot(GENDEP, aes(x = week, y = asec13wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec13wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC13 over time

plot38 <- ggplot(GENDEP, aes(x = factor(week), y = asec13wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec13wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC14 change over time for each individual

plot39 <- ggplot(GENDEP, aes(x = week, y = asec14wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec14wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC14 over time

plot40 <- ggplot(GENDEP, aes(x = factor(week), y = asec14wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec14wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC15 change over time for each individual

plot41 <- ggplot(GENDEP, aes(x = week, y = asec15wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec15wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC15 over time

plot42 <- ggplot(GENDEP, aes(x = factor(week), y = asec15wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec15wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC16 change over time for each individual

plot43 <- ggplot(GENDEP, aes(x = week, y = asec16wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec16wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC16 over time

plot44 <- ggplot(GENDEP, aes(x = factor(week), y = asec16wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec16wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC17 change over time for each individual

plot45 <- ggplot(GENDEP, aes(x = week, y = asec17wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec17wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC17 over time

plot46 <- ggplot(GENDEP, aes(x = factor(week), y = asec17wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec17wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC18 change over time for each individual

plot47 <- ggplot(GENDEP, aes(x = week, y = asec18wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec18wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC18 over time

plot48 <- ggplot(GENDEP, aes(x = factor(week), y = asec18wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec18wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC19 change over time for each individual

plot49 <- ggplot(GENDEP, aes(x = week, y = asec19wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec19wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC19 over time

plot50 <- ggplot(GENDEP, aes(x = factor(week), y = asec19wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec19wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC20 change over time for each individual

plot51 <- ggplot(GENDEP, aes(x = week, y = asec20wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec20wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC20 over time

plot52 <- ggplot(GENDEP, aes(x = factor(week), y = asec20wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec20wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Plot ASEC21 change over time for each individual

plot53 <- ggplot(GENDEP, aes(x = week, y = asec21wk)) +
  geom_jitter(aes(group = subjectid, colour = Drug), alpha = 0.8) +
  scale_color_manual(values=c("lightpink1", "greenyellow"))+
  geom_smooth(se = FALSE, size = 1, colour = "black") +
  xlab('Week') +
  ylab('asec21wk') +
  theme_bw(base_size = 14)+
  theme(legend.position = "none")

# Visualize distributions of ASEC21 over time

plot54 <- ggplot(GENDEP, aes(x = factor(week), y = asec21wk)) + 
  geom_violin(aes(fill=Drug))+ 
  scale_fill_manual(values=c("lightpink1", "greenyellow"))+
  xlab("Week") + 
  ylab("asec21wk") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "none")

# Produce pdf

pdf( "Exploratory plots.pdf" )

grid.arrange( plot1, plot2 )
grid.arrange( plot3, plot4 )
grid.arrange( plot5, plot6 )
grid.arrange( plot7, plot8 )
grid.arrange( plot9, plot10 )
grid.arrange( plot11, plot12 )
grid.arrange( plot13, plot14 )
grid.arrange( plot15, plot16 )
grid.arrange( plot17, plot18 )
grid.arrange( plot19, plot20 )
grid.arrange( plot21, plot22 )
grid.arrange( plot23, plot24 )
grid.arrange( plot25, plot26 )
grid.arrange( plot27, plot28 )
grid.arrange( plot29, plot30 )
grid.arrange( plot31, plot32 )
grid.arrange( plot33, plot34 )
grid.arrange( plot35, plot36 )
grid.arrange( plot37, plot38 )
grid.arrange( plot39, plot40 )
grid.arrange( plot41, plot42 )
grid.arrange( plot43, plot44 )
grid.arrange( plot45, plot46 )
grid.arrange( plot47, plot48 )
grid.arrange( plot49, plot50 )
grid.arrange( plot51, plot52 )
grid.arrange( plot53, plot54 )

dev.off()

##### Correlation matrix of PRSs #####

PolSco <- cbind(GEND_end$AZ_0_0001,
                GEND_end$AZ_0_01,
                GEND_end$AZ_0_05,
                GEND_end$AZ_0_1,
                GEND_end$AZ_0_5,
                GEND_end$AZ_1,
                GEND_end$CZ_0_0001,
                GEND_end$CZ_0_01,
                GEND_end$CZ_0_05,
                GEND_end$CZ_0_1,
                GEND_end$CZ_0_5,
                GEND_end$CZ_1,
                GEND_end$NZ_0_0001,
                GEND_end$NZ_0_01,
                GEND_end$NZ_0_05,
                GEND_end$NZ_0_1,
                GEND_end$NZ_0_5,
                GEND_end$NZ_1,
                GEND_end$SZ_0_0001,
                GEND_end$SZ_0_01,
                GEND_end$SZ_0_05,
                GEND_end$SZ_0_1,
                GEND_end$SZ_0_5,
                GEND_end$SZ_1,
                GEND_end$HZ_0_0001,
                GEND_end$HZ_0_01,
                GEND_end$HZ_0_05,
                GEND_end$HZ_0_1,
                GEND_end$HZ_0_5,
                GEND_end$HZ_1
)
colnames(PolSco)<-PRSs

correl <- cor(PolSco)
diag(correl) <- NA
range(correl, na.rm = TRUE)
scale_limits <- range(correl, na.rm = TRUE)
ggcorrplot(correl,
           type="upper", # get rid of redundancy
           hc.order = TRUE # hierarchical clustering
)+
  scale_fill_gradient2(
    limits = scale_limits,
    low = "blue",
    mid = "white",
    high = "red"
  )

sink("PRS_correlations.csv")
print(correl)
sink()

##### Correlation matrix of symptom dimensions #####

Dims <- cbind(GENDEP$zf1score,
                GENDEP$zf2score,
                GENDEP$zf3score,
                GENDEP$suiscore
)
colnames(Dims)<-c("Mood","Cognitive","Neurovegetative","Suicidality")

correl2 <- cor(Dims)
diag(correl2) <- NA
range(correl2, na.rm = TRUE)
scale_limits <- range(correl2, na.rm = TRUE)
ggcorrplot(correl2,
           type="upper", # get rid of redundancy
           hc.order = TRUE # hierarchical clustering
)+
  scale_fill_gradient2(
    limits = scale_limits,
    low = "blue",
    mid = "white",
    high = "red"
  )

sink("Symptom_correlations.csv")
print(correl2)
sink()

# Looking at dims over time (change GENDEP row indices to get a different
# subject)

ggplot(GENDEP[52:61,],aes(x=week, y=madrs, group=subjectid))+
  geom_line(aes(y=f1score,col="f1score"))+
  geom_line(aes(y=f2score,col="f2score"))+
  geom_line(aes(y=f3score,col="f3score"))+
  geom_line(aes(y=suiscore,col="suiscore"))

##### Some quick linear regressions between PS and phenotype #####

# Neuroticism PS predicts cognitive score

summary(lm(GENDEP$f2score~GENDEP$NZ_0_0001))
summary(lm(GENDEP$f2score~GENDEP$NZ_1)) # 2nd best p
summary(lm(GENDEP$f2score~GENDEP$NZ_0_5)) # 3rd best p

# Anorexia PS does NOT predict bmi_out

summary(lm(GEND_end$bmi_out~GEND_end$AZ_0_0001))

# Height PS predicts height

summary(lm(GEND_end$height~GEND_end$HZ_0_0001))


##### Get PS means for each group (rand/non; esc/nor) #####

# Keep one obs per individual

R_NR_esc <- GENDEP_escit %>% 
  distinct(subjectid, .keep_all = TRUE)

R_NR_nor <- GENDEP_nortrip %>% 
  distinct(subjectid, .keep_all = TRUE)

R_esc <- GEND_rand_e %>% 
  distinct(subjectid, .keep_all = TRUE)

R_nor <- GEND_rand_n %>% 
  distinct(subjectid, .keep_all = TRUE)

NR_esc <- GEND_nonrand_e %>% 
  distinct(subjectid, .keep_all = TRUE)

NR_nor <- GEND_nonrand_n %>% 
  distinct(subjectid, .keep_all = TRUE)

# Initialise empty data frame

PSs_by_group <- data.frame(
  matrix(
    nrow = 180,
    ncol = 4
  )
)

colnames( PSs_by_group ) <- c( "Group",
                               "PS",
                               "Mean",
                               "SD"
                               )

PSs_by_group$Group <- c(
  rep( "allEscit", 30),
  rep( "allNortr", 30),
  rep( "randEscit", 30),
  rep( "randNortr", 30),
  rep( "nonrandEscit", 30),
  rep( "nonrandNortr", 30)
)

PTs <- c(
  rep(c(
    "0.0001",
    "0.01",
    "0.05",
    "0.1",
    "0.5",
    "1.0"
    ),
    6
    )
)

PSs_by_group$PT <- c(
  rep( PTs, 5 )
)

Traits <- c(
  rep("Anorexia",6),
  rep("ReactTime",6),
  rep("Neuroticism",6),
  rep("Schizophrenia",6),
  rep("Height",6)
)

PSs_by_group$Trait <- c(
  rep( Traits, 6 )
)

# Get means and SDs

for (i in 1:length(PRSs)){
  
  PSs_by_group[i,3] <- mean( R_NR_esc[,PRSs[i]] )
  PSs_by_group[i,4] <- sd( R_NR_esc[,PRSs[i]] )
  PSs_by_group[i+30,3] <- mean( R_NR_nor[,PRSs[i]] )
  PSs_by_group[i+30,4] <- sd( R_NR_nor[,PRSs[i]] )
  PSs_by_group[i+60,3] <- mean( R_esc[,PRSs[i]] )
  PSs_by_group[i+60,4] <- sd( R_esc[,PRSs[i]] )
  PSs_by_group[i+90,3] <- mean( R_nor[,PRSs[i]] )
  PSs_by_group[i+90,4] <- sd( R_nor[,PRSs[i]] )
  PSs_by_group[i+120,3] <- mean( NR_esc[,PRSs[i]] )
  PSs_by_group[i+120,4] <- sd( NR_esc[,PRSs[i]] )
  PSs_by_group[i+150,3] <- mean( NR_nor[,PRSs[i]] )
  PSs_by_group[i+150,4] <- sd( NR_nor[,PRSs[i]] )

}

PSs_by_group <- PSs_by_group %>% 
  mutate( SDlow = Mean-SD)

PSs_by_group <- PSs_by_group %>% 
  mutate( SDhi = Mean+SD)

# Plot

plot1 <- ggplot( subset(PSs_by_group,Trait == "Anorexia"),
               aes( x = Group, y = Mean,
                           col = Group,
                           ymin = SDlow,
                           ymax = SDhi))+
  geom_point(size=4)+
  geom_linerange(size=1.5)+
  facet_wrap(vars(PT))+
  labs(title = "Z-scored anorexia polygenic score means, coloured by group, separated by p-value threshold")+
  theme( axis.text.x = element_blank(),
         axis.ticks = element_blank())+
  scale_color_manual(values=c(
          "#EC407A",
          "#66BB6A",
          "#D81B60",
          "#43A047",
          "#AD1457",
          "#2E7D32"
        )
  )

plot2 <- ggplot( subset(PSs_by_group,Trait == "ReactTime"),
        aes( x = Group, y = Mean,
             col = Group,
             ymin = SDlow,
             ymax = SDhi))+
  geom_point(size=4)+
  geom_linerange(size=1.5)+
  facet_wrap(vars(PT))+
  labs(title = "Z-scored reaction time polygenic score means, coloured by group, separated by p-value threshold")+
  theme( axis.text.x = element_blank(),
         axis.ticks = element_blank())+
  scale_color_manual(values=c(
    "#EC407A",
    "#66BB6A",
    "#D81B60",
    "#43A047",
    "#AD1457",
    "#2E7D32"
  )
  )

plot3 <- ggplot( subset(PSs_by_group,Trait == "Neuroticism"),
        aes( x = Group, y = Mean,
             col = Group,
             ymin = SDlow,
             ymax = SDhi))+
  geom_point(size=4)+
  geom_linerange(size=1.5)+
  facet_wrap(vars(PT))+
  labs(title = "Z-scored neuroticism polygenic score means, coloured by group, separated by p-value threshold")+
  theme( axis.text.x = element_blank(),
         axis.ticks = element_blank())+
  scale_color_manual(values=c(
    "#EC407A",
    "#66BB6A",
    "#D81B60",
    "#43A047",
    "#AD1457",
    "#2E7D32"
  )
  )

plot4 <- ggplot( subset(PSs_by_group,Trait == "Schizophrenia"),
        aes( x = Group, y = Mean,
             col = Group,
             ymin = SDlow,
             ymax = SDhi))+
  geom_point(size=4)+
  geom_linerange(size=1.5)+
  facet_wrap(vars(PT))+
  labs(title = "Z-scored schizophrenia polygenic score means, coloured by group, separated by p-value threshold")+
  theme( axis.text.x = element_blank(),
         axis.ticks = element_blank())+
  scale_color_manual(values=c(
    "#EC407A",
    "#66BB6A",
    "#D81B60",
    "#43A047",
    "#AD1457",
    "#2E7D32"
  )
  )

plot5 <- ggplot( subset(PSs_by_group,Trait == "Height"),
        aes( x = Group, y = Mean,
             col = Group,
             ymin = SDlow,
             ymax = SDhi))+
  geom_point(size=4)+
  geom_linerange(size=1.5)+
  facet_wrap(vars(PT))+
  labs(title = "Z-scored height polygenic score means, coloured by group, separated by p-value threshold")+
  theme( axis.text.x = element_blank(),
         axis.ticks = element_blank())+
  scale_color_manual(values=c(
    "#EC407A",
    "#66BB6A",
    "#D81B60",
    "#43A047",
    "#AD1457",
    "#2E7D32"
  )
  )

pdf( "Polygenic_score_means.pdf" )
plot1
plot2
plot3
plot4
plot5
dev.off()

##### 08) ANALYSE EFFECTS OF POLYGENIC SCORE ON VARIOUS BASELINES #####

# Dependent variables:

### blmd: baseline MADRS score
### bltotasec: baseline total ASEC side effect burden
### blfm: baseline mood score (f1)
### blfc: baseline cognitive score (f2)
### blfn: baseline neurovegetative score (f3)
### blsui: baseline suicidality score

##### 08a) Build a model for each baseline #####

mod1 <- lmer( blmd ~
                SZ_0_0001 +
                age + sex + (1|centreid), 
              data=GEND_end )

mod2 <- lmer( blmd ~
                SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                age + sex + (1|centreid), 
              data=GEND_end )

# Principal components cause "singular fit" issues
# Don't include PCs for blmd model

mod1 # n = 635

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( bltotasec ~
                SZ_0_0001 +
                age + sex + (1|centreid), 
              data=GEND_end )

mod2 <- lmer( bltotasec ~
                SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                age + sex + (1|centreid), 
              data=GEND_end )

anova(mod1, mod2) # Principal components improve model a bit
drop1(mod2) # Try with just zPC1

mod3 <- lmer( bltotasec ~
                SZ_0_0001 +
                zPC1 + 
                age + sex + (1|centreid), 
              data=GEND_end )

anova(mod1, mod2, mod3) # PC1 is worth adding

mod3 # n = 509

hist(residuals(mod3)) 
residplot(mod3, newwd=FALSE)

mod1 <- lmer( blfm ~
                SZ_0_0001 +
                age + sex + (1|centreid), 
              data=GEND_end )

mod2 <- lmer( blfm ~
                SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                age + sex + (1|centreid), 
              data=GEND_end )

anova(mod1, mod2) # Principal components improve model significantly

mod2 # n = 635

hist(residuals(mod2)) 
residplot(mod2, newwd=FALSE)

mod1 <- lmer( blfc ~
                SZ_0_0001 +
                age + sex + (1|centreid), 
              data=GEND_end )

mod2 <- lmer( blfc ~
                SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                age + sex + (1|centreid), 
              data=GEND_end )

anova(mod1, mod2) # Principal components make model worse

mod1 # n = 635

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( blfn ~
                SZ_0_0001 +
                age + sex + (1|centreid), 
              data=GEND_end )

mod2 <- lmer( blfn ~
                SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                age + sex + (1|centreid), 
              data=GEND_end )

anova(mod1, mod2) # Principal components improve model significantly

mod2 # n = 635

hist(residuals(mod2)) 
residplot(mod2, newwd=FALSE)

mod1 <- lmer( blsui ~
                SZ_0_0001 +
                age + sex + (1|centreid), 
              data=GEND_end )

mod2 <- lmer( blsui ~
                SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                age + sex + (1|centreid), 
              data=GEND_end )

anova(mod1, mod2) # Principal components improve model significantly

mod2 # n = 635

hist(residuals(mod2)) 
residplot(mod2, newwd=FALSE)

##### 08b) Fit mixed effect models for baselines #####

# Make lists to store models in

blmd_models <- vector(mode = "list", length = length(PRSs_31))
bltotasec_models <- vector(mode = "list", length = length(PRSs_31))
blfm_models <- vector(mode = "list", length = length(PRSs_31))
blfc_models <- vector(mode = "list", length = length(PRSs_31))
blfn_models <- vector(mode = "list", length = length(PRSs_31))
blsui_models <- vector(mode = "list", length = length(PRSs_31))

blmd_models_2 <- vector(mode = "list", length = length(PRSs_31))
bltotasec_models_2 <- vector(mode = "list", length = length(PRSs_31))
blfm_models_2 <- vector(mode = "list", length = length(PRSs_31))
blfc_models_2 <- vector(mode = "list", length = length(PRSs_31))
blfn_models_2 <- vector(mode = "list", length = length(PRSs_31))
blsui_models_2 <- vector(mode = "list", length = length(PRSs_31))

# Use "GEND_end" as it has one observation per participant
# Make a special subset for analysing baseline total ASEC

GEND_end_bltotasec <- subset( GEND_end, !is.na(bltotasec) )

# Covariates are age, sex and random effect: recruiting centre
# Some models have principal components, where this was beneficial

for ( i in 1:length( PRSs ) ){
  
  blmd_models[[i]] <- lmer( blmd ~
                              GEND_end[,PRSs[i]] +
                              # No principal components
                              age + sex + (1|centreid),
                            data=GEND_end, REML=FALSE )
  
  removd <- romr.fnc(blmd_models[[i]], GEND_end, trim = 2.5)$data # remove influential outliers
  
  blmd_models_2[[i]] <- lmer( blmd ~
                              removd[,PRSs[i]] +
                              # No principal components
                              age + sex + (1|centreid),
                            data=removd, REML=FALSE )  
  
  ####
  
  bltotasec_models[[i]] <- lmer( bltotasec ~ 
                                   GEND_end_bltotasec[,PRSs[i]] + 
                                   zPC1 + # Omit PC2 and PC3
                                   age + sex + (1|centreid), 
                                 data=GEND_end_bltotasec, REML=FALSE )
  
  removd <- romr.fnc(bltotasec_models[[i]], GEND_end_bltotasec, trim = 2.5)$data # remove influential outliers
  
  bltotasec_models_2[[i]] <- lmer( bltotasec ~ 
                                   removd[,PRSs[i]] + 
                                   zPC1 + # Omit PC2 and PC3
                                   age + sex + (1|centreid), 
                                 data=removd, REML=FALSE )
  
  ####
  
  blfm_models[[i]] <- lmer( blfm ~ 
                              GEND_end[,PRSs[i]] + 
                              zPC1 + zPC2 + zPC3 +
                              age + sex + (1|centreid), 
                            data=GEND_end, REML=FALSE )
  
  removd <- romr.fnc(blfm_models[[i]], GEND_end, trim = 2.5)$data # remove influential outliers
  
  blfm_models_2[[i]] <- lmer( blfm ~ 
                              removd[,PRSs[i]] + 
                              zPC1 + zPC2 + zPC3 +
                              age + sex + (1|centreid), 
                            data=removd, REML=FALSE )
  
  ####
  
  blfc_models[[i]] <- lmer( blfc ~ 
                              GEND_end[,PRSs[i]] + 
                              # No principal components
                              age + sex + (1|centreid), 
                            data=GEND_end, REML=FALSE )
  
  removd <- romr.fnc(blfc_models[[i]], GEND_end, trim = 2.5)$data # remove influential outliers
  
  blfc_models_2[[i]] <- lmer( blfc ~ 
                              removd[,PRSs[i]] + 
                              # No principal components
                              age + sex + (1|centreid), 
                            data=removd, REML=FALSE )
  
  ####
  
  blfn_models[[i]] <- lmer( blfn ~ 
                              GEND_end[,PRSs[i]] + 
                              zPC1 + zPC2 + zPC3 +
                              age + sex + (1|centreid), 
                            data=GEND_end, REML=FALSE )
  
  removd <- romr.fnc(blfn_models[[i]], GEND_end, trim = 2.5)$data # remove influential outliers
  
  blfn_models_2[[i]] <- lmer( blfn ~ 
                              removd[,PRSs[i]] + 
                              zPC1 + zPC2 + zPC3 +
                              age + sex + (1|centreid), 
                            data=removd, REML=FALSE )
  
  ####
  
  blsui_models[[i]] <- lmer( blsui ~ 
                               GEND_end[,PRSs[i]] + 
                               zPC1 + zPC2 + zPC3 +
                               age + sex + (1|centreid), 
                             data=GEND_end, REML=FALSE )
  
  removd <- romr.fnc(blsui_models[[i]], GEND_end, trim = 2.5)$data # remove influential outliers
  
  blsui_models_2[[i]] <- lmer( blsui ~ 
                               removd[,PRSs[i]] + 
                               zPC1 + zPC2 + zPC3 +
                               age + sex + (1|centreid), 
                             data=removd, REML=FALSE )
}

# Add non-PRS models for comparison

blmd_models[[31]] <- lmer( blmd ~
                             age + sex + (1|centreid), 
                           data=GEND_end, REML=FALSE )

removd <- romr.fnc(blmd_models[[31]], GEND_end, trim = 2.5)$data # remove influential outliers

blmd_models_2[[31]] <- lmer( blmd ~
                             age + sex + (1|centreid), 
                           data=removd, REML=FALSE )

####

bltotasec_models[[31]] <- lmer( bltotasec ~ 
                                  zPC1 + 
                                  age + sex + (1|centreid), 
                                data=GEND_end_bltotasec, REML=FALSE )

removd <- romr.fnc(bltotasec_models[[31]], GEND_end_bltotasec, trim = 2.5)$data # remove influential outliers

bltotasec_models_2[[31]] <- lmer( bltotasec ~
                                    zPC1 +
                                    age + sex + (1|centreid), 
                                  data=removd, REML=FALSE )

####

blfm_models[[31]] <- lmer( blfm ~ 
                             zPC1 + zPC2 + zPC3 +
                             age + sex + (1|centreid), 
                           data=GEND_end, REML=FALSE )

removd <- romr.fnc(blfm_models[[31]], GEND_end, trim = 2.5)$data # remove influential outliers

blfm_models_2[[31]] <- lmer( blfm ~
                               zPC1 + zPC2 + zPC3 +
                               age + sex + (1|centreid), 
                             data=removd, REML=FALSE )

####

blfc_models[[31]] <- lmer( blfc ~ 
                             age + sex + (1|centreid), 
                           data=GEND_end, REML=FALSE )

removd <- romr.fnc(blfc_models[[31]], GEND_end, trim = 2.5)$data # remove influential outliers

blfc_models_2[[31]] <- lmer( blfc ~
                               age + sex + (1|centreid), 
                             data=removd, REML=FALSE )

####

blfn_models[[31]] <- lmer( blfn ~ 
                             zPC1 + zPC2 + zPC3 +
                             age + sex + (1|centreid), 
                           data=GEND_end, REML=FALSE )

removd <- romr.fnc(blfn_models[[31]], GEND_end, trim = 2.5)$data # remove influential outliers

blfn_models_2[[31]] <- lmer( blfn ~
                               zPC1 + zPC2 + zPC3 +
                               age + sex + (1|centreid), 
                             data=removd, REML=FALSE )

####

blsui_models[[31]] <- lmer( blsui ~ 
                              zPC1 + zPC2 + zPC3 +
                              age + sex + (1|centreid), 
                            data=GEND_end, REML=FALSE )

removd <- romr.fnc(blsui_models[[31]], GEND_end, trim = 2.5)$data # remove influential outliers

blsui_models_2[[31]] <- lmer( blsui ~
                               zPC1 + zPC2 + zPC3 +
                               age + sex + (1|centreid), 
                             data=removd, REML=FALSE )

names(blmd_models) <- c(PRSs_31)
names(bltotasec_models) <- c(PRSs_31)
names(blfm_models) <- c(PRSs_31)
names(blfc_models) <- c(PRSs_31)
names(blfn_models) <- c(PRSs_31)
names(blsui_models) <- c(PRSs_31)

names(blmd_models_2) <- c(PRSs_31)
names(bltotasec_models_2) <- c(PRSs_31)
names(blfm_models_2) <- c(PRSs_31)
names(blfc_models_2) <- c(PRSs_31)
names(blfn_models_2) <- c(PRSs_31)
names(blsui_models_2) <- c(PRSs_31)

##### 08c) Get betas, standard errors and p-values for baselines #####

# broom.mixed::tidy(blmd_models[[1]], conf.int = TRUE)

blmd_betas <- vector(mode = "list", length = length(PRSs_31))
bltotasec_betas <- vector(mode = "list", length = length(PRSs_31))
blfm_betas <- vector(mode = "list", length = length(PRSs_31))
blfc_betas <- vector(mode = "list", length = length(PRSs_31))
blfn_betas <- vector(mode = "list", length = length(PRSs_31))
blsui_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(blmd_models)){
  blmd_betas[i] <- list(summary(blmd_models[[i]])$coefficients[2,])
}
for (i in 1:length(bltotasec_models)){
  bltotasec_betas[i] <- list(summary(bltotasec_models[[i]])$coefficients[2,])
}
for (i in 1:length(blfm_models)){
  blfm_betas[i] <- list(summary(blfm_models[[i]])$coefficients[2,])
}
for (i in 1:length(blfc_models)){
  blfc_betas[i] <- list(summary(blfc_models[[i]])$coefficients[2,])
}
for (i in 1:length(blfn_models)){
  blfn_betas[i] <- list(summary(blfn_models[[i]])$coefficients[2,])
}
for (i in 1:length(blsui_models)){
  blsui_betas[i] <- list(summary(blsui_models[[i]])$coefficients[2,])
}

names(blmd_betas) <- PRSs_31
names(bltotasec_betas) <- PRSs_31
names(blfm_betas) <- PRSs_31
names(blfc_betas) <- PRSs_31
names(blfn_betas) <- PRSs_31
names(blsui_betas) <- PRSs_31

blmd_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(blmd_betas))[,1],
    t(as.data.frame(blmd_betas))[,2],
    t(as.data.frame(blmd_betas))[,5]
  )
)

colnames( blmd_beta_matrix ) <- c(
  "Beta", "SE","pval"
)

bltotasec_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(bltotasec_betas))[,1],
    t(as.data.frame(bltotasec_betas))[,2],
    t(as.data.frame(bltotasec_betas))[,5]
  )
)

colnames( bltotasec_beta_matrix ) <- c(
  "Beta", "SE","pval"
)

blfm_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(blfm_betas))[,1],
    t(as.data.frame(blfm_betas))[,2],
    t(as.data.frame(blfm_betas))[,5]
  )
)

colnames( blfm_beta_matrix ) <- c(
  "Beta", "SE","pval"
)

blfc_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(blfc_betas))[,1],
    t(as.data.frame(blfc_betas))[,2],
    t(as.data.frame(blfc_betas))[,5]
  )
)

colnames( blfc_beta_matrix ) <- c(
  "Beta", "SE","pval"
)

blfn_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(blfn_betas))[,1],
    t(as.data.frame(blfn_betas))[,2],
    t(as.data.frame(blfn_betas))[,5]
  )
)

colnames( blfn_beta_matrix ) <- c(
  "Beta", "SE","pval"
)

blsui_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(blsui_betas))[,1],
    t(as.data.frame(blsui_betas))[,2],
    t(as.data.frame(blsui_betas))[,5]
  )
)

colnames( blsui_beta_matrix ) <- c(
  "Beta", "SE","pval"
)

blmd_beta_matrix <- blmd_beta_matrix[1:30,]
bltotasec_beta_matrix <- bltotasec_beta_matrix[1:30,]
blfm_beta_matrix <- blfm_beta_matrix[1:30,]
blfc_beta_matrix <- blfc_beta_matrix[1:30,]
blfn_beta_matrix <- blfn_beta_matrix[1:30,]
blsui_beta_matrix <- blsui_beta_matrix[1:30,]

####

blmd_betas_2 <- vector(mode = "list", length = length(PRSs_31))
bltotasec_betas_2 <- vector(mode = "list", length = length(PRSs_31))
blfm_betas_2 <- vector(mode = "list", length = length(PRSs_31))
blfc_betas_2 <- vector(mode = "list", length = length(PRSs_31))
blfn_betas_2 <- vector(mode = "list", length = length(PRSs_31))
blsui_betas_2 <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(blmd_models_2)){
  blmd_betas_2[i] <- list(summary(blmd_models_2[[i]])$coefficients[2,])
}
for (i in 1:length(bltotasec_models_2)){
  bltotasec_betas_2[i] <- list(summary(bltotasec_models_2[[i]])$coefficients[2,])
}
for (i in 1:length(blfm_models_2)){
  blfm_betas_2[i] <- list(summary(blfm_models_2[[i]])$coefficients[2,])
}
for (i in 1:length(blfc_models_2)){
  blfc_betas_2[i] <- list(summary(blfc_models_2[[i]])$coefficients[2,])
}
for (i in 1:length(blfn_models_2)){
  blfn_betas_2[i] <- list(summary(blfn_models_2[[i]])$coefficients[2,])
}
for (i in 1:length(blsui_models_2)){
  blsui_betas_2[i] <- list(summary(blsui_models_2[[i]])$coefficients[2,])
}

names(blmd_betas_2) <- PRSs_31
names(bltotasec_betas_2) <- PRSs_31
names(blfm_betas_2) <- PRSs_31
names(blfc_betas_2) <- PRSs_31
names(blfn_betas_2) <- PRSs_31
names(blsui_betas_2) <- PRSs_31

blmd_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(blmd_betas_2))[,1],
    t(as.data.frame(blmd_betas_2))[,2],
    t(as.data.frame(blmd_betas_2))[,5]
  )
)

colnames( blmd_beta_matrix_2 ) <- c(
  "Beta", "SE","pval"
)

bltotasec_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(bltotasec_betas_2))[,1],
    t(as.data.frame(bltotasec_betas_2))[,2],
    t(as.data.frame(bltotasec_betas_2))[,5]
  )
)

colnames( bltotasec_beta_matrix_2 ) <- c(
  "Beta", "SE","pval"
)

blfm_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(blfm_betas_2))[,1],
    t(as.data.frame(blfm_betas_2))[,2],
    t(as.data.frame(blfm_betas_2))[,5]
  )
)

colnames( blfm_beta_matrix_2 ) <- c(
  "Beta", "SE","pval"
)

blfc_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(blfc_betas_2))[,1],
    t(as.data.frame(blfc_betas_2))[,2],
    t(as.data.frame(blfc_betas_2))[,5]
  )
)

colnames( blfc_beta_matrix_2 ) <- c(
  "Beta", "SE","pval"
)

blfn_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(blfn_betas_2))[,1],
    t(as.data.frame(blfn_betas_2))[,2],
    t(as.data.frame(blfn_betas_2))[,5]
  )
)

colnames( blfn_beta_matrix_2 ) <- c(
  "Beta", "SE","pval"
)

blsui_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(blsui_betas_2))[,1],
    t(as.data.frame(blsui_betas_2))[,2],
    t(as.data.frame(blsui_betas_2))[,5]
  )
)

colnames( blsui_beta_matrix_2 ) <- c(
  "Beta", "SE","pval"
)

blmd_beta_matrix_2 <- blmd_beta_matrix_2[1:30,]
bltotasec_beta_matrix_2 <- bltotasec_beta_matrix_2[1:30,]
blfm_beta_matrix_2 <- blfm_beta_matrix_2[1:30,]
blfc_beta_matrix_2 <- blfc_beta_matrix_2[1:30,]
blfn_beta_matrix_2 <- blfn_beta_matrix_2[1:30,]
blsui_beta_matrix_2 <- blsui_beta_matrix_2[1:30,]

##### 08d) Get PRSs' explained variances for baselines #####

blmd_expvar <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
  c("R2m","R2c")))

bltotasec_expvar <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blfm_expvar <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blfc_expvar <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blfn_expvar <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blsui_expvar <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

for (i in 1:length(PRSs_31)){
  blmd_expvar[i,] <- r.squaredGLMM(blmd_models[[i]])
  bltotasec_expvar[i,] <- r.squaredGLMM(bltotasec_models[[i]])
  blfm_expvar[i,] <- r.squaredGLMM(blfm_models[[i]])
  blfc_expvar[i,] <- r.squaredGLMM(blfc_models[[i]])
  blfn_expvar[i,] <- r.squaredGLMM(blfn_models[[i]])
  blsui_expvar[i,] <- r.squaredGLMM(blsui_models[[i]])
}

####

blmd_expvar_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

bltotasec_expvar_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blfm_expvar_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blfc_expvar_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blfn_expvar_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

blsui_expvar_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(PRSs_31),
                c("R2m","R2c")))

for (i in 1:length(PRSs_31)){
  blmd_expvar_2[i,] <- r.squaredGLMM(blmd_models_2[[i]])
  bltotasec_expvar_2[i,] <- r.squaredGLMM(bltotasec_models_2[[i]])
  blfm_expvar_2[i,] <- r.squaredGLMM(blfm_models_2[[i]])
  blfc_expvar_2[i,] <- r.squaredGLMM(blfc_models_2[[i]])
  blfn_expvar_2[i,] <- r.squaredGLMM(blfn_models_2[[i]])
  blsui_expvar_2[i,] <- r.squaredGLMM(blsui_models_2[[i]])
}

##### 08e) Make output files for baselines #####

sink( "blmd.csv")

cat("blmd betas\n,")
write.table( blmd_beta_matrix, col.names=TRUE, sep="," )
cat("\n,")
write.table( blmd_expvar, col.names=TRUE, sep="," )

sink()

sink( "bltotasec.csv")

cat("bltotasec betas\n,")
write.table( bltotasec_beta_matrix, col.names=TRUE, sep="," )
cat("\n,")
write.table( bltotasec_expvar, col.names=TRUE, sep="," )

sink()

sink( "blfm.csv")

cat("blfm betas\n,")
write.table( blfm_beta_matrix, col.names=TRUE, sep="," )
cat("\n,")
write.table( blfm_expvar, col.names=TRUE, sep="," )

sink()

sink( "blfc.csv")

cat("blfc betas\n,")
write.table( blfc_beta_matrix, col.names=TRUE, sep="," )
cat("\n,")
write.table( blfc_expvar, col.names=TRUE, sep="," )

sink()

sink( "blfn.csv")

cat("blfn betas\n,")
write.table( blfn_beta_matrix, col.names=TRUE, sep="," )
cat("\n,")
write.table( blfn_expvar, col.names=TRUE, sep="," )

sink()

sink( "blsui.csv")

cat("blsui betas\n,")
write.table( blsui_beta_matrix, col.names=TRUE, sep="," )
cat("\n,")
write.table( blsui_expvar, col.names=TRUE, sep="," )

sink()

####

sink( "blmd_no_inf_outl.csv")

cat("blmd betas\n,")
write.table( blmd_beta_matrix_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( blmd_expvar_2, col.names=TRUE, sep="," )

sink()

sink( "bltotasec_no_inf_outl.csv")

cat("bltotasec betas\n,")
write.table( bltotasec_beta_matrix_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( bltotasec_expvar_2, col.names=TRUE, sep="," )

sink()

sink( "blfm_no_inf_outl.csv")

cat("blfm betas\n,")
write.table( blfm_beta_matrix_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( blfm_expvar_2, col.names=TRUE, sep="," )

sink()

sink( "blfc_no_inf_outl.csv")

cat("blfc betas\n,")
write.table( blfc_beta_matrix_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( blfc_expvar_2, col.names=TRUE, sep="," )

sink()

sink( "blfn_no_inf_outl.csv")

cat("blfn betas\n,")
write.table( blfn_beta_matrix_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( blfn_expvar_2, col.names=TRUE, sep="," )

sink()

sink( "blsui_no_inf_outl.csv")

cat("blsui betas\n,")
write.table( blsui_beta_matrix_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( blsui_expvar_2, col.names=TRUE, sep="," )

sink()

##### 08f) Make heatmap of effect sizes for baselines #####

baseline_betas_only <- cbind(
  blmd_beta_matrix[,1],
  bltotasec_beta_matrix[,1],
  blfm_beta_matrix[,1],
  blfc_beta_matrix[,1],
  blfn_beta_matrix[,1],
  blsui_beta_matrix[,1]
)

colnames(baseline_betas_only) <- c("blmd", "bltotasec", "blfm",
                                   "blfc", "blfn", "blsui")

baseline_SEs_only <- cbind(
  blmd_beta_matrix[,2],
  bltotasec_beta_matrix[,2],
  blfm_beta_matrix[,2],
  blfc_beta_matrix[,2],
  blfn_beta_matrix[,2],
  blsui_beta_matrix[,2]
)

colnames(baseline_SEs_only) <- c("blmd", "bltotasec", "blfm",
                                 "blfc", "blfn", "blsui")

baseline_ps_only <- cbind(
  blmd_beta_matrix[,3],
  bltotasec_beta_matrix[,3],
  blfm_beta_matrix[,3],
  blfc_beta_matrix[,3],
  blfn_beta_matrix[,3],
  blsui_beta_matrix[,3]
)

colnames(baseline_ps_only) <- c("blmd", "bltotasec", "blfm",
                                "blfc", "blfn", "blsui")

bl_lowest_ps <- matrix(nrow=30,ncol=6,dimnames=list(PRSs,
                                                    c("blmd", "bltotasec", "blfm",
                                                      "blfc", "blfn", "blsui")))

bl_lowest_ps[which(baseline_ps_only<=0.05)] <- 
  baseline_ps_only[which(baseline_ps_only<=0.05)]


heatmap.2( baseline_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(bl_lowest_ps,5),notecol="black",
           main="Baselines"
)

####

baseline_betas_only_2 <- cbind(
  blmd_beta_matrix_2[,1],
  bltotasec_beta_matrix_2[,1],
  blfm_beta_matrix_2[,1],
  blfc_beta_matrix_2[,1],
  blfn_beta_matrix_2[,1],
  blsui_beta_matrix_2[,1]
)

colnames(baseline_betas_only_2) <- c("blmd", "bltotasec", "blfm",
                                   "blfc", "blfn", "blsui")

baseline_SEs_only_2 <- cbind(
  blmd_beta_matrix_2[,2],
  bltotasec_beta_matrix_2[,2],
  blfm_beta_matrix_2[,2],
  blfc_beta_matrix_2[,2],
  blfn_beta_matrix_2[,2],
  blsui_beta_matrix_2[,2]
)

colnames(baseline_SEs_only_2) <- c("blmd", "bltotasec", "blfm",
                                 "blfc", "blfn", "blsui")

baseline_ps_only_2 <- cbind(
  blmd_beta_matrix_2[,3],
  bltotasec_beta_matrix_2[,3],
  blfm_beta_matrix_2[,3],
  blfc_beta_matrix_2[,3],
  blfn_beta_matrix_2[,3],
  blsui_beta_matrix_2[,3]
)

colnames(baseline_ps_only_2) <- c("blmd", "bltotasec", "blfm",
                                "blfc", "blfn", "blsui")

bl_lowest_ps_2 <- matrix(nrow=30,ncol=6,dimnames=list(PRSs,
                                                    c("blmd", "bltotasec", "blfm",
                                                      "blfc", "blfn", "blsui")))

bl_lowest_ps_2[which(baseline_ps_only_2<=0.05)] <- 
  baseline_ps_only_2[which(baseline_ps_only_2<=0.05)]


heatmap.2( baseline_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(bl_lowest_ps_2,5),notecol="black",
           main="Baselines, influential outliers removed"
)

##### 09) ANALYSE EFFECTS OF POLYGENIC SCORE ON DOSE #####

##### 09a) Build a model for dose #####

# I want to include a lot of covariates for theoretical reasons
# but inclusion of some, with missing data, causes loss of subjects.
# I am interested to see whether the baseline and time-varying
# depression and side effect scores are influential on dose.

length(levels(factor(GENDEP$subjectid))) # n = 647
length(levels(factor(subset(GENDEP,!is.na(zmadrs))$subjectid))) # n = 647
length(levels(factor(subset(GENDEP,!is.na(zblmd))$subjectid))) # n = 647
length(levels(factor(subset(GENDEP,!is.na(totasec))$subjectid))) # n = 638
length(levels(factor(subset(GENDEP,!is.na(bltotasec))$subjectid))) # n = 518

# Don't use bltotasec as a covariate

test <- subset(GENDEP, !is.na(zmadrs))
test <- subset(test, !is.na(zblmd))
test <- subset(test, !is.na(totasec))

# Compare models with and without time-varying depression
# and side-effect scores, and baseline depression

mod1 <- lmer( zdose ~ 
                SZ_0_0001 +
                drug +
                week + week2 + sex + cage + zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = test, REML=FALSE )

drop1(mod1) # Week2 is significant

mod2 <- lmer( zdose ~ 
                SZ_0_0001 +
                drug +
                zblmd + zmadrs + ztotasec + 
                week + week2 + sex + cage + zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = test, REML=FALSE )

drop1(mod2) # Keep zmadrs, ztotasec and zblmd

anova(mod1, mod2) # They make the model better, but it is still
# controversial to include depression score which could have a 
# reciprocal relationship on dose (a high dose might decrease 
# depression, but high depression might increase dose). If I 
# include as many covariates as possible, I will be more sure that
# the "signal" for the polygenic score is genuine - the polygenic
# score will be having an independent effect on dose somehow

mod2 # n = 572

hist(residuals(mod2)) 
residplot(mod2, newwd=FALSE)

##### 09b) Fit mixed effect models for dose #####

# Make lists to store models in

dose_models_both <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor <- vector(mode = "list", length = length(PRSs_31))

dose_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))

GENDEP_d <- subset(GENDEP, !is.na(zdose))
GENDEP_d <- subset(GENDEP_d, !is.na(zmadrs))
GENDEP_d <- subset(GENDEP_d, !is.na(totasec))
GENDEP_d <- subset(GENDEP_d, !is.na(zblmd))

GENDEP_d_e <- subset(GENDEP_d, drug==0)
GENDEP_d_n <- subset(GENDEP_d, drug==1)

# Fit models for dose as response variable

for ( i in 1:length( PRSs ) ){
  
  # First, fit the model for both drugs together
  
  dose_models_both[[i]] <- lmer( zdose ~
                                   GENDEP_d[,PRSs[i]] +
                                   drug +
                                   zblmd + zmadrs + ztotasec +
                                   week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                   (1|centreid) + (1|subjectid),
                                 data = GENDEP_d, REML=FALSE )
  
  removd <- romr.fnc(dose_models_both[[i]], GENDEP_d, trim = 2.5)$data # remove influential outliers
  
  dose_models_both_2[[i]] <- lmer( zdose ~
                                   removd[,PRSs[i]] +
                                   drug +
                                   zblmd + zmadrs + ztotasec +
                                   week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                   (1|centreid) + (1|subjectid),
                                 data = removd, REML=FALSE )
  
  # Then, fit the model for each drug separately
  
  dose_models_esc[[i]] <- lmer( zdose ~
                                  GENDEP_d_e[,PRSs[i]] +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = GENDEP_d_e, REML=FALSE )
  
  removd <- romr.fnc(dose_models_esc[[i]], GENDEP_d_e, trim = 2.5)$data # remove influential outliers
  
  dose_models_esc_2[[i]] <- lmer( zdose ~
                                  removd[,PRSs[i]] +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )

  dose_models_nor[[i]] <- lmer( zdose ~
                                  GENDEP_d_n[,PRSs[i]] +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = GENDEP_d_n, REML=FALSE )
  
  removd <- romr.fnc(dose_models_nor[[i]], GENDEP_d_n, trim = 2.5)$data # remove influential outliers
  
  dose_models_nor_2[[i]] <- lmer( zdose ~
                                  removd[,PRSs[i]] +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )
  
}

# Two models failed to converge (escit: HZ_0_0001, both(influential outliers removed): HZ_0_5)

# Add non-PRS models for comparison

dose_models_both[[31]] <- lmer( zdose ~
                                  drug +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = GENDEP_d, REML=FALSE )

removd <- romr.fnc(dose_models_both[[31]], GENDEP_d, trim = 2.5)$data # remove influential outliers

dose_models_both_2[[31]] <- lmer( zdose ~
                                  drug +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )

dose_models_esc[[31]] <- lmer( zdose ~
                                 zblmd + zmadrs + ztotasec +
                                 week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                 (1|centreid) + (1|subjectid),
                               data = GENDEP_d_e, REML=FALSE )

removd <- romr.fnc(dose_models_esc[[31]], GENDEP_d_e, trim = 2.5)$data # remove influential outliers

dose_models_esc_2[[31]] <- lmer( zdose ~
                                 zblmd + zmadrs + ztotasec +
                                 week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

dose_models_nor[[31]] <- lmer( zdose ~
                                 zblmd + zmadrs + ztotasec +
                                 week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                 (1|centreid) + (1|subjectid),
                               data = GENDEP_d_n, REML=FALSE )

removd <- romr.fnc(dose_models_nor[[31]], GENDEP_d_n, trim = 2.5)$data # remove influential outliers

dose_models_nor_2[[31]] <- lmer( zdose ~
                                 zblmd + zmadrs + ztotasec +
                                 week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

names(dose_models_both) <- PRSs_31
names(dose_models_esc) <- PRSs_31
names(dose_models_nor) <- PRSs_31

names(dose_models_both_2) <- PRSs_31
names(dose_models_esc_2) <- PRSs_31
names(dose_models_nor_2) <- PRSs_31

##### 09c) Get betas, standard errors and p-values for dose #####

dose_models_both_betas <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc_betas <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(dose_models_both)){
  dose_models_both_betas[i] <- list(summary(dose_models_both[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_esc)){
  dose_models_esc_betas[i] <- list(summary(dose_models_esc[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_nor)){
  dose_models_nor_betas[i] <- list(summary(dose_models_nor[[i]])$coefficients[2,])
}

names(dose_models_both_betas) <- PRSs_31
names(dose_models_esc_betas) <- PRSs_31
names(dose_models_nor_betas) <- PRSs_31

dose_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(dose_models_both_betas))[,1],
    t(as.data.frame(dose_models_both_betas))[,2],
    t(as.data.frame(dose_models_both_betas))[,5],
    t(as.data.frame(dose_models_esc_betas))[,1],
    t(as.data.frame(dose_models_esc_betas))[,2],
    t(as.data.frame(dose_models_esc_betas))[,5],
    t(as.data.frame(dose_models_nor_betas))[,1],
    t(as.data.frame(dose_models_nor_betas))[,2],
    t(as.data.frame(dose_models_nor_betas))[,5]
  )
)

colnames( dose_beta_matrix ) <- c(
  "dose_both_B", "dose_both_SE", "dose_both_p",
  "dose_esc_B", "dose_esc_SE", "dose_esc_p", 
  "dose_nor_B", "dose_nor_SE", "dose_nor_p"
)

####

dose_models_both_betas_2 <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc_betas_2 <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor_betas_2 <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(dose_models_both_2)){
  dose_models_both_betas_2[i] <- list(summary(dose_models_both_2[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_esc_2)){
  dose_models_esc_betas_2[i] <- list(summary(dose_models_esc_2[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_nor_2)){
  dose_models_nor_betas_2[i] <- list(summary(dose_models_nor_2[[i]])$coefficients[2,])
}

names(dose_models_both_betas_2) <- PRSs_31
names(dose_models_esc_betas_2) <- PRSs_31
names(dose_models_nor_betas_2) <- PRSs_31

dose_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(dose_models_both_betas_2))[,1],
    t(as.data.frame(dose_models_both_betas_2))[,2],
    t(as.data.frame(dose_models_both_betas_2))[,5],
    t(as.data.frame(dose_models_esc_betas_2))[,1],
    t(as.data.frame(dose_models_esc_betas_2))[,2],
    t(as.data.frame(dose_models_esc_betas_2))[,5],
    t(as.data.frame(dose_models_nor_betas_2))[,1],
    t(as.data.frame(dose_models_nor_betas_2))[,2],
    t(as.data.frame(dose_models_nor_betas_2))[,5]
  )
)

colnames( dose_beta_matrix_2 ) <- c(
  "dose_both_B", "dose_both_SE", "dose_both_p",
  "dose_esc_B", "dose_esc_SE", "dose_esc_p", 
  "dose_nor_B", "dose_nor_SE", "dose_nor_p"
)


##### 09d) Get PRSs' explained variances for dose #####

dose_matrix_both <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_esc <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_nor <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))

for (i in 1:31){
  dose_matrix_both[i,] <- r.squaredGLMM(dose_models_both[[i]])
  dose_matrix_esc[i,] <- r.squaredGLMM(dose_models_esc[[i]])
  dose_matrix_nor[i,] <- r.squaredGLMM(dose_models_nor[[i]])
}

####

dose_matrix_both_2 <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_esc_2 <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_nor_2 <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))

for (i in 1:31){
  dose_matrix_both_2[i,] <- r.squaredGLMM(dose_models_both_2[[i]])
  dose_matrix_esc_2[i,] <- r.squaredGLMM(dose_models_esc_2[[i]])
  dose_matrix_nor_2[i,] <- r.squaredGLMM(dose_models_nor_2[[i]])
}

##### 09e) Make output files for dose #####

sink( "dose_betas.csv")

cat("Dose betas\n,")
write.table( dose_beta_matrix[1:30,], col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_both, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_esc, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_nor, col.names=TRUE, sep="," )

sink()

####

sink( "dose_betas_inf_out_removed.csv")

cat("Dose betas\n,")
write.table( dose_beta_matrix_2[1:30,], col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_both_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_esc_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_nor_2, col.names=TRUE, sep="," )

sink()

##### 09f) Make heatmap of effect sizes for dose #####

dose_betas_only <- cbind(
  dose_beta_matrix[,1],
  dose_beta_matrix[,4],
  dose_beta_matrix[,7]
)

colnames(dose_betas_only) <- c("dose_both_B",
                                 "dose_esc_B",
                                 "dose_nor_B")

dose_betas_only <- dose_betas_only[1:30,]

dose_SEs_only <- cbind(
  dose_beta_matrix[,2],
  dose_beta_matrix[,5],
  dose_beta_matrix[,8]
)

colnames(dose_SEs_only) <- c("dose_both_SE",
                               "dose_esc_SE",
                               "dose_nor_SE")

dose_SEs_only <- dose_SEs_only[1:30,]

dose_ps_only <- cbind(
  dose_beta_matrix[,3],
  dose_beta_matrix[,6],
  dose_beta_matrix[,9]
)

colnames(dose_ps_only) <- c("dose_both_p",
                              "dose_esc_p",
                              "dose_nor_p")

dose_ps_only <- dose_ps_only[1:30,]

dose_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                      c("both", "esc", "nor")))

dose_lowest_ps[which(dose_ps_only<=0.05)] <- 
  dose_ps_only[which(dose_ps_only<=0.05)]


heatmap.2( dose_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(dose_lowest_ps,5),notecol="black",
           main="Dose"
)

####

dose_betas_only_2 <- cbind(
  dose_beta_matrix_2[,1],
  dose_beta_matrix_2[,4],
  dose_beta_matrix_2[,7]
)

colnames(dose_betas_only_2) <- c("dose_both_B",
                               "dose_esc_B",
                               "dose_nor_B")

dose_betas_only_2 <- dose_betas_only_2[1:30,]

dose_SEs_only_2 <- cbind(
  dose_beta_matrix_2[,2],
  dose_beta_matrix_2[,5],
  dose_beta_matrix_2[,8]
)

colnames(dose_SEs_only_2) <- c("dose_both_SE",
                             "dose_esc_SE",
                             "dose_nor_SE")

dose_SEs_only_2 <- dose_SEs_only_2[1:30,]

dose_ps_only_2 <- cbind(
  dose_beta_matrix_2[,3],
  dose_beta_matrix_2[,6],
  dose_beta_matrix_2[,9]
)

colnames(dose_ps_only_2) <- c("dose_both_p",
                            "dose_esc_p",
                            "dose_nor_p")

dose_ps_only_2 <- dose_ps_only_2[1:30,]

dose_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                      c("both", "esc", "nor")))

dose_lowest_ps_2[which(dose_ps_only_2<=0.05)] <- 
  dose_ps_only_2[which(dose_ps_only_2<=0.05)]


heatmap.2( dose_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(dose_lowest_ps_2,5),notecol="black",
           main="Dose, influential outliers removed"
)

##### 10) ANALYSE EFFECTS OF POLYGENIC SCORE ON DOSE: random only #####

##### 10a) Build a model for dose #####

# Randomized only!

test <- GEND_rand

# I want to include a lot of covariates for theoretical reasons
# but inclusion of some, with missing data, causes loss of subjects.
# I am interested to see whether the baseline and time-varying
# depression and side effect scores are influential on dose.

length(levels(factor(test$subjectid))) # n = 399
length(levels(factor(subset(test,!is.na(zmadrs))$subjectid))) # n = 399
length(levels(factor(subset(test,!is.na(zblmd))$subjectid))) # n = 399
length(levels(factor(subset(test,!is.na(totasec))$subjectid))) # n = 393
length(levels(factor(subset(test,!is.na(bltotasec))$subjectid))) # n = 349

# Don't use bltotasec as a covariate

test <- subset(GEND_rand, !is.na(zmadrs))
test <- subset(test, !is.na(zblmd))
test <- subset(test, !is.na(totasec))
test <- subset(test, !is.na(zdose))

# Compare models with and without time-varying depression
# and side-effect scores, and baseline depression

mod1 <- lmer( zdose ~ 
                SZ_0_0001 +
                drug +
                week + week2 + sex + cage + zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = test, REML=FALSE )

drop1(mod1) # Week2 is significant

mod2 <- lmer( zdose ~ 
                SZ_0_0001 +
                drug +
                zblmd + zmadrs + ztotasec + 
                week + week2 + sex + cage + zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = test, REML=FALSE )

drop1(mod2) # Keep zmadrs, ztotasec and zblmd, even though zblmd is not 
# significant, just to be the same as the "randomised and non-randomised"
# models

anova(mod1, mod2) # They make the model better, but it is still
# controversial to include depression score which could have a 
# reciprocal relationship on dose (a high dose might decrease 
# depression, but high depression might increase dose). If I 
# include as many covariates as possible, I will be more sure that
# the "signal" for the polygenic score is genuine - the polygenic
# score will be having an independent effect on dose somehow

mod2 # n = 370

hist(residuals(mod2)) 
predictmeans::residplot(mod2, newwd=FALSE)

##### 10b) Fit mixed effect models for dose #####

# Make lists to store models in

dose_models_both_r <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc_r <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor_r <- vector(mode = "list", length = length(PRSs_31))

dose_models_both_r_2 <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc_r_2 <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor_r_2 <- vector(mode = "list", length = length(PRSs_31))

GENDEP_d <- subset(GEND_rand, !is.na(zdose))
GENDEP_d <- subset(GENDEP_d, !is.na(zmadrs))
GENDEP_d <- subset(GENDEP_d, !is.na(totasec))
GENDEP_d <- subset(GENDEP_d, !is.na(zblmd))

GENDEP_d_e <- subset(GENDEP_d, drug==0)
GENDEP_d_n <- subset(GENDEP_d, drug==1)

# Fit models for dose as response variable

for ( i in 1:length( PRSs ) ){
  
  # First, fit the model for both drugs together
  
  dose_models_both_r[[i]] <- lmer( zdose ~
                                   GENDEP_d[,PRSs[i]] +
                                   drug +
                                   zblmd + zmadrs + ztotasec +
                                   week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                   (1|centreid) + (1|subjectid),
                                 data = GENDEP_d, REML=FALSE )
  
  removd <- romr.fnc(dose_models_both_r[[i]], GENDEP_d, trim = 2.5)$data # remove influential outliers
  
  dose_models_both_r_2[[i]] <- lmer( zdose ~
                                     removd[,PRSs[i]] +
                                     drug +
                                     zblmd + zmadrs + ztotasec +
                                     week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                     (1|centreid) + (1|subjectid),
                                   data = removd, REML=FALSE )
  
  # Then, fit the model for each drug separately
  
  dose_models_esc_r[[i]] <- lmer( zdose ~
                                  GENDEP_d_e[,PRSs[i]] +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = GENDEP_d_e, REML=FALSE )
  
  removd <- romr.fnc(dose_models_esc_r[[i]], GENDEP_d_e, trim = 2.5)$data # remove influential outliers
  
  dose_models_esc_r_2[[i]] <- lmer( zdose ~
                                    removd[,PRSs[i]] +
                                    zblmd + zmadrs + ztotasec +
                                    week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                    (1|centreid) + (1|subjectid),
                                  data = removd, REML=FALSE )
  
  dose_models_nor_r[[i]] <- lmer( zdose ~
                                  GENDEP_d_n[,PRSs[i]] +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = GENDEP_d_n, REML=FALSE )
  
  removd <- romr.fnc(dose_models_nor_r[[i]], GENDEP_d_n, trim = 2.5)$data # remove influential outliers
  
  dose_models_nor_r_2[[i]] <- lmer( zdose ~
                                    removd[,PRSs[i]] +
                                    zblmd + zmadrs + ztotasec +
                                    week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                    (1|centreid) + (1|subjectid),
                                  data = removd, REML=FALSE )
  
}

# Failure to converge (both: NZ_1, SZ_0_05, esc: AZ_0_5, nor(influential outliers
# removed): NZ_0_01, NZ_0_05)

# Add non-PRS models for comparison

dose_models_both_r[[31]] <- lmer( zdose ~
                                  drug +
                                  zblmd + zmadrs + ztotasec +
                                  week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                  (1|centreid) + (1|subjectid),
                                data = GENDEP_d, REML=FALSE )

removd <- romr.fnc(dose_models_both_r[[31]], GENDEP_d, trim = 2.5)$data # remove influential outliers

dose_models_both_r_2[[31]] <- lmer( zdose ~
                                    drug +
                                    zblmd + zmadrs + ztotasec +
                                    week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                    (1|centreid) + (1|subjectid),
                                  data = removd, REML=FALSE )

dose_models_esc_r[[31]] <- lmer( zdose ~
                                 zblmd + zmadrs + ztotasec +
                                 week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                 (1|centreid) + (1|subjectid),
                               data = GENDEP_d_e, REML=FALSE )

removd <- romr.fnc(dose_models_esc_r[[31]], GENDEP_d_e, trim = 2.5)$data # remove influential outliers

dose_models_esc_r_2[[31]] <- lmer( zdose ~
                                   zblmd + zmadrs + ztotasec +
                                   week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                   (1|centreid) + (1|subjectid),
                                 data = removd, REML=FALSE )

dose_models_nor_r[[31]] <- lmer( zdose ~
                                 zblmd + zmadrs + ztotasec +
                                 week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                 (1|centreid) + (1|subjectid),
                               data = GENDEP_d_n, REML=FALSE )

removd <- romr.fnc(dose_models_nor_r[[31]], GENDEP_d_n, trim = 2.5)$data # remove influential outliers

dose_models_nor_r_2[[31]] <- lmer( zdose ~
                                   zblmd + zmadrs + ztotasec +
                                   week + week2 + sex + cage +  zPC1 + zPC2 + zPC3 +
                                   (1|centreid) + (1|subjectid),
                                 data = removd, REML=FALSE )

names(dose_models_both_r) <- PRSs_31
names(dose_models_esc_r) <- PRSs_31
names(dose_models_nor_r) <- PRSs_31

names(dose_models_both_r_2) <- PRSs_31
names(dose_models_esc_r_2) <- PRSs_31
names(dose_models_nor_r_2) <- PRSs_31

##### 10c) Get betas, standard errors and p-values for dose #####

dose_models_both_r_betas <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc_r_betas <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor_r_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(dose_models_both_r)){
  dose_models_both_r_betas[i] <- list(summary(dose_models_both_r[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_esc_r)){
  dose_models_esc_r_betas[i] <- list(summary(dose_models_esc_r[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_nor_r)){
  dose_models_nor_r_betas[i] <- list(summary(dose_models_nor_r[[i]])$coefficients[2,])
}

names(dose_models_both_r_betas) <- PRSs_31
names(dose_models_esc_r_betas) <- PRSs_31
names(dose_models_nor_r_betas) <- PRSs_31

dose_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(dose_models_both_r_betas))[,1],
    t(as.data.frame(dose_models_both_r_betas))[,2],
    t(as.data.frame(dose_models_both_r_betas))[,5],
    t(as.data.frame(dose_models_esc_r_betas))[,1],
    t(as.data.frame(dose_models_esc_r_betas))[,2],
    t(as.data.frame(dose_models_esc_r_betas))[,5],
    t(as.data.frame(dose_models_nor_r_betas))[,1],
    t(as.data.frame(dose_models_nor_r_betas))[,2],
    t(as.data.frame(dose_models_nor_r_betas))[,5]
  )
)

colnames( dose_beta_matrix ) <- c(
  "dose_both_B", "dose_both_SE", "dose_both_p",
  "dose_esc_B", "dose_esc_SE", "dose_esc_p", 
  "dose_nor_B", "dose_nor_SE", "dose_nor_p"
)

####

dose_models_both_r_2_betas <- vector(mode = "list", length = length(PRSs_31))
dose_models_esc_r_2_betas <- vector(mode = "list", length = length(PRSs_31))
dose_models_nor_r_2_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(dose_models_both_r_2)){
  dose_models_both_r_2_betas[i] <- list(summary(dose_models_both_r_2[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_esc_r_2)){
  dose_models_esc_r_2_betas[i] <- list(summary(dose_models_esc_r_2[[i]])$coefficients[2,])
}
for (i in 1:length(dose_models_nor_r_2)){
  dose_models_nor_r_2_betas[i] <- list(summary(dose_models_nor_r_2[[i]])$coefficients[2,])
}

names(dose_models_both_r_2_betas) <- PRSs_31
names(dose_models_esc_r_2_betas) <- PRSs_31
names(dose_models_nor_r_2_betas) <- PRSs_31

dose_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(dose_models_both_r_2_betas))[,1],
    t(as.data.frame(dose_models_both_r_2_betas))[,2],
    t(as.data.frame(dose_models_both_r_2_betas))[,5],
    t(as.data.frame(dose_models_esc_r_2_betas))[,1],
    t(as.data.frame(dose_models_esc_r_2_betas))[,2],
    t(as.data.frame(dose_models_esc_r_2_betas))[,5],
    t(as.data.frame(dose_models_nor_r_2_betas))[,1],
    t(as.data.frame(dose_models_nor_r_2_betas))[,2],
    t(as.data.frame(dose_models_nor_r_2_betas))[,5]
  )
)

colnames( dose_beta_matrix_2 ) <- c(
  "dose_both_B", "dose_both_SE", "dose_both_p",
  "dose_esc_B", "dose_esc_SE", "dose_esc_p", 
  "dose_nor_B", "dose_nor_SE", "dose_nor_p"
)

##### 10d) Get PRSs' explained variances for dose #####

dose_matrix_both <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_esc <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_nor <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))

for (i in 1:31){
  dose_matrix_both[i,] <- r.squaredGLMM(dose_models_both_r[[i]])
  dose_matrix_esc[i,] <- r.squaredGLMM(dose_models_esc_r[[i]])
  dose_matrix_nor[i,] <- r.squaredGLMM(dose_models_nor_r[[i]])
}

####

dose_matrix_both_2 <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_esc_2 <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))
dose_matrix_nor_2 <- matrix(nrow=31,ncol=2,dimnames=list(PRSs_31,c("R2m","R2c")))

for (i in 1:31){
  dose_matrix_both_2[i,] <- r.squaredGLMM(dose_models_both_r_2[[i]])
  dose_matrix_esc_2[i,] <- r.squaredGLMM(dose_models_esc_r_2[[i]])
  dose_matrix_nor_2[i,] <- r.squaredGLMM(dose_models_nor_r_2[[i]])
}

##### 10e) Make output files for dose #####

sink( "dose_betas_rand_only.csv")

cat("Dose betas\n,")
write.table( dose_beta_matrix[1:30,], col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_both, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_esc, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_nor, col.names=TRUE, sep="," )

sink()

####

sink( "dose_betas_rand_only_inf_outliers_removed.csv")

cat("Dose betas\n,")
write.table( dose_beta_matrix_2[1:30,], col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_both_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_esc_2, col.names=TRUE, sep="," )
cat("\n,")
write.table( dose_matrix_nor_2, col.names=TRUE, sep="," )

sink()

##### 10f) Make heatmap of effect sizes for dose #####

dose_r_betas_only <- cbind(
  dose_beta_matrix[,1],
  dose_beta_matrix[,4],
  dose_beta_matrix[,7]
)

colnames(dose_r_betas_only) <- c("dose_both_B_r",
                                 "dose_esc_B_r",
                                 "dose_nor_B_r")

dose_r_betas_only <- dose_r_betas_only[1:30,]

dose_r_SEs_only <- cbind(
  dose_beta_matrix[,2],
  dose_beta_matrix[,5],
  dose_beta_matrix[,8]
)

colnames(dose_r_SEs_only) <- c("dose_both_SE_r",
                               "dose_esc_SE_r",
                               "dose_nor_SE_r")

dose_r_SEs_only <- dose_r_SEs_only[1:30,]

dose_r_ps_only <- cbind(
  dose_beta_matrix[,3],
  dose_beta_matrix[,6],
  dose_beta_matrix[,9]
)

colnames(dose_r_ps_only) <- c("dose_both_p_r",
                              "dose_esc_p_r",
                              "dose_nor_p_r")

dose_r_ps_only <- dose_r_ps_only[1:30,]

dose_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

dose_lowest_ps[which(dose_r_ps_only<=0.05)] <- 
  dose_r_ps_only[which(dose_r_ps_only<=0.05)]


heatmap.2( dose_r_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(dose_lowest_ps,5),notecol="black",
           main="Dose, randomized-participants-only"
)

####

dose_r_betas_only_2 <- cbind(
  dose_beta_matrix_2[,1],
  dose_beta_matrix_2[,4],
  dose_beta_matrix_2[,7]
)

colnames(dose_r_betas_only_2) <- c("dose_both_B_r",
                                 "dose_esc_B_r",
                                 "dose_nor_B_r")

dose_r_betas_only_2 <- dose_r_betas_only_2[1:30,]

dose_r_SEs_only_2 <- cbind(
  dose_beta_matrix_2[,2],
  dose_beta_matrix_2[,5],
  dose_beta_matrix_2[,8]
)

colnames(dose_r_SEs_only_2) <- c("dose_both_SE_r",
                               "dose_esc_SE_r",
                               "dose_nor_SE_r")

dose_r_SEs_only_2 <- dose_r_SEs_only_2[1:30,]

dose_r_ps_only_2 <- cbind(
  dose_beta_matrix_2[,3],
  dose_beta_matrix_2[,6],
  dose_beta_matrix_2[,9]
)

colnames(dose_r_ps_only_2) <- c("dose_both_p_r",
                              "dose_esc_p_r",
                              "dose_nor_p_r")

dose_r_ps_only_2 <- dose_r_ps_only_2[1:30,]

dose_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                      c("both", "esc", "nor")))

dose_lowest_ps_2[which(dose_r_ps_only_2<=0.05)] <- 
  dose_r_ps_only_2[which(dose_r_ps_only_2<=0.05)]


heatmap.2( dose_r_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(dose_lowest_ps_2,5),notecol="black",
           main="Dose, randomized-participants-only, influential outliers removed"
)

##### 11) ANALYSE EFFECTS OF POLYGENIC SCORE ON DEPRESSION (MADRS) #####

##### 11a) Build a model for ZMADRS #####

# Generally all variables are included for theoretical reasons
# and because they have been included in previous studies.
# PCs 1-3 were important in screeplot when QC-ing SNP data.
# I did test using dose as a covariate but ultimately I 
# decided a) depression score and dose influence each other
# and b) due to missing data, 66 people become excluded.
# Below code verifies that week2 is important.

mod1 <- lmer( zmadrs ~ 
                zblmd + week + 
                sex + cage + drug + AZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

mod2 <- lmer( zmadrs ~ 
                zblmd + week + week2 +
                sex + cage + drug + AZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

anova(mod1, mod2) # week2 improves the model

drop1(mod2) # week2 is significant

hist(residuals(mod2)) 
residplot(mod2, newwd=FALSE)

##### 11b) Fit mixed effects models for ZMADRS #####

# Make lists to store models in

z_models_both <- vector(mode = "list", length = length(PRSs_31))
z_models_esc <- vector(mode = "list", length = length(PRSs_31))
z_models_nor <- vector(mode = "list", length = length(PRSs_31))

z_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))

# Make new data set without missing data in important variables

GENDEP_z <- subset( GENDEP, !is.na(zmadrs) )
GENDEP_z_e <- subset( GENDEP_z, drug==0 )
GENDEP_z_n <- subset( GENDEP_z, drug==1 )

# Make a vector to note models which failed to converge

ftc_zmadrs <- vector(mode = "list", length = 0)
ftc_zmadrs_2 <- vector(mode = "list", length = 0)

# Fit models

for ( i in 1:length( PRSs ) ){
  
  z_mod_both <- lmer( zmadrs ~ zblmd + week + week2 +
                           sex + cage + drug + GENDEP_z[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_z, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_both)[-1]))==1){
    ftc_zmadrs[length(ftc_zmadrs)+1] <- z_mod_both
    names(ftc_zmadrs)[length(ftc_zmadrs)] <- PRSs[i]
  }
  z_models_both[[i]] <- z_mod_both
  
  removd <- romr.fnc(z_models_both[[i]], GENDEP_z, trim = 2.5)$data # remove influential outliers
  
  z_mod_both <- lmer( zmadrs ~ zblmd + week + week2 +
                        sex + cage + drug + removd[,PRSs[i]] +
                        zPC1 + zPC2 + zPC3 +
                        (1|centreid) + (1|subjectid),
                      data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_both)[-1]))==1){
    ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- z_mod_both
    names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- PRSs[i]
  }
  z_models_both_2[[i]] <- z_mod_both

  # Then, fit the models for each drug separately
  
  z_mod_esc <- lmer( zmadrs ~ zblmd + week + week2 +
                          sex + cage + GENDEP_z_e[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 +
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_z_e, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_esc)[-1]))==1){
    ftc_zmadrs[length(ftc_zmadrs)+1] <- z_mod_esc
    names(ftc_zmadrs)[length(ftc_zmadrs)] <- PRSs[i]
  }
  z_models_esc[[i]] <- z_mod_esc
  
  removd <- romr.fnc(z_models_esc[[i]], GENDEP_z_e, trim = 2.5)$data # remove influential outliers
  
  z_mod_esc <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + removd[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_esc)[-1]))==1){
    ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- z_mod_esc
    names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- PRSs[i]
  }
  z_models_esc_2[[i]] <- z_mod_esc
  
  z_mod_nor <- lmer( zmadrs ~ zblmd + week + week2 +
                          sex + cage + GENDEP_z_n[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 +
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_z_n, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_nor)[-1]))==1){
    ftc_zmadrs[length(ftc_zmadrs)+1] <- z_mod_nor
    names(ftc_zmadrs)[length(ftc_zmadrs)] <- PRSs[i]
  }
  z_models_nor[[i]] <- z_mod_nor
  
  removd <- romr.fnc(z_models_nor[[i]], GENDEP_z_n, trim = 2.5)$data # remove influential outliers
  
  z_mod_nor <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + removd[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_nor)[-1]))==1){
    ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- z_mod_nor
    names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- PRSs[i]
  }
  z_models_nor_2[[i]] <- z_mod_nor
  
}

names(z_models_both) <- PRSs_31
names(z_models_esc) <- PRSs_31
names(z_models_nor) <- PRSs_31

names(z_models_both_2) <- PRSs_31
names(z_models_esc_2) <- PRSs_31
names(z_models_nor_2) <- PRSs_31

# GENDEP_z, 5930, 647
# removd, 5810, 647
# GENDEP_z_e, 3517, 367
# removd, 3439, 367
# GENDEP_z_n, 2413, 280
# removd, 2368, 280

# One model ftc: nortriptyline, with influential outliers removed: NZ_0_01

# Look at zmadrs models without polygenic score

noPRS_both_z <- lmer( zmadrs ~ zblmd + week + week2 +
                        sex + cage + drug + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_z, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_z)[-1]))==1){
  ftc_zmadrs[length(ftc_zmadrs)+1] <- noPRS_both_z
  names(ftc_zmadrs)[length(ftc_zmadrs)] <- "No_PRS"
}
z_models_both[[31]] <- noPRS_both_z

removd <- romr.fnc(z_models_both[[31]], GENDEP_z, trim = 2.5)$data # remove influential outliers

noPRS_both_z <- lmer( zmadrs ~ zblmd + week + week2 +
                        sex + cage + drug + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_z)[-1]))==1){
  ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- noPRS_both_z
  names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- "No_PRS"
}
z_models_both_2[[31]] <- noPRS_both_z

####

noPRS_esc_z <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_z_e, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_z)[-1]))==1){
  ftc_zmadrs[length(ftc_zmadrs)+1] <- noPRS_esc_z
  names(ftc_zmadrs)[length(ftc_zmadrs)] <- "No_PRS"
}
z_models_esc[[31]] <- noPRS_esc_z

removd <- romr.fnc(z_models_esc[[31]], GENDEP_z_e, trim = 2.5)$data # remove influential outliers

noPRS_esc_z <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_z)[-1]))==1){
  ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- noPRS_esc_z
  names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- "No_PRS"
}
z_models_esc_2[[31]] <- noPRS_esc_z

####

noPRS_nor_z <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_z_n, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_z)[-1]))==1){
  ftc_zmadrs[length(ftc_zmadrs)+1] <- noPRS_nor_z
  names(ftc_zmadrs)[length(ftc_zmadrs)] <- "No_PRS"
}
z_models_nor[[31]] <- noPRS_nor_z

removd <- romr.fnc(z_models_nor[[31]], GENDEP_z_n, trim = 2.5)$data # remove influential outliers

noPRS_nor_z <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_z)[-1]))==1){
  ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- noPRS_nor_z
  names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- "No_PRS"
}
z_models_nor_2[[31]] <- noPRS_nor_z

##### 11c) Get betas, standard errors and p-values for zmadrs #####

z_models_both_betas <- vector(mode = "list", length = length(PRSs_31))
z_models_esc_betas <- vector(mode = "list", length = length(PRSs_31))
z_models_nor_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  z_models_both_betas[i] <- list(summary(z_models_both[[i]])$coefficients[8,])
}
for (i in 1:length(PRSs_31)){
  z_models_esc_betas[i] <- list(summary(z_models_esc[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  z_models_nor_betas[i] <- list(summary(z_models_nor[[i]])$coefficients[7,])
}

names(z_models_both_betas) <- PRSs_31
names(z_models_esc_betas) <- PRSs_31
names(z_models_nor_betas) <- PRSs_31

zmadrs_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(z_models_both_betas))[,1],
    t(as.data.frame(z_models_both_betas))[,2],
    t(as.data.frame(z_models_both_betas))[,5],
    t(as.data.frame(z_models_esc_betas))[,1],
    t(as.data.frame(z_models_esc_betas))[,2],
    t(as.data.frame(z_models_esc_betas))[,5],
    t(as.data.frame(z_models_nor_betas))[,1],
    t(as.data.frame(z_models_nor_betas))[,2],
    t(as.data.frame(z_models_nor_betas))[,5]
  )
)

colnames( zmadrs_beta_matrix ) <- c(
  "zmadrs_both_B", "zmadrs_both_SE","zmadrs_both_p",
  "zmadrs_esc_B", "zmadrs_esc_SE", "zmadrs_esc_p",
  "zmadrs_nor_B", "zmadrs_nor_SE", "zmadrs_nor_p"
)

####

z_models_both_betas_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_esc_betas_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_nor_betas_2 <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  z_models_both_betas_2[i] <- list(summary(z_models_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(PRSs_31)){
  z_models_esc_betas_2[i] <- list(summary(z_models_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  z_models_nor_betas_2[i] <- list(summary(z_models_nor_2[[i]])$coefficients[7,])
}

names(z_models_both_betas_2) <- PRSs_31
names(z_models_esc_betas_2) <- PRSs_31
names(z_models_nor_betas_2) <- PRSs_31

zmadrs_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(z_models_both_betas_2))[,1],
    t(as.data.frame(z_models_both_betas_2))[,2],
    t(as.data.frame(z_models_both_betas_2))[,5],
    t(as.data.frame(z_models_esc_betas_2))[,1],
    t(as.data.frame(z_models_esc_betas_2))[,2],
    t(as.data.frame(z_models_esc_betas_2))[,5],
    t(as.data.frame(z_models_nor_betas_2))[,1],
    t(as.data.frame(z_models_nor_betas_2))[,2],
    t(as.data.frame(z_models_nor_betas_2))[,5]
  )
)

colnames( zmadrs_beta_matrix_2 ) <- c(
  "zmadrs_both_B", "zmadrs_both_SE","zmadrs_both_p",
  "zmadrs_esc_B", "zmadrs_esc_SE", "zmadrs_esc_p",
  "zmadrs_nor_B", "zmadrs_nor_SE", "zmadrs_nor_p"
)

##### 11d) Get PRSs' explained variances for zmadrs #####

z_matrix_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

z_matrix_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

z_matrix_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  z_matrix_both[i,] <- r.squaredGLMM(z_models_both[[i]])
  z_matrix_esc[i,] <- r.squaredGLMM(z_models_esc[[i]])
  z_matrix_nor[i,] <- r.squaredGLMM(z_models_nor[[i]])
}


####

z_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

z_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

z_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  z_matrix_both_2[i,] <- r.squaredGLMM(z_models_both_2[[i]])
  z_matrix_esc_2[i,] <- r.squaredGLMM(z_models_esc_2[[i]])
  z_matrix_nor_2[i,] <- r.squaredGLMM(z_models_nor_2[[i]])
}

##### 11e) Make output files for zmadrs #####

sink( "zmadrs.csv")

# Betas and p-values

cat("Zmadrs_betas\n,")
write.table( zmadrs_beta_matrix[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(z_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(z_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(z_matrix_nor,
            col.names=TRUE, sep=",")

sink()

####

sink( "zmadrs_inf_outliers_removed.csv")

# Betas and p-values

cat("Zmadrs_betas\n,")
write.table( zmadrs_beta_matrix_2[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(z_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(z_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(z_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()

##### 11f) Make heatmap of effect sizes for zmadrs #####

z_betas_only <- cbind(
  zmadrs_beta_matrix[,1],
  zmadrs_beta_matrix[,4],
  zmadrs_beta_matrix[,7]
)

colnames(z_betas_only) <- c("z_both_B_r",
                                 "z_esc_B_r",
                                 "z_nor_B_r")

z_betas_only <- z_betas_only[1:30,]

z_SEs_only <- cbind(
  zmadrs_beta_matrix[,2],
  zmadrs_beta_matrix[,5],
  zmadrs_beta_matrix[,8]
)

colnames(z_SEs_only) <- c("z_both_SE_r",
                               "z_esc_SE_r",
                               "z_nor_SE_r")

z_SEs_only <- z_SEs_only[1:30,]

z_ps_only <- cbind(
  zmadrs_beta_matrix[,3],
  zmadrs_beta_matrix[,6],
  zmadrs_beta_matrix[,9]
)

colnames(z_ps_only) <- c("z_both_p_r",
                              "z_esc_p_r",
                              "z_nor_p_r")

z_ps_only <- z_ps_only[1:30,]

z_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                      c("both", "esc", "nor")))

z_lowest_ps[which(z_ps_only<=0.05)] <- 
  z_ps_only[which(z_ps_only<=0.05)]

heatmap.2( z_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(z_lowest_ps,5),notecol="black",
           main="MADRS"
)

####

z_betas_only_2 <- cbind(
  zmadrs_beta_matrix_2[,1],
  zmadrs_beta_matrix_2[,4],
  zmadrs_beta_matrix_2[,7]
)

colnames(z_betas_only_2) <- c("z_both_B_r",
                            "z_esc_B_r",
                            "z_nor_B_r")

z_betas_only_2 <- z_betas_only_2[1:30,]

z_SEs_only_2 <- cbind(
  zmadrs_beta_matrix_2[,2],
  zmadrs_beta_matrix_2[,5],
  zmadrs_beta_matrix_2[,8]
)

colnames(z_SEs_only_2) <- c("z_both_SE_r",
                          "z_esc_SE_r",
                          "z_nor_SE_r")

z_SEs_only_2 <- z_SEs_only_2[1:30,]

z_ps_only_2 <- cbind(
  zmadrs_beta_matrix_2[,3],
  zmadrs_beta_matrix_2[,6],
  zmadrs_beta_matrix_2[,9]
)

colnames(z_ps_only_2) <- c("z_both_p_r",
                         "z_esc_p_r",
                         "z_nor_p_r")

z_ps_only_2 <- z_ps_only_2[1:30,]

z_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                   c("both", "esc", "nor")))

z_lowest_ps_2[which(z_ps_only_2<=0.05)] <- 
  z_ps_only_2[which(z_ps_only_2<=0.05)]

heatmap.2( z_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(z_lowest_ps_2,5),notecol="black",
           main="MADRS, influential outliers removed"
)

##### 12) ANALYSE EFFECTS OF POLYGENIC SCORE ON DEPRESSION (MADRS): random only #####

##### 12a) Build a model for ZMADRS #####

# Generally all variables are included for theoretical reasons
# and because they have been included in previous studies.
# PCs 1-3 were important in screeplot when QC-ing SNP data.
# I did test using dose as a covariate but ultimately I 
# decided a) depression score and dose influence each other
# and b) due to missing data, 66 people become excluded.
# Below code verifies that week2 is important.

mod1 <- lmer( zmadrs ~ 
                zblmd + week + 
                sex + cage + drug + AZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

mod2 <- lmer( zmadrs ~ 
                zblmd + week + week2 +
                sex + cage + drug + AZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

anova(mod1, mod2) #  week2 improves the model

drop1(mod2) # week2 is significant

hist(residuals(mod2)) 
residplot(mod2, newwd=FALSE)

##### 12b) Fit mixed effects models for ZMADRS #####

# Make lists to store models in

z_models_both_r <- vector(mode = "list", length = length(PRSs_31))
z_models_esc_r <- vector(mode = "list", length = length(PRSs_31))
z_models_nor_r <- vector(mode = "list", length = length(PRSs_31))

z_models_both_r_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_esc_r_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_nor_r_2 <- vector(mode = "list", length = length(PRSs_31))

# Make new data set without missing data in important variables

GENDEP_z <- subset(GEND_rand, !is.na(zmadrs))

GENDEP_z_e <- subset(GENDEP_z, drug==0)
GENDEP_z_n <- subset(GENDEP_z, drug==1)

# Make a vector to note models which failed to converge

ftc_zmadrs <- vector(mode = "list", length = 0)
ftc_zmadrs_2 <- vector(mode = "list", length = 0)

# Fit models

for ( i in 1:length( PRSs ) ){
  
  z_mod_both <- lmer( zmadrs ~ zblmd + week + week2 +
                        sex + cage + drug + GENDEP_z[,PRSs[i]] +
                        zPC1 + zPC2 + zPC3 +
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_z, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_both)[-1]))==1){
    ftc_zmadrs[length(ftc_zmadrs)+1] <- z_mod_both
    names(ftc_zmadrs)[length(ftc_zmadrs)] <- PRSs[i]
  }
  z_models_both_r[[i]] <- z_mod_both
  
  removd <- romr.fnc(z_models_both_r[[i]], GENDEP_z, trim = 2.5)$data # remove influential outliers
  
  z_mod_both <- lmer( zmadrs ~ zblmd + week + week2 +
                        sex + cage + drug + removd[,PRSs[i]] +
                        zPC1 + zPC2 + zPC3 +
                        (1|centreid) + (1|subjectid),
                      data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_both)[-1]))==1){
    ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- z_mod_both
    names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- PRSs[i]
  }
  z_models_both_r_2[[i]] <- z_mod_both
  
  # Then, fit the models for each drug separately
  
  z_mod_esc <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + GENDEP_z_e[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_z_e, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_esc)[-1]))==1){
    ftc_zmadrs[length(ftc_zmadrs)+1] <- z_mod_esc
    names(ftc_zmadrs)[length(ftc_zmadrs)] <- PRSs[i]
  }
  z_models_esc_r[[i]] <- z_mod_esc
  
  removd <- romr.fnc(z_models_esc_r[[i]], GENDEP_z_e, trim = 2.5)$data # remove influential outliers
  
  z_mod_esc <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + removd[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_esc)[-1]))==1){
    ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- z_mod_esc
    names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- PRSs[i]
  }
  z_models_esc_r_2[[i]] <- z_mod_esc
  
  z_mod_nor <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + GENDEP_z_n[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_z_n, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_nor)[-1]))==1){
    ftc_zmadrs[length(ftc_zmadrs)+1] <- z_mod_nor
    names(ftc_zmadrs)[length(ftc_zmadrs)] <- PRSs[i]
  }
  z_models_nor_r[[i]] <- z_mod_nor
  
  removd <- romr.fnc(z_models_nor_r[[i]], GENDEP_z_n, trim = 2.5)$data # remove influential outliers
  
  z_mod_nor <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + removd[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(z_mod_nor)[-1]))==1){
    ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- z_mod_nor
    names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- PRSs[i]
  }
  z_models_nor_r_2[[i]] <- z_mod_nor
  
}

# GENDEP_z, 3626, 399
# removd, 3559, 399
# GENDEP_z_e, 1912, 198
# removd, 1873, 198
# GENDEP_z_n, 1714, 201
# removd, 1685, 201

# Look at zmadrs models without polygenic score

noPRS_both_r_z <- lmer( zmadrs ~ zblmd + week + week2 +
                        sex + cage + drug + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_z, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_r_z)[-1]))==1){
  ftc_zmadrs[length(ftc_zmadrs)+1] <- noPRS_both_r_z
  names(ftc_zmadrs)[length(ftc_zmadrs)] <- "No_PRS"
}
z_models_both_r[[31]] <- noPRS_both_r_z

removd <- romr.fnc(z_models_both_r[[31]], GENDEP_z, trim = 2.5)$data # remove influential outliers

noPRS_both_r_z <- lmer( zmadrs ~ zblmd + week + week2 +
                          sex + cage + drug + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_r_z)[-1]))==1){
  ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- noPRS_both_r_z
  names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- "No_PRS"
}
z_models_both_r_2[[31]] <- noPRS_both_r_z

####

noPRS_esc_r_z <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_z_e, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_r_z)[-1]))==1){
  ftc_zmadrs[length(ftc_zmadrs)+1] <- noPRS_esc_r_z
  names(ftc_zmadrs)[length(ftc_zmadrs)] <- "No_PRS"
}
z_models_esc_r[[31]] <- noPRS_esc_r_z

removd <- romr.fnc(z_models_esc_r[[31]], GENDEP_z_e, trim = 2.5)$data # remove influential outliers

noPRS_esc_r_z <- lmer( zmadrs ~ zblmd + week + week2 +
                         sex + cage + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_r_z)[-1]))==1){
  ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- noPRS_esc_r_z
  names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- "No_PRS"
}
z_models_esc_r_2[[31]] <- noPRS_esc_r_z

####

noPRS_nor_r_z <- lmer( zmadrs ~ zblmd + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_z_n, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_r_z)[-1]))==1){
  ftc_zmadrs[length(ftc_zmadrs)+1] <- noPRS_nor_r_z
  names(ftc_zmadrs)[length(ftc_zmadrs)] <- "No_PRS"
}
z_models_nor_r[[31]] <- noPRS_nor_r_z

removd <- romr.fnc(z_models_nor_r[[31]], GENDEP_z_n, trim = 2.5)$data # remove influential outliers

noPRS_nor_r_z <- lmer( zmadrs ~ zblmd + week + week2 +
                         sex + cage + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_r_z)[-1]))==1){
  ftc_zmadrs_2[length(ftc_zmadrs_2)+1] <- noPRS_nor_r_z
  names(ftc_zmadrs_2)[length(ftc_zmadrs_2)] <- "No_PRS"
}
z_models_nor_r_2[[31]] <- noPRS_nor_r_z

names(z_models_both_r) <- PRSs_31
names(z_models_esc_r) <- PRSs_31
names(z_models_nor_r) <- PRSs_31

names(z_models_both_r_2) <- PRSs_31
names(z_models_esc_r_2) <- PRSs_31
names(z_models_nor_r_2) <- PRSs_31

##### 12c) Get betas, standard errors and p-values for zmadrs #####

z_models_both_r_betas <- vector(mode = "list", length = length(PRSs_31))
z_models_esc_r_betas <- vector(mode = "list", length = length(PRSs_31))
z_models_nor_r_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  z_models_both_r_betas[i] <- list(summary(z_models_both_r[[i]])$coefficients[8,])
}
for (i in 1:length(PRSs_31)){
  z_models_esc_r_betas[i] <- list(summary(z_models_esc_r[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  z_models_nor_r_betas[i] <- list(summary(z_models_nor_r[[i]])$coefficients[7,])
}

names(z_models_both_r_betas) <- PRSs_31
names(z_models_esc_r_betas) <- PRSs_31
names(z_models_nor_r_betas) <- PRSs_31

zmadrs_beta_matrix_r <- as.matrix(
  cbind(
    t(as.data.frame(z_models_both_r_betas))[,1],
    t(as.data.frame(z_models_both_r_betas))[,2],
    t(as.data.frame(z_models_both_r_betas))[,5],
    t(as.data.frame(z_models_esc_r_betas))[,1],
    t(as.data.frame(z_models_esc_r_betas))[,2],
    t(as.data.frame(z_models_esc_r_betas))[,5],
    t(as.data.frame(z_models_nor_r_betas))[,1],
    t(as.data.frame(z_models_nor_r_betas))[,2],
    t(as.data.frame(z_models_nor_r_betas))[,5]
  )
)


colnames( zmadrs_beta_matrix_r ) <- c(
  "zmadrs_both_r_B", "zmadrs_both_r_SE","zmadrs_both_r_p",
  "zmadrs_esc_r_B", "zmadrs_esc_r_SE", "zmadrs_esc_r_p",
  "zmadrs_nor_r_B", "zmadrs_nor_r_SE", "zmadrs_nor_r_p"
)

####

z_models_both_r_betas_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_esc_r_betas_2 <- vector(mode = "list", length = length(PRSs_31))
z_models_nor_r_betas_2 <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  z_models_both_r_betas_2[i] <- list(summary(z_models_both_r_2[[i]])$coefficients[8,])
}
for (i in 1:length(PRSs_31)){
  z_models_esc_r_betas_2[i] <- list(summary(z_models_esc_r_2[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  z_models_nor_r_betas_2[i] <- list(summary(z_models_nor_r_2[[i]])$coefficients[7,])
}

names(z_models_both_r_betas_2) <- PRSs_31
names(z_models_esc_r_betas_2) <- PRSs_31
names(z_models_nor_r_betas_2) <- PRSs_31

zmadrs_beta_matrix_r_2 <- as.matrix(
  cbind(
    t(as.data.frame(z_models_both_r_betas_2))[,1],
    t(as.data.frame(z_models_both_r_betas_2))[,2],
    t(as.data.frame(z_models_both_r_betas_2))[,5],
    t(as.data.frame(z_models_esc_r_betas_2))[,1],
    t(as.data.frame(z_models_esc_r_betas_2))[,2],
    t(as.data.frame(z_models_esc_r_betas_2))[,5],
    t(as.data.frame(z_models_nor_r_betas_2))[,1],
    t(as.data.frame(z_models_nor_r_betas_2))[,2],
    t(as.data.frame(z_models_nor_r_betas_2))[,5]
  )
)

colnames( zmadrs_beta_matrix_r_2 ) <- c(
  "zmadrs_both_r_B", "zmadrs_both_r_SE","zmadrs_both_r_p",
  "zmadrs_esc_r_B", "zmadrs_esc_r_SE", "zmadrs_esc_r_p",
  "zmadrs_nor_r_B", "zmadrs_nor_r_SE", "zmadrs_nor_r_p"
)

##### 12d) Get PRSs' explained variances for zmadrs #####

z_matrix_both_r <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

z_matrix_esc_r <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

z_matrix_nor_r <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  z_matrix_both_r[i,] <- r.squaredGLMM(z_models_both_r[[i]])
  z_matrix_esc_r[i,] <- r.squaredGLMM(z_models_esc_r[[i]])
  z_matrix_nor_r[i,] <- r.squaredGLMM(z_models_nor_r[[i]])
}

####

z_matrix_both_r_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

z_matrix_esc_r_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

z_matrix_nor_r_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  z_matrix_both_r_2[i,] <- r.squaredGLMM(z_models_both_r_2[[i]])
  z_matrix_esc_r_2[i,] <- r.squaredGLMM(z_models_esc_r_2[[i]])
  z_matrix_nor_r_2[i,] <- r.squaredGLMM(z_models_nor_r_2[[i]])
}


##### 12e) Make output files for zmadrs #####

sink( "zmadrs_rand_only.csv")

# Betas and p-values

cat("Zmadrs_betas\n,")
write.table( zmadrs_beta_matrix_r[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(z_matrix_both_r,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(z_matrix_esc_r,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(z_matrix_nor_r,
            col.names=TRUE, sep=",")

sink()

####

sink( "zmadrs_rand_only_inf_outliers_removed.csv")

# Betas and p-values

cat("Zmadrs_betas\n,")
write.table( zmadrs_beta_matrix_r_2[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(z_matrix_both_r_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(z_matrix_esc_r_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(z_matrix_nor_r_2,
            col.names=TRUE, sep=",")

sink()

##### 12f) Make heatmap of effect sizes for zmadrs #####

z_r_betas_only <- cbind(
  zmadrs_beta_matrix_r[,1],
  zmadrs_beta_matrix_r[,4],
  zmadrs_beta_matrix_r[,7]
)

colnames(z_r_betas_only) <- c("z_both_r_B_r",
                            "z_esc_r_B_r",
                            "z_nor_r_B_r")

z_r_betas_only <- z_r_betas_only[1:30,]

z_r_SEs_only <- cbind(
  zmadrs_beta_matrix_r[,2],
  zmadrs_beta_matrix_r[,5],
  zmadrs_beta_matrix_r[,8]
)

colnames(z_r_SEs_only) <- c("z_both_r_SE_r",
                          "z_esc_r_SE_r",
                          "z_nor_r_SE_r")

z_r_SEs_only <- z_r_SEs_only[1:30,]

z_r_ps_only <- cbind(
  zmadrs_beta_matrix_r[,3],
  zmadrs_beta_matrix_r[,6],
  zmadrs_beta_matrix_r[,9]
)

colnames(z_r_ps_only) <- c("z_both_r_p_r",
                         "z_esc_r_p_r",
                         "z_nor_r_p_r")

z_r_ps_only <- z_r_ps_only[1:30,]

z_r_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                   c("both", "esc", "nor")))

z_r_lowest_ps[which(z_r_ps_only<=0.05)] <- 
  z_r_ps_only[which(z_r_ps_only<=0.05)]

heatmap.2( z_r_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(z_r_lowest_ps,5),notecol="black",
           main="MADRS, randomized-participants-only"
)

####

z_r_betas_only_2 <- cbind(
  zmadrs_beta_matrix_r_2[,1],
  zmadrs_beta_matrix_r_2[,4],
  zmadrs_beta_matrix_r_2[,7]
)

colnames(z_r_betas_only_2) <- c("z_both_r_B_r",
                              "z_esc_r_B_r",
                              "z_nor_r_B_r")

z_r_betas_only_2 <- z_r_betas_only_2[1:30,]

z_r_SEs_only_2 <- cbind(
  zmadrs_beta_matrix_r_2[,2],
  zmadrs_beta_matrix_r_2[,5],
  zmadrs_beta_matrix_r_2[,8]
)

colnames(z_r_SEs_only_2) <- c("z_both_r_SE_r",
                            "z_esc_r_SE_r",
                            "z_nor_r_SE_r")

z_r_SEs_only_2 <- z_r_SEs_only_2[1:30,]

z_r_ps_only_2 <- cbind(
  zmadrs_beta_matrix_r_2[,3],
  zmadrs_beta_matrix_r_2[,6],
  zmadrs_beta_matrix_r_2[,9]
)

colnames(z_r_ps_only_2) <- c("z_both_r_p_r",
                           "z_esc_r_p_r",
                           "z_nor_r_p_r")

z_r_ps_only_2 <- z_r_ps_only_2[1:30,]

z_r_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                     c("both", "esc", "nor")))

z_r_lowest_ps_2[which(z_r_ps_only_2<=0.05)] <- 
  z_r_ps_only_2[which(z_r_ps_only_2<=0.05)]

heatmap.2( z_r_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(z_r_lowest_ps_2,5),notecol="black",
           main="MADRS, randomized-participants-only, influential outliers removed"
)


##### 13) ANALYSE EFFECTS OF POLYGENIC SCORE ON SYMPTOM SUBSETS (F1-F3, suicidality) #####

# Dependent variables:

### zf1score: mood score
### zf2score: cognitive score
### zf3score: neurovegetative score
### suiscore: suicidality score

##### 13a) Build models for symptom subsets #####

# I know that I want to include many covariates on principle
# I'm interested to see whether the relevant symptom sub-sets'
# baselines, and week2, are significant

mod1 <- lmer( zf1score ~ blfm + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( zf2score ~ blfc + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( zf3score ~ blfn + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( suiscore ~ blsui + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

##### 13b) Fit mixed effects models for symptom subsets #####

# Make datasests with no missing required data

GENDEP_f1 <- subset( GENDEP, !is.na(zf1score) )
GENDEP_f1_e <- subset( GENDEP_f1, drug==0 )
GENDEP_f1_n <- subset( GENDEP_f1, drug==1 )

GENDEP_f2 <- subset( GENDEP, !is.na(zf2score) )
GENDEP_f2_e <- subset( GENDEP_f2, drug==0 )
GENDEP_f2_n <- subset( GENDEP_f2, drug==1 )

GENDEP_f3 <- subset( GENDEP, !is.na(zf3score) )
GENDEP_f3_e <- subset( GENDEP_f3, drug==0 )
GENDEP_f3_n <- subset( GENDEP_f3, drug==1 )

GENDEP_sui <- subset( GENDEP, !is.na(suiscore) )
GENDEP_sui_e <- subset( GENDEP_sui, drug==0 )
GENDEP_sui_n <- subset( GENDEP_sui, drug==1 )

# Make lists to store models in, for both drugs and separate drugs
# When fitting models, ensure REML=FALSE for assessing fixed effects

f1_models_both <- vector(mode = "list", length = length(PRSs_31))
f1_models_esc <- vector(mode = "list", length = length(PRSs_31))
f1_models_nor <- vector(mode = "list", length = length(PRSs_31))
f2_models_both <- vector(mode = "list", length = length(PRSs_31))
f2_models_esc <- vector(mode = "list", length = length(PRSs_31))
f2_models_nor <- vector(mode = "list", length = length(PRSs_31))
f3_models_both <- vector(mode = "list", length = length(PRSs_31))
f3_models_esc <- vector(mode = "list", length = length(PRSs_31))
f3_models_nor <- vector(mode = "list", length = length(PRSs_31))
sui_models_both <- vector(mode = "list", length = length(PRSs_31))
sui_models_esc <- vector(mode = "list", length = length(PRSs_31))
sui_models_nor <- vector(mode = "list", length = length(PRSs_31))

f1_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
f1_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
f1_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))
f2_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
f2_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
f2_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))
f3_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
f3_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
f3_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))
sui_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
sui_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
sui_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))

# Make a vector to note f1 models which failed to converge

ftc_zf1score <- vector(mode = "list", length = 0)
ftc_zf1score_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  f1_both <- lmer( zf1score ~ bl1fm + week + week2 +
                            sex + cage + drug + GENDEP_f1[,PRSs[i]] +
                            zPC1 + zPC2 + zPC3 +
                            (1|centreid) + (1|subjectid),
                          data = GENDEP_f1, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_both)[-1]))==1){
    ftc_zf1score[length(ftc_zf1score)+1] <- f1_both
    names(ftc_zf1score)[length(ftc_zf1score)] <- PRSs[i]
  }
  f1_models_both[[i]] <- f1_both
  
  removd <- romr.fnc(f1_models_both[[i]], GENDEP_f1, trim = 2.5)$data
  
  f1_both <- lmer( zf1score ~ bl1fm + week + week2 +
                     sex + cage + drug + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_both)[-1]))==1){
    ftc_zf1score_2[length(ftc_zf1score_2)+1] <- f1_both
    names(ftc_zf1score_2)[length(ftc_zf1score_2)] <- PRSs[i]
  }
  f1_models_both_2[[i]] <- f1_both  
  
  # Then, fit the models for each drug separately
  
  f1_esc <- lmer( zf1score ~ bl1fm + week + week2 +
                           sex + cage + GENDEP_f1_e[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f1_e, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_esc)[-1]))==1){
    ftc_zf1score[length(ftc_zf1score)+1] <- f1_esc
    names(ftc_zf1score)[length(ftc_zf1score)] <- PRSs[i]
  }
  f1_models_esc[[i]] <- f1_esc
  
  removd <- romr.fnc(f1_models_esc[[i]], GENDEP_f1_e, trim = 2.5)$data
  
  f1_esc <- lmer( zf1score ~ bl1fm + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_esc)[-1]))==1){
    ftc_zf1score_2[length(ftc_zf1score_2)+1] <- f1_esc
    names(ftc_zf1score_2)[length(ftc_zf1score_2)] <- PRSs[i]
  }
  f1_models_esc_2[[i]] <- f1_esc
 
  f1_nor <- lmer( zf1score ~ bl1fm + week + week2 +
                           sex + cage + GENDEP_f1_n[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f1_n, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_nor)[-1]))==1){
    ftc_zf1score[length(ftc_zf1score)+1] <- f1_nor
    names(ftc_zf1score)[length(ftc_zf1score)] <- PRSs[i]
  }
  f1_models_nor[[i]] <- f1_nor
  
  removd <- romr.fnc(f1_models_nor[[i]], GENDEP_f1_n, trim = 2.5)$data
  
  f1_nor <- lmer( zf1score ~ bl1fm + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_nor)[-1]))==1){
    ftc_zf1score_2[length(ftc_zf1score_2)+1] <- f1_nor
    names(ftc_zf1score_2)[length(ftc_zf1score_2)] <- PRSs[i]
  }
  f1_models_nor_2[[i]] <- f1_nor

}

# GENDEP_f1, 5996, 647
# removd, 5860, 647
# GENDEP_f1_e, 1912, 198
# removd, 1873, 198
# GENDEP_f1_n, 1714, 201
# removd, 1676, 201

# Make a vector to note f2 models which failed to converge

ftc_zf2score <- vector(mode = "list", length = 0)
ftc_zf2score_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  f2_both <- lmer( zf2score ~ bl1fc + week + week2 +
                            sex + cage + drug + GENDEP_f2[,PRSs[i]] +
                            zPC1 + zPC2 + zPC3 +
                            (1|centreid) + (1|subjectid),
                          data = GENDEP_f2, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_both)[-1]))==1){
    ftc_zf2score[length(ftc_zf2score)+1] <- f2_both
    names(ftc_zf2score)[length(ftc_zf2score)] <- PRSs[i]
  }
  f2_models_both[[i]] <- f2_both
  
  removd <- romr.fnc(f2_models_both[[i]], GENDEP_f2, trim = 2.5)$data
  
  f2_both <- lmer( zf2score ~ bl1fc + week + week2 +
                     sex + cage + drug + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_both)[-1]))==1){
    ftc_zf2score_2[length(ftc_zf2score_2)+1] <- f2_both
    names(ftc_zf2score_2)[length(ftc_zf2score_2)] <- PRSs[i]
  }
  f2_models_both_2[[i]] <- f2_both
  
  # Then, fit the models for each drug separately
  
  f2_esc <- lmer( zf2score ~ bl1fc + week + week2 +
                           sex + cage + GENDEP_f2_e[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f2_e, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_esc)[-1]))==1){
    ftc_zf2score[length(ftc_zf2score)+1] <- f2_esc
    names(ftc_zf2score)[length(ftc_zf2score)] <- PRSs[i]
  }
  f2_models_esc[[i]] <- f2_esc
  
  removd <- romr.fnc(f2_models_esc[[i]], GENDEP_f2_e, trim = 2.5)$data
  
  f2_esc <- lmer( zf2score ~ bl1fc + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_esc)[-1]))==1){
    ftc_zf2score_2[length(ftc_zf2score_2)+1] <- f2_esc
    names(ftc_zf2score_2)[length(ftc_zf2score_2)] <- PRSs[i]
  }
  f2_models_esc_2[[i]] <- f2_esc
  
  f2_nor <- lmer( zf2score ~ bl1fc + week + week2 +
                           sex + cage + GENDEP_f2_n[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f2_n, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_nor)[-1]))==1){
    ftc_zf2score[length(ftc_zf2score)+1] <- f2_nor
    names(ftc_zf2score)[length(ftc_zf2score)] <- PRSs[i]
  }
  f2_models_nor[[i]] <- f2_nor
  
  removd <- romr.fnc(f2_models_nor[[i]], GENDEP_f2_n, trim = 2.5)$data
  
  f2_nor <- lmer( zf2score ~ bl1fc + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_nor)[-1]))==1){
    ftc_zf2score_2[length(ftc_zf2score_2)+1] <- f2_nor
    names(ftc_zf2score_2)[length(ftc_zf2score_2)] <- PRSs[i]
  }
  f2_models_nor_2[[i]] <- f2_nor
  
}

# GENDEP_f2, 5996, 647
# removd, 5857, 647
# GENDEP_f2_e, 1912, 198
# removd, 1875, 198
# GENDEP_f2_n, 1714, 201
# removd, 1671, 201

# Make a vector to note f3 models which failed to converge

ftc_zf3score <- vector(mode = "list", length = 0)
ftc_zf3score_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  f3_both <- lmer( zf3score ~ bl1fn + week + week2 +
                            sex + cage + drug + GENDEP_f3[,PRSs[i]] +
                            zPC1 + zPC2 + zPC3 +
                            (1|centreid) + (1|subjectid),
                          data = GENDEP_f3, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_both)[-1]))==1){
    ftc_zf3score[length(ftc_zf3score)+1] <- f3_both
    names(ftc_zf3score)[length(ftc_zf3score)] <- PRSs[i]
  }
  f3_models_both[[i]] <- f3_both
  
  removd <- romr.fnc(f3_models_both[[i]], GENDEP_f3, trim = 2.5)$data
  
  f3_both <- lmer( zf3score ~ bl1fn + week + week2 +
                     sex + cage + drug + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_both)[-1]))==1){
    ftc_zf3score_2[length(ftc_zf3score_2)+1] <- f3_both
    names(ftc_zf3score_2)[length(ftc_zf3score_2)] <- PRSs[i]
  }
  f3_models_both_2[[i]] <- f3_both
  
  # Then, fit the models for each drug separately
  
  f3_esc <- lmer( zf3score ~ bl1fn + week + week2 +
                           sex + cage + GENDEP_f3_e[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f3_e, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_esc)[-1]))==1){
    ftc_zf3score[length(ftc_zf3score)+1] <- f3_esc
    names(ftc_zf3score)[length(ftc_zf3score)] <- PRSs[i]
  }
  f3_models_esc[[i]] <- f3_esc
  
  removd <- romr.fnc(f3_models_esc[[i]], GENDEP_f3_e, trim = 2.5)$data
  
  f3_esc <- lmer( zf3score ~ bl1fn + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_esc)[-1]))==1){
    ftc_zf3score_2[length(ftc_zf3score_2)+1] <- f3_esc
    names(ftc_zf3score_2)[length(ftc_zf3score_2)] <- PRSs[i]
  }
  f3_models_esc_2[[i]] <- f3_esc
  
  #######
  
  f3_nor<- lmer( zf3score ~ bl1fn + week + week2 +
                           sex + cage + GENDEP_f3_n[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f3_n, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_nor)[-1]))==1){
    ftc_zf3score[length(ftc_zf3score)+1] <- f3_nor
    names(ftc_zf3score)[length(ftc_zf3score)] <- PRSs[i]
  }
  f3_models_nor[[i]] <- f3_nor
  
  removd <- romr.fnc(f3_models_nor[[i]], GENDEP_f3_n, trim = 2.5)$data
  
  f3_nor<- lmer( zf3score ~ bl1fn + week + week2 +
                   sex + cage + removd[,PRSs[i]] +
                   zPC1 + zPC2 + zPC3 +
                   (1|centreid) + (1|subjectid),
                 data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_nor)[-1]))==1){
    ftc_zf3score_2[length(ftc_zf3score_2)+1] <- f3_nor
    names(ftc_zf3score_2)[length(ftc_zf3score_2)] <- PRSs[i]
  }
  f3_models_nor_2[[i]] <- f3_nor

}

# GENDEP_f3, 5996, 647
# removd, 5858, 647
# GENDEP_f3_e, 3551, 367
# removd, 3463, 367
# GENDEP_f3_n, 2445, 280
# removd, 2393, 280

# Make a vector to note sui models which failed to converge

ftc_suiscore <- vector(mode = "list", length = 0)
ftc_suiscore_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){

  sui_both <- lmer( suiscore ~ bl1sui + week + week2 +
                             sex + cage + drug + GENDEP_sui[,PRSs[i]] +
                             zPC1 + zPC2 + zPC3 +
                             (1|centreid) + (1|subjectid),
                           data = GENDEP_sui, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_both)[-1]))==1){
    ftc_suiscore[length(ftc_suiscore)+1] <- sui_both
    names(ftc_suiscore)[length(ftc_suiscore)] <- PRSs[i]
  }
  sui_models_both[[i]] <- sui_both
  
  removd <- romr.fnc(sui_models_both[[i]], GENDEP_sui, trim = 2.5)$data
  
  sui_both <- lmer( suiscore ~ bl1sui + week + week2 +
                      sex + cage + drug + removd[,PRSs[i]] +
                      zPC1 + zPC2 + zPC3 +
                      (1|centreid) + (1|subjectid),
                    data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_both)[-1]))==1){
    ftc_suiscore_2[length(ftc_suiscore_2)+1] <- sui_both
    names(ftc_suiscore_2)[length(ftc_suiscore_2)] <- PRSs[i]
  }
  sui_models_both_2[[i]] <- sui_both
  
  # Then, fit the models for each drug separately
  
  sui_esc <- lmer( suiscore ~ bl1sui + week + week2 +
                            sex + cage + GENDEP_sui_e[,PRSs[i]] +
                            zPC1 + zPC2 + zPC3 +
                            (1|centreid) + (1|subjectid),
                          data = GENDEP_sui_e, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_esc)[-1]))==1){
    ftc_suiscore[length(ftc_suiscore)+1] <- sui_esc
    names(ftc_suiscore)[length(ftc_suiscore)] <- PRSs[i]
  }
  sui_models_esc[[i]] <- sui_esc
  
  removd <- romr.fnc(sui_models_esc[[i]], GENDEP_sui_e, trim = 2.5)$data
  
  sui_esc <- lmer( suiscore ~ bl1sui + week + week2 +
                     sex + cage + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_esc)[-1]))==1){
    ftc_suiscore_2[length(ftc_suiscore_2)+1] <- sui_esc
    names(ftc_suiscore_2)[length(ftc_suiscore_2)] <- PRSs[i]
  }
  sui_models_esc_2[[i]] <- sui_esc
  
  
  #######
  
  sui_nor <- lmer( suiscore ~ bl1sui + week + week2 +
                            sex + cage + GENDEP_sui_n[,PRSs[i]] +
                            zPC1 + zPC2 + zPC3 +
                            (1|centreid) + (1|subjectid),
                          data = GENDEP_sui_n, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_nor)[-1]))==1){
    ftc_suiscore[length(ftc_suiscore)+1] <- sui_nor
    names(ftc_suiscore)[length(ftc_suiscore)] <- PRSs[i]
  }
  sui_models_nor[[i]] <- sui_nor
  
  removd <- romr.fnc(sui_models_nor[[i]], GENDEP_sui_n, trim = 2.5)$data
  
  sui_nor <- lmer( suiscore ~ bl1sui + week + week2 +
                     sex + cage + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_nor)[-1]))==1){
    ftc_suiscore_2[length(ftc_suiscore_2)+1] <- sui_nor
    names(ftc_suiscore_2)[length(ftc_suiscore_2)] <- PRSs[i]
  }
  sui_models_nor_2[[i]] <- sui_nor
  
  
}

# GENDEP_sui, 5996, 647
# removd, 5840, 647
# GENDEP_sui_e, 3551, 367
# removd, 3451, 367
# GENDEP_sui_n, 2445, 280
# removd, 2392, 280

# Models which failed to converge:

# f1score, escit-only, CZ_0_01
# f2score, both drugs, AZ_1
# f3score, escit-only, NZ_0_1
# f3score, both drugs, SZ_0_1
# suiscore, nortrip-only, NZ_0_01


# Look at models without polygenic score

f1_models_both[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                           sex + cage + drug + 
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f1, REML=FALSE )

removd <- romr.fnc(f1_models_both[[31]], GENDEP_f1, trim = 2.5)$data

f1_models_both_2[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                                sex + cage + drug + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = removd, REML=FALSE )

##

f1_models_esc[[31]] <-  lmer( zf1score ~ bl1fm + week + week2 +
                          sex + cage + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_f1_e, REML=FALSE )

removd <- romr.fnc(f1_models_esc[[31]], GENDEP_f1_e, trim = 2.5)$data

f1_models_esc_2[[31]] <-  lmer( zf1score ~ bl1fm + week + week2 +
                                sex + cage + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = removd, REML=FALSE )

##

f1_models_nor[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                          sex + cage + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_f1_n, REML=FALSE )

removd <- romr.fnc(f1_models_nor[[31]], GENDEP_f1_n, trim = 2.5)$data

f1_models_nor_2[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = removd, REML=FALSE )

##

f2_models_both[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                           sex + cage + drug + 
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f2, REML=FALSE )

removd <- romr.fnc(f2_models_both[[31]], GENDEP_f2, trim = 2.5)$data

f2_models_both_2[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                                sex + cage + drug + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = removd, REML=FALSE )

##

f2_models_esc[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                          sex + cage + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_f2_e, REML=FALSE )

removd <- romr.fnc(f2_models_esc[[31]], GENDEP_f2_e, trim = 2.5)$data

f2_models_esc_2[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = removd, REML=FALSE )

##

f2_models_nor[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                          sex + cage + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_f2_n, REML=FALSE )

removd <- romr.fnc(f2_models_nor[[31]], GENDEP_f2_n, trim = 2.5)$data

f2_models_nor_2[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = removd, REML=FALSE )

##

f3_models_both[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                           sex + cage + drug + 
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_f3, REML=FALSE )

removd <- romr.fnc(f3_models_both[[31]], GENDEP_f3, trim = 2.5)$data

f3_models_both_2[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                                sex + cage + drug + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = removd, REML=FALSE )

##

f3_models_esc[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                          sex + cage + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_f3_e, REML=FALSE )

removd <- romr.fnc(f3_models_esc[[31]], GENDEP_f3_e, trim = 2.5)$data

f3_models_esc_2[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = removd, REML=FALSE )

##

f3_models_nor[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                          sex + cage + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_f3_n, REML=FALSE )

removd <- romr.fnc(f3_models_nor[[31]], GENDEP_f3_n, trim = 2.5)$data

f3_models_nor_2[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = removd, REML=FALSE )

##

sui_models_both[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                            sex + cage + drug + 
                            zPC1 + zPC2 + zPC3 + 
                            (1|centreid) + (1|subjectid),
                          data = GENDEP_sui, REML=FALSE )

removd <- romr.fnc(sui_models_both[[31]], GENDEP_sui, trim = 2.5)$data

sui_models_both_2[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                 sex + cage + drug + 
                                 zPC1 + zPC2 + zPC3 + 
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

##

sui_models_esc[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                           sex + cage + 
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_sui_e, REML=FALSE )

removd <- romr.fnc(sui_models_esc[[31]], GENDEP_sui_e, trim = 2.5)$data

sui_models_esc_2[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                sex + cage + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = removd, REML=FALSE )

##

sui_models_nor[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                           sex + cage + 
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_sui_n, REML=FALSE )

removd <- romr.fnc(sui_models_nor[[31]], GENDEP_sui_n, trim = 2.5)$data

sui_models_nor_2[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                sex + cage + 
                                zPC1 + zPC2 + zPC3 +
                                (1|centreid) + (1|subjectid),
                              data = removd, REML=FALSE )

##

names(f1_models_both) <- PRSs_31
names(f1_models_esc) <- PRSs_31
names(f1_models_nor) <- PRSs_31

names(f2_models_both) <- PRSs_31
names(f2_models_esc) <- PRSs_31
names(f2_models_nor) <- PRSs_31

names(f3_models_both) <- PRSs_31
names(f3_models_esc) <- PRSs_31
names(f3_models_nor) <- PRSs_31

names(sui_models_both) <- PRSs_31
names(sui_models_esc) <- PRSs_31
names(sui_models_nor) <- PRSs_31

names(f1_models_both_2) <- PRSs_31
names(f1_models_esc_2) <- PRSs_31
names(f1_models_nor_2) <- PRSs_31

names(f2_models_both_2) <- PRSs_31
names(f2_models_esc_2) <- PRSs_31
names(f2_models_nor_2) <- PRSs_31

names(f3_models_both_2) <- PRSs_31
names(f3_models_esc_2) <- PRSs_31
names(f3_models_nor_2) <- PRSs_31

names(sui_models_both_2) <- PRSs_31
names(sui_models_esc_2) <- PRSs_31
names(sui_models_nor_2) <- PRSs_31

##### 13c) Get betas, standard errors and p-values for symptom subsets #####

f1_models_both_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_esc_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_nor_betas <- vector(mode = "list", length = length(PRSs_31))

f2_models_both_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_esc_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_nor_betas <- vector(mode = "list", length = length(PRSs_31))

f3_models_both_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_esc_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_nor_betas <- vector(mode = "list", length = length(PRSs_31))

sui_models_both_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_esc_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_nor_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(f1_models_both)){
  f1_models_both_betas[i] <- list(summary(f1_models_both[[i]])$coefficients[8,])
}
for (i in 1:length(f1_models_esc)){
  f1_models_esc_betas[i] <- list(summary(f1_models_esc[[i]])$coefficients[7,])
}
for (i in 1:length(f1_models_nor)){
  f1_models_nor_betas[i] <- list(summary(f1_models_nor[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_both)){
  f2_models_both_betas[i] <- list(summary(f2_models_both[[i]])$coefficients[8,])
}
for (i in 1:length(f2_models_esc)){
  f2_models_esc_betas[i] <- list(summary(f2_models_esc[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_nor)){
  f2_models_nor_betas[i] <- list(summary(f2_models_nor[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_both)){
  f3_models_both_betas[i] <- list(summary(f3_models_both[[i]])$coefficients[8,])
}
for (i in 1:length(f3_models_esc)){
  f3_models_esc_betas[i] <- list(summary(f3_models_esc[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_nor)){
  f3_models_nor_betas[i] <- list(summary(f3_models_nor[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_both)){
  sui_models_both_betas[i] <- list(summary(sui_models_both[[i]])$coefficients[8,])
}
for (i in 1:length(sui_models_esc)){
  sui_models_esc_betas[i] <- list(summary(sui_models_esc[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_nor)){
  sui_models_nor_betas[i] <- list(summary(sui_models_nor[[i]])$coefficients[7,])
}

names(f1_models_both_betas) <- PRSs_31
names(f1_models_esc_betas) <- PRSs_31
names(f1_models_nor_betas) <- PRSs_31
names(f2_models_both_betas) <- PRSs_31
names(f2_models_esc_betas) <- PRSs_31
names(f2_models_nor_betas) <- PRSs_31
names(f3_models_both_betas) <- PRSs_31
names(f3_models_esc_betas) <- PRSs_31
names(f3_models_nor_betas) <- PRSs_31
names(sui_models_both_betas) <- PRSs_31
names(sui_models_esc_betas) <- PRSs_31
names(sui_models_nor_betas) <- PRSs_31

f1score_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(f1_models_both_betas))[,1],
    t(as.data.frame(f1_models_both_betas))[,2],
    t(as.data.frame(f1_models_both_betas))[,5],
    t(as.data.frame(f1_models_esc_betas))[,1],
    t(as.data.frame(f1_models_esc_betas))[,2],
    t(as.data.frame(f1_models_esc_betas))[,5],
    t(as.data.frame(f1_models_nor_betas))[,1],
    t(as.data.frame(f1_models_nor_betas))[,2],
    t(as.data.frame(f1_models_nor_betas))[,5]
  )
)

f2score_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(f2_models_both_betas))[,1],
    t(as.data.frame(f2_models_both_betas))[,2],
    t(as.data.frame(f2_models_both_betas))[,5],
    t(as.data.frame(f2_models_esc_betas))[,1],
    t(as.data.frame(f2_models_esc_betas))[,2],
    t(as.data.frame(f2_models_esc_betas))[,5],
    t(as.data.frame(f2_models_nor_betas))[,1],
    t(as.data.frame(f2_models_nor_betas))[,2],
    t(as.data.frame(f2_models_nor_betas))[,5]
  )
)

f3score_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(f3_models_both_betas))[,1],
    t(as.data.frame(f3_models_both_betas))[,2],
    t(as.data.frame(f3_models_both_betas))[,5],
    t(as.data.frame(f3_models_esc_betas))[,1],
    t(as.data.frame(f3_models_esc_betas))[,2],
    t(as.data.frame(f3_models_esc_betas))[,5],
    t(as.data.frame(f3_models_nor_betas))[,1],
    t(as.data.frame(f3_models_nor_betas))[,2],
    t(as.data.frame(f3_models_nor_betas))[,5]
  )
)

suiscore_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(sui_models_both_betas))[,1],
    t(as.data.frame(sui_models_both_betas))[,2],
    t(as.data.frame(sui_models_both_betas))[,5],
    t(as.data.frame(sui_models_esc_betas))[,1],
    t(as.data.frame(sui_models_esc_betas))[,2],
    t(as.data.frame(sui_models_esc_betas))[,5],
    t(as.data.frame(sui_models_nor_betas))[,1],
    t(as.data.frame(sui_models_nor_betas))[,2],
    t(as.data.frame(sui_models_nor_betas))[,5]
  )
)

colnames( f1score_beta_matrix ) <- c(
  "f1score_both_B", "f1score_both_SE","f1score_both_p",
  "f1score_esc_B", "f1score_esc_SE", "f1score_esc_p",
  "f1score_nor_B", "f1score_nor_SE", "f1score_nor_p"
)

colnames( f2score_beta_matrix ) <- c(
  "f2score_both_B", "f2score_both_SE","f2score_both_p",
  "f2score_esc_B", "f2score_esc_SE", "f2score_esc_p",
  "f2score_nor_B", "f2score_nor_SE", "f2score_nor_p"
)

colnames( f3score_beta_matrix ) <- c(
  "f3score_both_B", "f3score_both_SE","f3score_both_p",
  "f3score_esc_B", "f3score_esc_SE", "f3score_esc_p",
  "f3score_nor_B", "f3score_nor_SE", "f3score_nor_p"
)

colnames( suiscore_beta_matrix ) <- c(
  "suiscore_both_B", "suiscore_both_SE","suiscore_both_p",
  "suiscore_esc_B", "suiscore_esc_SE", "suiscore_esc_p",
  "suiscore_nor_B", "suiscore_nor_SE", "suiscore_nor_p"
)


f1_models_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

f2_models_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

f3_models_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

sui_models_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(f1_models_both_2)){
  f1_models_both_2_betas[i] <- list(summary(f1_models_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(f1_models_esc_2)){
  f1_models_esc_2_betas[i] <- list(summary(f1_models_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(f1_models_nor_2)){
  f1_models_nor_2_betas[i] <- list(summary(f1_models_nor_2[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_both_2)){
  f2_models_both_2_betas[i] <- list(summary(f2_models_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(f2_models_esc_2)){
  f2_models_esc_2_betas[i] <- list(summary(f2_models_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_nor_2)){
  f2_models_nor_2_betas[i] <- list(summary(f2_models_nor_2[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_both_2)){
  f3_models_both_2_betas[i] <- list(summary(f3_models_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(f3_models_esc_2)){
  f3_models_esc_2_betas[i] <- list(summary(f3_models_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_nor_2)){
  f3_models_nor_2_betas[i] <- list(summary(f3_models_nor_2[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_both_2)){
  sui_models_both_2_betas[i] <- list(summary(sui_models_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(sui_models_esc_2)){
  sui_models_esc_2_betas[i] <- list(summary(sui_models_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_nor_2)){
  sui_models_nor_2_betas[i] <- list(summary(sui_models_nor_2[[i]])$coefficients[7,])
}

names(f1_models_both_2_betas) <- PRSs_31
names(f1_models_esc_2_betas) <- PRSs_31
names(f1_models_nor_2_betas) <- PRSs_31
names(f2_models_both_2_betas) <- PRSs_31
names(f2_models_esc_2_betas) <- PRSs_31
names(f2_models_nor_2_betas) <- PRSs_31
names(f3_models_both_2_betas) <- PRSs_31
names(f3_models_esc_2_betas) <- PRSs_31
names(f3_models_nor_2_betas) <- PRSs_31
names(sui_models_both_2_betas) <- PRSs_31
names(sui_models_esc_2_betas) <- PRSs_31
names(sui_models_nor_2_betas) <- PRSs_31

f1score_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(f1_models_both_2_betas))[,1],
    t(as.data.frame(f1_models_both_2_betas))[,2],
    t(as.data.frame(f1_models_both_2_betas))[,5],
    t(as.data.frame(f1_models_esc_2_betas))[,1],
    t(as.data.frame(f1_models_esc_2_betas))[,2],
    t(as.data.frame(f1_models_esc_2_betas))[,5],
    t(as.data.frame(f1_models_nor_2_betas))[,1],
    t(as.data.frame(f1_models_nor_2_betas))[,2],
    t(as.data.frame(f1_models_nor_2_betas))[,5]
  )
)

f2score_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(f2_models_both_2_betas))[,1],
    t(as.data.frame(f2_models_both_2_betas))[,2],
    t(as.data.frame(f2_models_both_2_betas))[,5],
    t(as.data.frame(f2_models_esc_2_betas))[,1],
    t(as.data.frame(f2_models_esc_2_betas))[,2],
    t(as.data.frame(f2_models_esc_2_betas))[,5],
    t(as.data.frame(f2_models_nor_2_betas))[,1],
    t(as.data.frame(f2_models_nor_2_betas))[,2],
    t(as.data.frame(f2_models_nor_2_betas))[,5]
  )
)

f3score_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(f3_models_both_2_betas))[,1],
    t(as.data.frame(f3_models_both_2_betas))[,2],
    t(as.data.frame(f3_models_both_2_betas))[,5],
    t(as.data.frame(f3_models_esc_2_betas))[,1],
    t(as.data.frame(f3_models_esc_2_betas))[,2],
    t(as.data.frame(f3_models_esc_2_betas))[,5],
    t(as.data.frame(f3_models_nor_2_betas))[,1],
    t(as.data.frame(f3_models_nor_2_betas))[,2],
    t(as.data.frame(f3_models_nor_2_betas))[,5]
  )
)

suiscore_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(sui_models_both_2_betas))[,1],
    t(as.data.frame(sui_models_both_2_betas))[,2],
    t(as.data.frame(sui_models_both_2_betas))[,5],
    t(as.data.frame(sui_models_esc_2_betas))[,1],
    t(as.data.frame(sui_models_esc_2_betas))[,2],
    t(as.data.frame(sui_models_esc_2_betas))[,5],
    t(as.data.frame(sui_models_nor_2_betas))[,1],
    t(as.data.frame(sui_models_nor_2_betas))[,2],
    t(as.data.frame(sui_models_nor_2_betas))[,5]
  )
)

colnames( f1score_beta_matrix_2 ) <- c(
  "f1score_both_B", "f1score_both_SE","f1score_both_p",
  "f1score_esc_B", "f1score_esc_SE", "f1score_esc_p",
  "f1score_nor_B", "f1score_nor_SE", "f1score_nor_p"
)

colnames( f2score_beta_matrix_2 ) <- c(
  "f2score_both_B", "f2score_both_SE","f2score_both_p",
  "f2score_esc_B", "f2score_esc_SE", "f2score_esc_p",
  "f2score_nor_B", "f2score_nor_SE", "f2score_nor_p"
)

colnames( f3score_beta_matrix_2 ) <- c(
  "f3score_both_B", "f3score_both_SE","f3score_both_p",
  "f3score_esc_B", "f3score_esc_SE", "f3score_esc_p",
  "f3score_nor_B", "f3score_nor_SE", "f3score_nor_p"
)

colnames( suiscore_beta_matrix_2 ) <- c(
  "suiscore_both_B", "suiscore_both_SE","suiscore_both_p",
  "suiscore_esc_B", "suiscore_esc_SE", "suiscore_esc_p",
  "suiscore_nor_B", "suiscore_nor_SE", "suiscore_nor_p"
)

##### 13d) Get PRSs' explained variances for symptom subsets #####

f1_matrix_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

for (i in 1:31){
  f1_matrix_both[i,] <- r.squaredGLMM(f1_models_both[[i]])
  f1_matrix_esc[i,] <- r.squaredGLMM(f1_models_esc[[i]])
  f1_matrix_nor[i,] <- r.squaredGLMM(f1_models_nor[[i]])
}

for (i in 1:31){
  f2_matrix_both[i,] <- r.squaredGLMM(f2_models_both[[i]])
  f2_matrix_esc[i,] <- r.squaredGLMM(f2_models_esc[[i]])
  f2_matrix_nor[i,] <- r.squaredGLMM(f2_models_nor[[i]])
}

for (i in 1:31){
  f3_matrix_both[i,] <- r.squaredGLMM(f3_models_both[[i]])
  f3_matrix_esc[i,] <- r.squaredGLMM(f3_models_esc[[i]])
  f3_matrix_nor[i,] <- r.squaredGLMM(f3_models_nor[[i]])
}

for (i in 1:31){
  sui_matrix_both[i,] <- r.squaredGLMM(sui_models_both[[i]])
  sui_matrix_esc[i,] <- r.squaredGLMM(sui_models_esc[[i]])
  sui_matrix_nor[i,] <- r.squaredGLMM(sui_models_nor[[i]])
}


f1_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

for (i in 1:31){
  f1_matrix_both_2[i,] <- r.squaredGLMM(f1_models_both_2[[i]])
  f1_matrix_esc_2[i,] <- r.squaredGLMM(f1_models_esc_2[[i]])
  f1_matrix_nor_2[i,] <- r.squaredGLMM(f1_models_nor_2[[i]])
}

for (i in 1:31){
  f2_matrix_both_2[i,] <- r.squaredGLMM(f2_models_both_2[[i]])
  f2_matrix_esc_2[i,] <- r.squaredGLMM(f2_models_esc_2[[i]])
  f2_matrix_nor_2[i,] <- r.squaredGLMM(f2_models_nor_2[[i]])
}

for (i in 1:31){
  f3_matrix_both_2[i,] <- r.squaredGLMM(f3_models_both_2[[i]])
  f3_matrix_esc_2[i,] <- r.squaredGLMM(f3_models_esc_2[[i]])
  f3_matrix_nor_2[i,] <- r.squaredGLMM(f3_models_nor_2[[i]])
}

for (i in 1:31){
  sui_matrix_both_2[i,] <- r.squaredGLMM(sui_models_both_2[[i]])
  sui_matrix_esc_2[i,] <- r.squaredGLMM(sui_models_esc_2[[i]])
  sui_matrix_nor_2[i,] <- r.squaredGLMM(sui_models_nor_2[[i]])
}

##### 13e) Make output files for symptom subsets #####

sink( "f1score_betas.csv")

# Betas and p-values

cat("f1score_betas\n,")
write.table( f1score_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(f1_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(f1_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(f1_matrix_nor,
            col.names=TRUE, sep=",")

sink()

sink( "f2score_betas.csv")

# Betas and p-values

cat("f2score_betas\n,")
write.table( f2score_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(f2_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(f2_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(f2_matrix_nor,
            col.names=TRUE, sep=",")

sink()

sink( "f3score_betas.csv")

# Betas and p-values

cat("f3score_betas\n,")
write.table( f3score_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(f3_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(f3_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(f3_matrix_nor,
            col.names=TRUE, sep=",")

sink()

sink( "suiscore_betas.csv")

# Betas and p-values

cat("suiscore_betas\n,")
write.table( suiscore_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(sui_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(sui_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(sui_matrix_nor,
            col.names=TRUE, sep=",")

sink()


sink( "f1score_betas_inf_out_removed.csv")

# Betas and p-values

cat("f1score_betas\n,")
write.table( f1score_beta_matrix_2, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(f1_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(f1_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(f1_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()

sink( "f2score_betas_inf_out_removed.csv")

# Betas and p-values

cat("f2score_betas\n,")
write.table( f2score_beta_matrix_2, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(f2_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(f2_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(f2_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()

sink( "f3score_betas_inf_out_removed.csv")

# Betas and p-values

cat("f3score_betas\n,")
write.table( f3score_beta_matrix_2, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(f3_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(f3_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(f3_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()

sink( "suiscore_betas_inf_out_removed.csv")

# Betas and p-values

cat("suiscore_betas\n,")
write.table( suiscore_beta_matrix_2, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both_REML\n,")
write.table(sui_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram_REML\n,")
write.table(sui_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline_REML\n,")
write.table(sui_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()


##### 13f) Make heatmaps of effect sizes for symptom subsets #####

f1score_betas_only <- cbind(f1score_beta_matrix[,1],
                            f1score_beta_matrix[,4],
                            f1score_beta_matrix[,7])

colnames(f1score_betas_only) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_betas_only <- f1score_betas_only[1:30,]

f1score_ses_only <- cbind(f1score_beta_matrix[,2],
                          f1score_beta_matrix[,5],
                          f1score_beta_matrix[,8])

colnames(f1score_ses_only) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_ses_only <- f1score_ses_only[1:30,]

f1score_p_only <- cbind(f1score_beta_matrix[,3],
                        f1score_beta_matrix[,6],
                        f1score_beta_matrix[,9])

colnames(f1score_p_only) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_p_only <- f1score_p_only[1:30,]

f1_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                      c("both", "esc", "nor")))

f1_lowest_ps[which(f1score_p_only<=0.05)] <- 
  f1score_p_only[which(f1score_p_only<=0.05)]

heatmap.2( f1score_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f1_lowest_ps,5),notecol="black",
           main="F1 (mood) score"
)

##

f1score_betas_only_2 <- cbind(f1score_beta_matrix_2[,1],
                            f1score_beta_matrix_2[,4],
                            f1score_beta_matrix_2[,7])

colnames(f1score_betas_only_2) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_betas_only_2 <- f1score_betas_only_2[1:30,]

f1score_ses_only_2 <- cbind(f1score_beta_matrix_2[,2],
                          f1score_beta_matrix_2[,5],
                          f1score_beta_matrix_2[,8])

colnames(f1score_ses_only_2) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_ses_only_2 <- f1score_ses_only_2[1:30,]

f1score_p_only_2 <- cbind(f1score_beta_matrix_2[,3],
                        f1score_beta_matrix_2[,6],
                        f1score_beta_matrix_2[,9])

colnames(f1score_p_only_2) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_p_only_2 <- f1score_p_only_2[1:30,]

f1_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f1_lowest_ps_2[which(f1score_p_only_2<=0.05)] <- 
  f1score_p_only_2[which(f1score_p_only_2<=0.05)]

heatmap.2( f1score_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f1_lowest_ps_2,5),notecol="black",
           main="F1 (mood) score, inf out rem"
)

##

f2score_betas_only <- cbind(f2score_beta_matrix[,1],
                            f2score_beta_matrix[,4],
                            f2score_beta_matrix[,7])

colnames(f2score_betas_only) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_betas_only <- f2score_betas_only[1:30,]

f2score_ses_only <- cbind(f2score_beta_matrix[,2],
                          f2score_beta_matrix[,5],
                          f2score_beta_matrix[,8])

colnames(f2score_ses_only) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_ses_only <- f2score_ses_only[1:30,]

f2score_p_only <- cbind(f2score_beta_matrix[,3],
                        f2score_beta_matrix[,6],
                        f2score_beta_matrix[,9])

colnames(f2score_p_only) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_p_only <- f2score_p_only[1:30,]

f2_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f2_lowest_ps[which(f2score_p_only<=0.05)] <- 
  f2score_p_only[which(f2score_p_only<=0.05)]

heatmap.2( f2score_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f2_lowest_ps,5),notecol="black",
           main="F2 (cognitive) score"
)

##

f2score_betas_only_2 <- cbind(f2score_beta_matrix_2[,1],
                            f2score_beta_matrix_2[,4],
                            f2score_beta_matrix_2[,7])

colnames(f2score_betas_only_2) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_betas_only_2 <- f2score_betas_only_2[1:30,]

f2score_ses_only_2 <- cbind(f2score_beta_matrix_2[,2],
                          f2score_beta_matrix_2[,5],
                          f2score_beta_matrix_2[,8])

colnames(f2score_ses_only_2) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_ses_only_2 <- f2score_ses_only_2[1:30,]

f2score_p_only_2 <- cbind(f2score_beta_matrix_2[,3],
                        f2score_beta_matrix_2[,6],
                        f2score_beta_matrix_2[,9])

colnames(f2score_p_only_2) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_p_only_2 <- f2score_p_only_2[1:30,]

f2_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f2_lowest_ps_2[which(f2score_p_only_2<=0.05)] <- 
  f2score_p_only_2[which(f2score_p_only_2<=0.05)]

heatmap.2( f2score_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f2_lowest_ps_2,5),notecol="black",
           main="F2 (cognitive) score, inf out rem"
)


##

f3score_betas_only <- cbind(f3score_beta_matrix[,1],
                            f3score_beta_matrix[,4],
                            f3score_beta_matrix[,7])

colnames(f3score_betas_only) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_betas_only <- f3score_betas_only[1:30,]

f3score_ses_only <- cbind(f3score_beta_matrix[,2],
                          f3score_beta_matrix[,5],
                          f3score_beta_matrix[,8])

colnames(f3score_ses_only) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_ses_only <- f3score_ses_only[1:30,]

f3score_p_only <- cbind(f3score_beta_matrix[,3],
                        f3score_beta_matrix[,6],
                        f3score_beta_matrix[,9])

colnames(f3score_p_only) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_p_only <- f3score_p_only[1:30,]

f3_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f3_lowest_ps[which(f3score_p_only<=0.05)] <- 
  f3score_p_only[which(f3score_p_only<=0.05)]

heatmap.2( f3score_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f3_lowest_ps,5),notecol="black",
           main="F3 (neurovegetative) score"
)

##

f3score_betas_only_2 <- cbind(f3score_beta_matrix_2[,1],
                            f3score_beta_matrix_2[,4],
                            f3score_beta_matrix_2[,7])

colnames(f3score_betas_only_2) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_betas_only_2 <- f3score_betas_only_2[1:30,]

f3score_ses_only_2 <- cbind(f3score_beta_matrix_2[,2],
                          f3score_beta_matrix_2[,5],
                          f3score_beta_matrix_2[,8])

colnames(f3score_ses_only_2) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_ses_only_2 <- f3score_ses_only_2[1:30,]

f3score_p_only_2 <- cbind(f3score_beta_matrix_2[,3],
                        f3score_beta_matrix_2[,6],
                        f3score_beta_matrix_2[,9])

colnames(f3score_p_only_2) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_p_only_2 <- f3score_p_only_2[1:30,]

f3_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f3_lowest_ps_2[which(f3score_p_only_2<=0.05)] <- 
  f3score_p_only_2[which(f3score_p_only_2<=0.05)]

heatmap.2( f3score_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f3_lowest_ps_2,5),notecol="black",
           main="F3 (neurovegetative) score, inf out rem"
)


##

suiscore_betas_only <- cbind(suiscore_beta_matrix[,1],
                            suiscore_beta_matrix[,4],
                            suiscore_beta_matrix[,7])

colnames(suiscore_betas_only) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_betas_only <- suiscore_betas_only[1:30,]


suiscore_ses_only <- cbind(suiscore_beta_matrix[,2],
                          suiscore_beta_matrix[,5],
                          suiscore_beta_matrix[,8])

colnames(suiscore_ses_only) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_ses_only <- suiscore_ses_only[1:30,]

suiscore_p_only <- cbind(suiscore_beta_matrix[,3],
                        suiscore_beta_matrix[,6],
                        suiscore_beta_matrix[,9])

colnames(suiscore_p_only) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_p_only <- suiscore_p_only[1:30,]

sui_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

sui_lowest_ps[which(suiscore_p_only<=0.05)] <- 
  suiscore_p_only[which(suiscore_p_only<=0.05)]

heatmap.2( suiscore_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(sui_lowest_ps,5),notecol="black",
           main="Suicidality score"
)

##

suiscore_betas_only_2 <- cbind(suiscore_beta_matrix_2[,1],
                             suiscore_beta_matrix_2[,4],
                             suiscore_beta_matrix_2[,7])

colnames(suiscore_betas_only_2) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_betas_only_2 <- suiscore_betas_only_2[1:30,]


suiscore_ses_only_2 <- cbind(suiscore_beta_matrix_2[,2],
                           suiscore_beta_matrix_2[,5],
                           suiscore_beta_matrix_2[,8])

colnames(suiscore_ses_only_2) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_ses_only_2 <- suiscore_ses_only_2[1:30,]

suiscore_p_only_2 <- cbind(suiscore_beta_matrix_2[,3],
                         suiscore_beta_matrix_2[,6],
                         suiscore_beta_matrix_2[,9])

colnames(suiscore_p_only_2) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_p_only_2 <- suiscore_p_only_2[1:30,]

sui_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                     c("both", "esc", "nor")))

sui_lowest_ps_2[which(suiscore_p_only_2<=0.05)] <- 
  suiscore_p_only_2[which(suiscore_p_only_2<=0.05)]

heatmap.2( suiscore_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(sui_lowest_ps_2,5),notecol="black",
           main="Suicidality score, inf out rem"
)

##### 14) ANALYSE EFFECTS OF POLYGENIC SCORE ON SYMPTOM SUBSETS (F1-F3, suicidality): random only #####

# Dependent variables:

### zf1score: mood score
### zf2score: cognitive score
### zf3score: neurovegetative score
### suiscore: suicidality score

##### 14a) Build models for symptom subsets #####

# I know that I want to include many covariates on principle
# I'm interested to see whether the relevant symptom sub-sets'
# baselines, and week2, are significant

mod1 <- lmer( zf1score ~ blfm + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( zf2score ~ blfc + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( zf3score ~ blfn + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

mod1 <- lmer( suiscore ~ blsui + week + week2 +
                sex + cage + drug + NZ_1 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

drop1(mod1) # They are

hist(residuals(mod1)) 
residplot(mod1, newwd=FALSE)

# I subsequently found singular fit issues with the suiscore model in 
# nortriptyline only, randomized-only...

drop1(
  lmer( suiscore ~ blsui + week + week2 +
          sex + cage + NZ_1 +
          zPC1 + zPC2 + zPC3 +
          (1|centreid) + (1|subjectid),
        data = GEND_rand_n, REML=FALSE )
)

drop1(
  lmer( suiscore ~ blsui + week + week2 +
          sex + cage + NZ_1 +
          zPC1 + zPC2 + zPC3 +
           (1|subjectid),
        data = GEND_rand_n, REML=FALSE )
) # Removing centre helps

##### 14b) Fit mixed effects models for symptom subsets #####

# Make datasests with no missing required data

GENDEP_f1 <- subset( GEND_rand, !is.na(zf1score) )
GENDEP_f1_e <- subset( GENDEP_f1, drug==0 )
GENDEP_f1_n <- subset( GENDEP_f1, drug==1 )

GENDEP_f2 <- subset( GEND_rand, !is.na(zf2score) )
GENDEP_f2_e <- subset( GENDEP_f2, drug==0 )
GENDEP_f2_n <- subset( GENDEP_f2, drug==1 )

GENDEP_f3 <- subset( GEND_rand, !is.na(zf3score) )
GENDEP_f3_e <- subset( GENDEP_f3, drug==0 )
GENDEP_f3_n <- subset( GENDEP_f3, drug==1 )

GENDEP_sui <- subset( GEND_rand, !is.na(suiscore) )
GENDEP_sui_e <- subset( GENDEP_sui, drug==0 )
GENDEP_sui_n <- subset( GENDEP_sui, drug==1 )

# Make lists to store models in, for both drugs and separate drugs
# When fitting models, ensure REML=FALSE for assessing fixed effects

f1_models_r_both <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_esc <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_nor <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_both <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_esc <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_nor <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_both <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_esc <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_nor <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_both <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_esc <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_nor <- vector(mode = "list", length = length(PRSs_31))

f1_models_r_both_2 <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_esc_2 <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_nor_2 <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_both_2 <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_esc_2 <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_nor_2 <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_both_2 <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_esc_2 <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_nor_2 <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_both_2 <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_esc_2 <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_nor_2 <- vector(mode = "list", length = length(PRSs_31))

# Make a vector to note f1 models_r which failed to converge

ftc_zf1score <- vector(mode = "list", length = 0)
ftc_zf1score_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  f1_both <- lmer( zf1score ~ bl1fm + week + week2 +
                     sex + cage + drug + GENDEP_f1[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = GENDEP_f1, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_both)[-1]))==1){
    ftc_zf1score[length(ftc_zf1score)+1] <- f1_both
    names(ftc_zf1score)[length(ftc_zf1score)] <- PRSs[i]
  }
  f1_models_r_both[[i]] <- f1_both
  
  removd <- romr.fnc(f1_models_r_both[[i]], GENDEP_f1, trim = 2.5)$data
  
  f1_both <- lmer( zf1score ~ bl1fm + week + week2 +
                     sex + cage + drug + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_both)[-1]))==1){
    ftc_zf1score_2[length(ftc_zf1score_2)+1] <- f1_both
    names(ftc_zf1score_2)[length(ftc_zf1score_2)] <- PRSs[i]
  }
  f1_models_r_both_2[[i]] <- f1_both  
  
  # Then, fit the models_r for each drug separately
  
  f1_esc <- lmer( zf1score ~ bl1fm + week + week2 +
                    sex + cage + GENDEP_f1_e[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = GENDEP_f1_e, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_esc)[-1]))==1){
    ftc_zf1score[length(ftc_zf1score)+1] <- f1_esc
    names(ftc_zf1score)[length(ftc_zf1score)] <- PRSs[i]
  }
  f1_models_r_esc[[i]] <- f1_esc
  
  removd <- romr.fnc(f1_models_r_esc[[i]], GENDEP_f1_e, trim = 2.5)$data
  
  f1_esc <- lmer( zf1score ~ bl1fm + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_esc)[-1]))==1){
    ftc_zf1score_2[length(ftc_zf1score_2)+1] <- f1_esc
    names(ftc_zf1score_2)[length(ftc_zf1score_2)] <- PRSs[i]
  }
  f1_models_r_esc_2[[i]] <- f1_esc
  
  f1_nor <- lmer( zf1score ~ bl1fm + week + week2 +
                    sex + cage + GENDEP_f1_n[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = GENDEP_f1_n, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_nor)[-1]))==1){
    ftc_zf1score[length(ftc_zf1score)+1] <- f1_nor
    names(ftc_zf1score)[length(ftc_zf1score)] <- PRSs[i]
  }
  f1_models_r_nor[[i]] <- f1_nor
  
  removd <- romr.fnc(f1_models_r_nor[[i]], GENDEP_f1_n, trim = 2.5)$data
  
  f1_nor <- lmer( zf1score ~ bl1fm + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f1_nor)[-1]))==1){
    ftc_zf1score_2[length(ftc_zf1score_2)+1] <- f1_nor
    names(ftc_zf1score_2)[length(ftc_zf1score_2)] <- PRSs[i]
  }
  f1_models_r_nor_2[[i]] <- f1_nor
  
}

# Make a vector to note f2 models_r which failed to converge


ftc_zf2score <- vector(mode = "list", length = 0)
ftc_zf2score_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  f2_both <- lmer( zf2score ~ bl1fc + week + week2 +
                     sex + cage + drug + GENDEP_f2[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = GENDEP_f2, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_both)[-1]))==1){
    ftc_zf2score[length(ftc_zf2score)+1] <- f2_both
    names(ftc_zf2score)[length(ftc_zf2score)] <- PRSs[i]
  }
  f2_models_r_both[[i]] <- f2_both
  
  removd <- romr.fnc(f2_models_r_both[[i]], GENDEP_f2, trim = 2.5)$data
  
  f2_both <- lmer( zf2score ~ bl1fc + week + week2 +
                     sex + cage + drug + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_both)[-1]))==1){
    ftc_zf2score_2[length(ftc_zf2score_2)+1] <- f2_both
    names(ftc_zf2score_2)[length(ftc_zf2score_2)] <- PRSs[i]
  }
  f2_models_r_both_2[[i]] <- f2_both
  
  # Then, fit the models_r for each drug separately
  
  f2_esc <- lmer( zf2score ~ bl1fc + week + week2 +
                    sex + cage + GENDEP_f2_e[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = GENDEP_f2_e, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_esc)[-1]))==1){
    ftc_zf2score[length(ftc_zf2score)+1] <- f2_esc
    names(ftc_zf2score)[length(ftc_zf2score)] <- PRSs[i]
  }
  f2_models_r_esc[[i]] <- f2_esc
  
  removd <- romr.fnc(f2_models_r_esc[[i]], GENDEP_f2_e, trim = 2.5)$data
  
  f2_esc <- lmer( zf2score ~ bl1fc + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_esc)[-1]))==1){
    ftc_zf2score_2[length(ftc_zf2score_2)+1] <- f2_esc
    names(ftc_zf2score_2)[length(ftc_zf2score_2)] <- PRSs[i]
  }
  f2_models_r_esc_2[[i]] <- f2_esc
  
  f2_nor <- lmer( zf2score ~ bl1fc + week + week2 +
                    sex + cage + GENDEP_f2_n[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = GENDEP_f2_n, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_nor)[-1]))==1){
    ftc_zf2score[length(ftc_zf2score)+1] <- f2_nor
    names(ftc_zf2score)[length(ftc_zf2score)] <- PRSs[i]
  }
  f2_models_r_nor[[i]] <- f2_nor
  
  removd <- romr.fnc(f2_models_r_nor[[i]], GENDEP_f2_n, trim = 2.5)$data
  
  f2_nor <- lmer( zf2score ~ bl1fc + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f2_nor)[-1]))==1){
    ftc_zf2score_2[length(ftc_zf2score_2)+1] <- f2_nor
    names(ftc_zf2score_2)[length(ftc_zf2score_2)] <- PRSs[i]
  }
  f2_models_r_nor_2[[i]] <- f2_nor
  
}


# Make a vector to note f3 models_r which failed to converge

ftc_zf3score <- vector(mode = "list", length = 0)
ftc_zf3score_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  f3_both <- lmer( zf3score ~ bl1fn + week + week2 +
                     sex + cage + drug + GENDEP_f3[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = GENDEP_f3, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_both)[-1]))==1){
    ftc_zf3score[length(ftc_zf3score)+1] <- f3_both
    names(ftc_zf3score)[length(ftc_zf3score)] <- PRSs[i]
  }
  f3_models_r_both[[i]] <- f3_both
  
  removd <- romr.fnc(f3_models_r_both[[i]], GENDEP_f3, trim = 2.5)$data
  
  f3_both <- lmer( zf3score ~ bl1fn + week + week2 +
                     sex + cage + drug + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_both)[-1]))==1){
    ftc_zf3score_2[length(ftc_zf3score_2)+1] <- f3_both
    names(ftc_zf3score_2)[length(ftc_zf3score_2)] <- PRSs[i]
  }
  f3_models_r_both_2[[i]] <- f3_both
  
  # Then, fit the models_r for each drug separately
  
  f3_esc <- lmer( zf3score ~ bl1fn + week + week2 +
                    sex + cage + GENDEP_f3_e[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = GENDEP_f3_e, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_esc)[-1]))==1){
    ftc_zf3score[length(ftc_zf3score)+1] <- f3_esc
    names(ftc_zf3score)[length(ftc_zf3score)] <- PRSs[i]
  }
  f3_models_r_esc[[i]] <- f3_esc
  
  removd <- romr.fnc(f3_models_r_esc[[i]], GENDEP_f3_e, trim = 2.5)$data
  
  f3_esc <- lmer( zf3score ~ bl1fn + week + week2 +
                    sex + cage + removd[,PRSs[i]] +
                    zPC1 + zPC2 + zPC3 +
                    (1|centreid) + (1|subjectid),
                  data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_esc)[-1]))==1){
    ftc_zf3score_2[length(ftc_zf3score_2)+1] <- f3_esc
    names(ftc_zf3score_2)[length(ftc_zf3score_2)] <- PRSs[i]
  }
  f3_models_r_esc_2[[i]] <- f3_esc
  
  #######
  
  f3_nor<- lmer( zf3score ~ bl1fn + week + week2 +
                   sex + cage + GENDEP_f3_n[,PRSs[i]] +
                   zPC1 + zPC2 + zPC3 +
                   (1|centreid) + (1|subjectid),
                 data = GENDEP_f3_n, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_nor)[-1]))==1){
    ftc_zf3score[length(ftc_zf3score)+1] <- f3_nor
    names(ftc_zf3score)[length(ftc_zf3score)] <- PRSs[i]
  }
  f3_models_r_nor[[i]] <- f3_nor
  
  removd <- romr.fnc(f3_models_r_nor[[i]], GENDEP_f3_n, trim = 2.5)$data
  
  f3_nor<- lmer( zf3score ~ bl1fn + week + week2 +
                   sex + cage + removd[,PRSs[i]] +
                   zPC1 + zPC2 + zPC3 +
                   (1|centreid) + (1|subjectid),
                 data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(f3_nor)[-1]))==1){
    ftc_zf3score_2[length(ftc_zf3score_2)+1] <- f3_nor
    names(ftc_zf3score_2)[length(ftc_zf3score_2)] <- PRSs[i]
  }
  f3_models_r_nor_2[[i]] <- f3_nor
  
}

# Make a vector to note sui models_r_r which failed to converge

ftc_suiscore <- vector(mode = "list", length = 0)
ftc_suiscore_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  sui_both <- lmer( suiscore ~ bl1sui + week + week2 +
                      sex + cage + drug + GENDEP_sui[,PRSs[i]] +
                      zPC1 + zPC2 + zPC3 +
                      (1|centreid) + (1|subjectid),
                    data = GENDEP_sui, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_both)[-1]))==1){
    ftc_suiscore[length(ftc_suiscore)+1] <- sui_both
    names(ftc_suiscore)[length(ftc_suiscore)] <- PRSs[i]
  }
  sui_models_r_both[[i]] <- sui_both
  
  removd <- romr.fnc(sui_models_r_both[[i]], GENDEP_sui, trim = 2.5)$data
  
  sui_both <- lmer( suiscore ~ bl1sui + week + week2 +
                      sex + cage + drug + removd[,PRSs[i]] +
                      zPC1 + zPC2 + zPC3 +
                      (1|centreid) + (1|subjectid),
                    data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_both)[-1]))==1){
    ftc_suiscore_2[length(ftc_suiscore_2)+1] <- sui_both
    names(ftc_suiscore_2)[length(ftc_suiscore_2)] <- PRSs[i]
  }
  sui_models_r_both_2[[i]] <- sui_both
  
  # Then, fit the models_r for each drug separately
  
  sui_esc <- lmer( suiscore ~ bl1sui + week + week2 +
                     sex + cage + GENDEP_sui_e[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = GENDEP_sui_e, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_esc)[-1]))==1){
    ftc_suiscore[length(ftc_suiscore)+1] <- sui_esc
    names(ftc_suiscore)[length(ftc_suiscore)] <- PRSs[i]
  }
  sui_models_r_esc[[i]] <- sui_esc
  
  removd <- romr.fnc(sui_models_r_esc[[i]], GENDEP_sui_e, trim = 2.5)$data
  
  sui_esc <- lmer( suiscore ~ bl1sui + week + week2 +
                     sex + cage + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_esc)[-1]))==1){
    ftc_suiscore_2[length(ftc_suiscore_2)+1] <- sui_esc
    names(ftc_suiscore_2)[length(ftc_suiscore_2)] <- PRSs[i]
  }
  sui_models_r_esc_2[[i]] <- sui_esc
  
  
  #######
  
  sui_nor <- lmer( suiscore ~ bl1sui + week + week2 +
                     sex + cage + GENDEP_sui_n[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|subjectid),
                   data = GENDEP_sui_n, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_nor)[-1]))==1){
    ftc_suiscore[length(ftc_suiscore)+1] <- sui_nor
    names(ftc_suiscore)[length(ftc_suiscore)] <- PRSs[i]
  }
  sui_models_r_nor[[i]] <- sui_nor
  
  removd <- romr.fnc(sui_models_r_nor[[i]], GENDEP_sui_n, trim = 2.5)$data
  
  sui_nor <- lmer( suiscore ~ bl1sui + week + week2 +
                     sex + cage + removd[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|subjectid),
                   data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(sui_nor)[-1]))==1){
    ftc_suiscore_2[length(ftc_suiscore_2)+1] <- sui_nor
    names(ftc_suiscore_2)[length(ftc_suiscore_2)] <- PRSs[i]
  }
  sui_models_r_nor_2[[i]] <- sui_nor
  
  
}

# Models which failed to converge

# f3score, nortrip-only, NZ_0_5
# suiscore, nortrip-only, ALL THIRTY PRSs! Singular fit issue.

# Look at models without polygenic score

f1_models_r_both[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                                sex + cage + drug + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = GENDEP_f1, REML=FALSE )

removd <- romr.fnc(f1_models_r_both[[31]], GENDEP_f1, trim = 2.5)$data

f1_models_r_both_2[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                                  sex + cage + drug + 
                                  zPC1 + zPC2 + zPC3 + 
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )

##

f1_models_r_esc[[31]] <-  lmer( zf1score ~ bl1fm + week + week2 +
                                sex + cage + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = GENDEP_f1_e, REML=FALSE )

removd <- romr.fnc(f1_models_r_esc[[31]], GENDEP_f1_e, trim = 2.5)$data

f1_models_r_esc_2[[31]] <-  lmer( zf1score ~ bl1fm + week + week2 +
                                  sex + cage + 
                                  zPC1 + zPC2 + zPC3 + 
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )

##

f1_models_r_nor[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = GENDEP_f1_n, REML=FALSE )

removd <- romr.fnc(f1_models_r_nor[[31]], GENDEP_f1_n, trim = 2.5)$data

f1_models_r_nor_2[[31]] <- lmer( zf1score ~ bl1fm + week + week2 +
                                 sex + cage + 
                                 zPC1 + zPC2 + zPC3 + 
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

##

f2_models_r_both[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                                sex + cage + drug + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = GENDEP_f2, REML=FALSE )

removd <- romr.fnc(f2_models_r_both[[31]], GENDEP_f2, trim = 2.5)$data

f2_models_r_both_2[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                                  sex + cage + drug + 
                                  zPC1 + zPC2 + zPC3 + 
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )

##

f2_models_r_esc[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = GENDEP_f2_e, REML=FALSE )

removd <- romr.fnc(f2_models_r_esc[[31]], GENDEP_f2_e, trim = 2.5)$data

f2_models_r_esc_2[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                                 sex + cage + 
                                 zPC1 + zPC2 + zPC3 + 
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

##

f2_models_r_nor[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = GENDEP_f2_n, REML=FALSE )

removd <- romr.fnc(f2_models_r_nor[[31]], GENDEP_f2_n, trim = 2.5)$data

f2_models_r_nor_2[[31]] <- lmer( zf2score ~ bl1fc + week + week2 +
                                 sex + cage + 
                                 zPC1 + zPC2 + zPC3 + 
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

##

f3_models_r_both[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                                sex + cage + drug + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = GENDEP_f3, REML=FALSE )

removd <- romr.fnc(f3_models_r_both[[31]], GENDEP_f3, trim = 2.5)$data

f3_models_r_both_2[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                                  sex + cage + drug + 
                                  zPC1 + zPC2 + zPC3 + 
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )

##

f3_models_r_esc[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = GENDEP_f3_e, REML=FALSE )

removd <- romr.fnc(f3_models_r_esc[[31]], GENDEP_f3_e, trim = 2.5)$data

f3_models_r_esc_2[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                                 sex + cage + 
                                 zPC1 + zPC2 + zPC3 + 
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

##

f3_models_r_nor[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                               sex + cage + 
                               zPC1 + zPC2 + zPC3 + 
                               (1|centreid) + (1|subjectid),
                             data = GENDEP_f3_n, REML=FALSE )

removd <- romr.fnc(f3_models_r_nor[[31]], GENDEP_f3_n, trim = 2.5)$data

f3_models_r_nor_2[[31]] <- lmer( zf3score ~ bl1fn + week + week2 +
                                 sex + cage + 
                                 zPC1 + zPC2 + zPC3 + 
                                 (1|centreid) + (1|subjectid),
                               data = removd, REML=FALSE )

##

sui_models_r_both[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                 sex + cage + drug + 
                                 zPC1 + zPC2 + zPC3 + 
                                 (1|centreid) + (1|subjectid),
                               data = GENDEP_sui, REML=FALSE )

removd <- romr.fnc(sui_models_r_both[[31]], GENDEP_sui, trim = 2.5)$data

sui_models_r_both_2[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                   sex + cage + drug + 
                                   zPC1 + zPC2 + zPC3 + 
                                   (1|centreid) + (1|subjectid),
                                 data = removd, REML=FALSE )

##

sui_models_r_esc[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                sex + cage + 
                                zPC1 + zPC2 + zPC3 + 
                                (1|centreid) + (1|subjectid),
                              data = GENDEP_sui_e, REML=FALSE )

removd <- romr.fnc(sui_models_r_esc[[31]], GENDEP_sui_e, trim = 2.5)$data

sui_models_r_esc_2[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                  sex + cage + 
                                  zPC1 + zPC2 + zPC3 + 
                                  (1|centreid) + (1|subjectid),
                                data = removd, REML=FALSE )

##

sui_models_r_nor[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                sex + cage + 
                                zPC1 + zPC2 + zPC3 +
                                (1|subjectid),
                              data = GENDEP_sui_n, REML=FALSE )

removd <- romr.fnc(sui_models_r_nor[[31]], GENDEP_sui_n, trim = 2.5)$data

sui_models_r_nor_2[[31]] <- lmer( suiscore ~ bl1sui + week + week2 +
                                  sex + cage + 
                                  zPC1 + zPC2 + zPC3 +
                                  (1|subjectid),
                                data = removd, REML=FALSE )

names(f1_models_r_both) <- PRSs_31
names(f1_models_r_esc) <- PRSs_31
names(f1_models_r_nor) <- PRSs_31

names(f2_models_r_both) <- PRSs_31
names(f2_models_r_esc) <- PRSs_31
names(f2_models_r_nor) <- PRSs_31

names(f3_models_r_both) <- PRSs_31
names(f3_models_r_esc) <- PRSs_31
names(f3_models_r_nor) <- PRSs_31

names(sui_models_r_both) <- PRSs_31
names(sui_models_r_esc) <- PRSs_31
names(sui_models_r_nor) <- PRSs_31

names(f1_models_r_both_2) <- PRSs_31
names(f1_models_r_esc_2) <- PRSs_31
names(f1_models_r_nor_2) <- PRSs_31

names(f2_models_r_both_2) <- PRSs_31
names(f2_models_r_esc_2) <- PRSs_31
names(f2_models_r_nor_2) <- PRSs_31

names(f3_models_r_both_2) <- PRSs_31
names(f3_models_r_esc_2) <- PRSs_31
names(f3_models_r_nor_2) <- PRSs_31

names(sui_models_r_both_2) <- PRSs_31
names(sui_models_r_esc_2) <- PRSs_31
names(sui_models_r_nor_2) <- PRSs_31

##### 14c) Get betas, standard errors and p-values for symptom subsets #####

f1_models_r_both_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_esc_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_nor_betas <- vector(mode = "list", length = length(PRSs_31))

f2_models_r_both_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_esc_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_nor_betas <- vector(mode = "list", length = length(PRSs_31))

f3_models_r_both_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_esc_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_nor_betas <- vector(mode = "list", length = length(PRSs_31))

sui_models_r_both_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_esc_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_nor_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(f1_models_r_both)){
  f1_models_r_both_betas[i] <- list(summary(f1_models_r_both[[i]])$coefficients[8,])
}
for (i in 1:length(f1_models_r_esc)){
  f1_models_r_esc_betas[i] <- list(summary(f1_models_r_esc[[i]])$coefficients[7,])
}
for (i in 1:length(f1_models_r_nor)){
  f1_models_r_nor_betas[i] <- list(summary(f1_models_r_nor[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_r_both)){
  f2_models_r_both_betas[i] <- list(summary(f2_models_r_both[[i]])$coefficients[8,])
}
for (i in 1:length(f2_models_r_esc)){
  f2_models_r_esc_betas[i] <- list(summary(f2_models_r_esc[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_r_nor)){
  f2_models_r_nor_betas[i] <- list(summary(f2_models_r_nor[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_r_both)){
  f3_models_r_both_betas[i] <- list(summary(f3_models_r_both[[i]])$coefficients[8,])
}
for (i in 1:length(f3_models_r_esc)){
  f3_models_r_esc_betas[i] <- list(summary(f3_models_r_esc[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_r_nor)){
  f3_models_r_nor_betas[i] <- list(summary(f3_models_r_nor[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_r_both)){
  sui_models_r_both_betas[i] <- list(summary(sui_models_r_both[[i]])$coefficients[8,])
}
for (i in 1:length(sui_models_r_esc)){
  sui_models_r_esc_betas[i] <- list(summary(sui_models_r_esc[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_r_nor)){
  sui_models_r_nor_betas[i] <- list(summary(sui_models_r_nor[[i]])$coefficients[7,])
}

names(f1_models_r_both_betas) <- PRSs_31
names(f1_models_r_esc_betas) <- PRSs_31
names(f1_models_r_nor_betas) <- PRSs_31
names(f2_models_r_both_betas) <- PRSs_31
names(f2_models_r_esc_betas) <- PRSs_31
names(f2_models_r_nor_betas) <- PRSs_31
names(f3_models_r_both_betas) <- PRSs_31
names(f3_models_r_esc_betas) <- PRSs_31
names(f3_models_r_nor_betas) <- PRSs_31
names(sui_models_r_both_betas) <- PRSs_31
names(sui_models_r_esc_betas) <- PRSs_31
names(sui_models_r_nor_betas) <- PRSs_31

f1score_r_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(f1_models_r_both_betas))[,1],
    t(as.data.frame(f1_models_r_both_betas))[,2],
    t(as.data.frame(f1_models_r_both_betas))[,5],
    t(as.data.frame(f1_models_r_esc_betas))[,1],
    t(as.data.frame(f1_models_r_esc_betas))[,2],
    t(as.data.frame(f1_models_r_esc_betas))[,5],
    t(as.data.frame(f1_models_r_nor_betas))[,1],
    t(as.data.frame(f1_models_r_nor_betas))[,2],
    t(as.data.frame(f1_models_r_nor_betas))[,5]
  )
)

f2score_r_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(f2_models_r_both_betas))[,1],
    t(as.data.frame(f2_models_r_both_betas))[,2],
    t(as.data.frame(f2_models_r_both_betas))[,5],
    t(as.data.frame(f2_models_r_esc_betas))[,1],
    t(as.data.frame(f2_models_r_esc_betas))[,2],
    t(as.data.frame(f2_models_r_esc_betas))[,5],
    t(as.data.frame(f2_models_r_nor_betas))[,1],
    t(as.data.frame(f2_models_r_nor_betas))[,2],
    t(as.data.frame(f2_models_r_nor_betas))[,5]
  )
)

f3score_r_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(f3_models_r_both_betas))[,1],
    t(as.data.frame(f3_models_r_both_betas))[,2],
    t(as.data.frame(f3_models_r_both_betas))[,5],
    t(as.data.frame(f3_models_r_esc_betas))[,1],
    t(as.data.frame(f3_models_r_esc_betas))[,2],
    t(as.data.frame(f3_models_r_esc_betas))[,5],
    t(as.data.frame(f3_models_r_nor_betas))[,1],
    t(as.data.frame(f3_models_r_nor_betas))[,2],
    t(as.data.frame(f3_models_r_nor_betas))[,5]
  )
)

suiscore_r_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(sui_models_r_both_betas))[,1],
    t(as.data.frame(sui_models_r_both_betas))[,2],
    t(as.data.frame(sui_models_r_both_betas))[,5],
    t(as.data.frame(sui_models_r_esc_betas))[,1],
    t(as.data.frame(sui_models_r_esc_betas))[,2],
    t(as.data.frame(sui_models_r_esc_betas))[,5],
    t(as.data.frame(sui_models_r_nor_betas))[,1],
    t(as.data.frame(sui_models_r_nor_betas))[,2],
    t(as.data.frame(sui_models_r_nor_betas))[,5]
  )
)

colnames( f1score_r_beta_matrix ) <- c(
  "f1score_both_B", "f1score_both_SE","f1score_both_p",
  "f1score_esc_B", "f1score_esc_SE", "f1score_esc_p",
  "f1score_nor_B", "f1score_nor_SE", "f1score_nor_p"
)

colnames( f2score_r_beta_matrix ) <- c(
  "f2score_both_B", "f2score_both_SE","f2score_both_p",
  "f2score_esc_B", "f2score_esc_SE", "f2score_esc_p",
  "f2score_nor_B", "f2score_nor_SE", "f2score_nor_p"
)

colnames( f3score_r_beta_matrix ) <- c(
  "f3score_both_B", "f3score_both_SE","f3score_both_p",
  "f3score_esc_B", "f3score_esc_SE", "f3score_esc_p",
  "f3score_nor_B", "f3score_nor_SE", "f3score_nor_p"
)

colnames( suiscore_r_beta_matrix ) <- c(
  "suiscore_both_B", "suiscore_both_SE","suiscore_both_p",
  "suiscore_esc_B", "suiscore_esc_SE", "suiscore_esc_p",
  "suiscore_nor_B", "suiscore_nor_SE", "suiscore_nor_p"
)

##

f1_models_r_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
f1_models_r_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

f2_models_r_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
f2_models_r_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

f3_models_r_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
f3_models_r_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

sui_models_r_both_2_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_esc_2_betas <- vector(mode = "list", length = length(PRSs_31))
sui_models_r_nor_2_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(f1_models_r_both_2)){
  f1_models_r_both_2_betas[i] <- list(summary(f1_models_r_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(f1_models_r_esc_2)){
  f1_models_r_esc_2_betas[i] <- list(summary(f1_models_r_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(f1_models_r_nor_2)){
  f1_models_r_nor_2_betas[i] <- list(summary(f1_models_r_nor_2[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_r_both_2)){
  f2_models_r_both_2_betas[i] <- list(summary(f2_models_r_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(f2_models_r_esc_2)){
  f2_models_r_esc_2_betas[i] <- list(summary(f2_models_r_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(f2_models_r_nor_2)){
  f2_models_r_nor_2_betas[i] <- list(summary(f2_models_r_nor_2[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_r_both_2)){
  f3_models_r_both_2_betas[i] <- list(summary(f3_models_r_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(f3_models_r_esc_2)){
  f3_models_r_esc_2_betas[i] <- list(summary(f3_models_r_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(f3_models_r_nor_2)){
  f3_models_r_nor_2_betas[i] <- list(summary(f3_models_r_nor_2[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_r_both_2)){
  sui_models_r_both_2_betas[i] <- list(summary(sui_models_r_both_2[[i]])$coefficients[8,])
}
for (i in 1:length(sui_models_r_esc_2)){
  sui_models_r_esc_2_betas[i] <- list(summary(sui_models_r_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(sui_models_r_nor_2)){
  sui_models_r_nor_2_betas[i] <- list(summary(sui_models_r_nor_2[[i]])$coefficients[7,])
}

names(f1_models_r_both_2_betas) <- PRSs_31
names(f1_models_r_esc_2_betas) <- PRSs_31
names(f1_models_r_nor_2_betas) <- PRSs_31
names(f2_models_r_both_2_betas) <- PRSs_31
names(f2_models_r_esc_2_betas) <- PRSs_31
names(f2_models_r_nor_2_betas) <- PRSs_31
names(f3_models_r_both_2_betas) <- PRSs_31
names(f3_models_r_esc_2_betas) <- PRSs_31
names(f3_models_r_nor_2_betas) <- PRSs_31
names(sui_models_r_both_2_betas) <- PRSs_31
names(sui_models_r_esc_2_betas) <- PRSs_31
names(sui_models_r_nor_2_betas) <- PRSs_31

f1score_r_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(f1_models_r_both_2_betas))[,1],
    t(as.data.frame(f1_models_r_both_2_betas))[,2],
    t(as.data.frame(f1_models_r_both_2_betas))[,5],
    t(as.data.frame(f1_models_r_esc_2_betas))[,1],
    t(as.data.frame(f1_models_r_esc_2_betas))[,2],
    t(as.data.frame(f1_models_r_esc_2_betas))[,5],
    t(as.data.frame(f1_models_r_nor_2_betas))[,1],
    t(as.data.frame(f1_models_r_nor_2_betas))[,2],
    t(as.data.frame(f1_models_r_nor_2_betas))[,5]
  )
)

f2score_r_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(f2_models_r_both_2_betas))[,1],
    t(as.data.frame(f2_models_r_both_2_betas))[,2],
    t(as.data.frame(f2_models_r_both_2_betas))[,5],
    t(as.data.frame(f2_models_r_esc_2_betas))[,1],
    t(as.data.frame(f2_models_r_esc_2_betas))[,2],
    t(as.data.frame(f2_models_r_esc_2_betas))[,5],
    t(as.data.frame(f2_models_r_nor_2_betas))[,1],
    t(as.data.frame(f2_models_r_nor_2_betas))[,2],
    t(as.data.frame(f2_models_r_nor_2_betas))[,5]
  )
)

f3score_r_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(f3_models_r_both_2_betas))[,1],
    t(as.data.frame(f3_models_r_both_2_betas))[,2],
    t(as.data.frame(f3_models_r_both_2_betas))[,5],
    t(as.data.frame(f3_models_r_esc_2_betas))[,1],
    t(as.data.frame(f3_models_r_esc_2_betas))[,2],
    t(as.data.frame(f3_models_r_esc_2_betas))[,5],
    t(as.data.frame(f3_models_r_nor_2_betas))[,1],
    t(as.data.frame(f3_models_r_nor_2_betas))[,2],
    t(as.data.frame(f3_models_r_nor_2_betas))[,5]
  )
)

suiscore_r_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(sui_models_r_both_2_betas))[,1],
    t(as.data.frame(sui_models_r_both_2_betas))[,2],
    t(as.data.frame(sui_models_r_both_2_betas))[,5],
    t(as.data.frame(sui_models_r_esc_2_betas))[,1],
    t(as.data.frame(sui_models_r_esc_2_betas))[,2],
    t(as.data.frame(sui_models_r_esc_2_betas))[,5],
    t(as.data.frame(sui_models_r_nor_2_betas))[,1],
    t(as.data.frame(sui_models_r_nor_2_betas))[,2],
    t(as.data.frame(sui_models_r_nor_2_betas))[,5]
  )
)

colnames( f1score_r_beta_matrix_2 ) <- c(
  "f1score_both_2_B", "f1score_both_2_SE","f1score_both_2_p",
  "f1score_esc_2_B", "f1score_esc_2_SE", "f1score_esc_2_p",
  "f1score_nor_2_B", "f1score_nor_2_SE", "f1score_nor_2_p"
)

colnames( f2score_r_beta_matrix_2 ) <- c(
  "f2score_both_2_B", "f2score_both_2_SE","f2score_both_2_p",
  "f2score_esc_2_B", "f2score_esc_2_SE", "f2score_esc_2_p",
  "f2score_nor_2_B", "f2score_nor_2_SE", "f2score_nor_2_p"
)

colnames( f3score_r_beta_matrix_2 ) <- c(
  "f3score_both_2_B", "f3score_both_2_SE","f3score_both_2_p",
  "f3score_esc_2_B", "f3score_esc_2_SE", "f3score_esc_2_p",
  "f3score_nor_2_B", "f3score_nor_2_SE", "f3score_nor_2_p"
)

colnames( suiscore_r_beta_matrix_2 ) <- c(
  "suiscore_both_2_B", "suiscore_both_2_SE","suiscore_both_2_p",
  "suiscore_esc_2_B", "suiscore_esc_2_SE", "suiscore_esc_2_p",
  "suiscore_nor_2_B", "suiscore_nor_2_SE", "suiscore_nor_2_p"
)

##### 14d) Get PRSs' explained variances for symptom subsets #####

f1_matrix_r_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_r_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_r_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_r_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_r_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_r_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_r_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_r_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_r_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_r_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_r_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_r_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

for (i in 1:31){
  f1_matrix_r_both[i,] <- r.squaredGLMM(f1_models_r_both[[i]])
  f1_matrix_r_esc[i,] <- r.squaredGLMM(f1_models_r_esc[[i]])
  f1_matrix_r_nor[i,] <- r.squaredGLMM(f1_models_r_nor[[i]])
}

for (i in 1:31){
  f2_matrix_r_both[i,] <- r.squaredGLMM(f2_models_r_both[[i]])
  f2_matrix_r_esc[i,] <- r.squaredGLMM(f2_models_r_esc[[i]])
  f2_matrix_r_nor[i,] <- r.squaredGLMM(f2_models_r_nor[[i]])
}

for (i in 1:31){
  f3_matrix_r_both[i,] <- r.squaredGLMM(f3_models_r_both[[i]])
  f3_matrix_r_esc[i,] <- r.squaredGLMM(f3_models_r_esc[[i]])
  f3_matrix_r_nor[i,] <- r.squaredGLMM(f3_models_r_nor[[i]])
}

for (i in 1:31){
  sui_matrix_r_both[i,] <- r.squaredGLMM(sui_models_r_both[[i]])
  sui_matrix_r_esc[i,] <- r.squaredGLMM(sui_models_r_esc[[i]])
  sui_matrix_r_nor[i,] <- r.squaredGLMM(sui_models_r_nor[[i]])
}

##

f1_matrix_r_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_r_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f1_matrix_r_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_r_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_r_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f2_matrix_r_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_r_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_r_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

f3_matrix_r_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_r_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_r_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

sui_matrix_r_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

for (i in 1:31){
  f1_matrix_r_both_2[i,] <- r.squaredGLMM(f1_models_r_both_2[[i]])
  f1_matrix_r_esc_2[i,] <- r.squaredGLMM(f1_models_r_esc_2[[i]])
  f1_matrix_r_nor_2[i,] <- r.squaredGLMM(f1_models_r_nor_2[[i]])
}

for (i in 1:31){
  f2_matrix_r_both_2[i,] <- r.squaredGLMM(f2_models_r_both_2[[i]])
  f2_matrix_r_esc_2[i,] <- r.squaredGLMM(f2_models_r_esc_2[[i]])
  f2_matrix_r_nor_2[i,] <- r.squaredGLMM(f2_models_r_nor_2[[i]])
}

for (i in 1:31){
  f3_matrix_r_both_2[i,] <- r.squaredGLMM(f3_models_r_both_2[[i]])
  f3_matrix_r_esc_2[i,] <- r.squaredGLMM(f3_models_r_esc_2[[i]])
  f3_matrix_r_nor_2[i,] <- r.squaredGLMM(f3_models_r_nor_2[[i]])
}

for (i in 1:31){
  sui_matrix_r_both_2[i,] <- r.squaredGLMM(sui_models_r_both_2[[i]])
  sui_matrix_r_esc_2[i,] <- r.squaredGLMM(sui_models_r_esc_2[[i]])
  sui_matrix_r_nor_2[i,] <- r.squaredGLMM(sui_models_r_nor_2[[i]])
}


##### 14e) Make output files for symptom subsets #####

sink( "f1score_betas_rand_only.csv")

# Betas and p-values

cat("f1score_betas\n,")
write.table( f1score_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(f1_matrix_r_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(f1_matrix_r_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(f1_matrix_r_nor,
            col.names=TRUE, sep=",")

sink()

sink( "f2score_betas_rand_only.csv")

# Betas and p-values

cat("f2score_betas\n,")
write.table( f2score_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(f2_matrix_r_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(f2_matrix_r_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(f2_matrix_r_nor,
            col.names=TRUE, sep=",")

sink()

sink( "f3score_betas_rand_only.csv")

# Betas and p-values

cat("f3score_betas\n,")
write.table( f3score_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(f3_matrix_r_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(f3_matrix_r_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(f3_matrix_r_nor,
            col.names=TRUE, sep=",")

sink()

sink( "suiscore_betas_rand_only.csv")

# Betas and p-values

cat("suiscore_betas\n,")
write.table( suiscore_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(sui_matrix_r_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(sui_matrix_r_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(sui_matrix_r_nor,
            col.names=TRUE, sep=",")

sink()

###

sink( "f1score_betas_rand_onlyinfl_outl_removed.csv")

# Betas and p-values

cat("f1score_betas\n,")
write.table( f1score_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(f1_matrix_r_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(f1_matrix_r_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(f1_matrix_r_nor_2,
            col.names=TRUE, sep=",")

sink()

sink( "f2score_betas_rand_onlyinfl_outl_removed.csv")

# Betas and p-values

cat("f2score_betas\n,")
write.table( f2score_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(f2_matrix_r_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(f2_matrix_r_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(f2_matrix_r_nor_2,
            col.names=TRUE, sep=",")

sink()

sink( "f3score_betas_rand_onlyinfl_outl_removed.csv")

# Betas and p-values

cat("f3score_betas\n,")
write.table( f3score_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(f3_matrix_r_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(f3_matrix_r_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(f3_matrix_r_nor_2,
            col.names=TRUE, sep=",")

sink()

sink( "suiscore_betas_rand_onlyinfl_outl_removed.csv")

# Betas and p-values

cat("suiscore_betas\n,")
write.table( suiscore_r_beta_matrix, col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(sui_matrix_r_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(sui_matrix_r_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(sui_matrix_r_nor_2,
            col.names=TRUE, sep=",")

sink()

##### 14f) Make heatmaps of effect sizes for symptom subsets #####

f1score_betas_only <- cbind(f1score_r_beta_matrix[,1],
                            f1score_r_beta_matrix[,4],
                            f1score_r_beta_matrix[,7])

colnames(f1score_betas_only) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_betas_only <- f1score_betas_only[1:30,]

f1score_ses_only <- cbind(f1score_r_beta_matrix[,2],
                          f1score_r_beta_matrix[,5],
                          f1score_r_beta_matrix[,8])

colnames(f1score_ses_only) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_ses_only <- f1score_ses_only[1:30,]

f1score_p_only <- cbind(f1score_r_beta_matrix[,3],
                        f1score_r_beta_matrix[,6],
                        f1score_r_beta_matrix[,9])

colnames(f1score_p_only) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_p_only <- f1score_p_only[1:30,]

f1_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f1_lowest_ps[which(f1score_p_only<=0.05)] <- 
  f1score_p_only[which(f1score_p_only<=0.05)]

heatmap.2( f1score_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f1_lowest_ps,5),notecol="black",
           main="F1 (mood) score, randomized-participants-only"
)

##

f1score_betas_only_2 <- cbind(f1score_r_beta_matrix_2[,1],
                            f1score_r_beta_matrix_2[,4],
                            f1score_r_beta_matrix_2[,7])

colnames(f1score_betas_only_2) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_betas_only_2 <- f1score_betas_only_2[1:30,]

f1score_ses_only_2 <- cbind(f1score_r_beta_matrix_2[,2],
                          f1score_r_beta_matrix_2[,5],
                          f1score_r_beta_matrix_2[,8])

colnames(f1score_ses_only_2) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_ses_only_2 <- f1score_ses_only_2[1:30,]

f1score_p_only_2 <- cbind(f1score_r_beta_matrix_2[,3],
                        f1score_r_beta_matrix_2[,6],
                        f1score_r_beta_matrix_2[,9])

colnames(f1score_p_only_2) <- c("f1score_both", "f1score_esc","f1score_nor")

f1score_p_only_2 <- f1score_p_only_2[1:30,]

f1_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f1_lowest_ps_2[which(f1score_p_only_2<=0.05)] <- 
  f1score_p_only_2[which(f1score_p_only_2<=0.05)]

heatmap.2( f1score_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f1_lowest_ps_2,5),notecol="black",
           main="F1 (mood) score, randomized-participants-only\ninf out rem"
)

##

f2score_betas_only <- cbind(f2score_r_beta_matrix[,1],
                            f2score_r_beta_matrix[,4],
                            f2score_r_beta_matrix[,7])

colnames(f2score_betas_only) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_betas_only <- f2score_betas_only[1:30,]

f2score_ses_only <- cbind(f2score_r_beta_matrix[,2],
                          f2score_r_beta_matrix[,5],
                          f2score_r_beta_matrix[,8])

colnames(f2score_ses_only) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_ses_only <- f2score_ses_only[1:30,]

f2score_p_only <- cbind(f2score_r_beta_matrix[,3],
                        f2score_r_beta_matrix[,6],
                        f2score_r_beta_matrix[,9])

colnames(f2score_p_only) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_p_only <- f2score_p_only[1:30,]

f2_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f2_lowest_ps[which(f2score_p_only<=0.05)] <- 
  f2score_p_only[which(f2score_p_only<=0.05)]

heatmap.2( f2score_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f2_lowest_ps,5),notecol="black",
           main="F2 (cognitive) score,\nrandomized-participants-only"
)

##

f2score_betas_only_2 <- cbind(f2score_r_beta_matrix_2[,1],
                            f2score_r_beta_matrix_2[,4],
                            f2score_r_beta_matrix_2[,7])

colnames(f2score_betas_only_2) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_betas_only_2 <- f2score_betas_only_2[1:30,]

f2score_ses_only_2 <- cbind(f2score_r_beta_matrix_2[,2],
                          f2score_r_beta_matrix_2[,5],
                          f2score_r_beta_matrix_2[,8])

colnames(f2score_ses_only_2) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_ses_only_2 <- f2score_ses_only_2[1:30,]

f2score_p_only_2 <- cbind(f2score_r_beta_matrix_2[,3],
                        f2score_r_beta_matrix_2[,6],
                        f2score_r_beta_matrix_2[,9])

colnames(f2score_p_only_2) <- c("f2score_both", "f2score_esc","f2score_nor")

f2score_p_only_2 <- f2score_p_only_2[1:30,]

f2_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f2_lowest_ps_2[which(f2score_p_only_2<=0.05)] <- 
  f2score_p_only_2[which(f2score_p_only_2<=0.05)]

heatmap.2( f2score_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f2_lowest_ps_2,5),notecol="black",
           main="F2 (cognitive) score,\nrandomized-participants-only\ninf out rem"
)

##

f3score_betas_only <- cbind(f3score_r_beta_matrix[,1],
                            f3score_r_beta_matrix[,4],
                            f3score_r_beta_matrix[,7])

colnames(f3score_betas_only) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_betas_only <- f3score_betas_only[1:30,]

f3score_ses_only <- cbind(f3score_r_beta_matrix[,2],
                          f3score_r_beta_matrix[,5],
                          f3score_r_beta_matrix[,8])

colnames(f3score_ses_only) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_ses_only <- f3score_ses_only[1:30,]

f3score_p_only <- cbind(f3score_r_beta_matrix[,3],
                        f3score_r_beta_matrix[,6],
                        f3score_r_beta_matrix[,9])

colnames(f3score_p_only) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_p_only <- f3score_p_only[1:30,]

f3_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f3_lowest_ps[which(f3score_p_only<=0.05)] <- 
  f3score_p_only[which(f3score_p_only<=0.05)]

heatmap.2( f3score_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f3_lowest_ps,5),notecol="black",
           main="F3 (neurovegetative) score,\nrandomized-participants-only"
)

##

f3score_betas_only_2 <- cbind(f3score_r_beta_matrix_2[,1],
                            f3score_r_beta_matrix_2[,4],
                            f3score_r_beta_matrix_2[,7])

colnames(f3score_betas_only_2) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_betas_only_2 <- f3score_betas_only_2[1:30,]

f3score_ses_only_2 <- cbind(f3score_r_beta_matrix_2[,2],
                          f3score_r_beta_matrix_2[,5],
                          f3score_r_beta_matrix_2[,8])

colnames(f3score_ses_only_2) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_ses_only_2 <- f3score_ses_only_2[1:30,]

f3score_p_only_2 <- cbind(f3score_r_beta_matrix_2[,3],
                        f3score_r_beta_matrix_2[,6],
                        f3score_r_beta_matrix_2[,9])

colnames(f3score_p_only_2) <- c("f3score_both", "f3score_esc","f3score_nor")

f3score_p_only_2 <- f3score_p_only_2[1:30,]

f3_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                    c("both", "esc", "nor")))

f3_lowest_ps_2[which(f3score_p_only_2<=0.05)] <- 
  f3score_p_only_2[which(f3score_p_only_2<=0.05)]

heatmap.2( f3score_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(f3_lowest_ps_2,5),notecol="black",
           main="F3 (neurovegetative) score,\nrandomized-participants-only\ninf out rem"
)

##

suiscore_betas_only <- cbind(suiscore_r_beta_matrix[,1],
                             suiscore_r_beta_matrix[,4],
                             suiscore_r_beta_matrix[,7])

colnames(suiscore_betas_only) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_betas_only <- suiscore_betas_only[1:30,]

suiscore_ses_only <- cbind(suiscore_r_beta_matrix[,2],
                           suiscore_r_beta_matrix[,5],
                           suiscore_r_beta_matrix[,8])

colnames(suiscore_ses_only) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_ses_only <- suiscore_ses_only[1:30,]

suiscore_p_only <- cbind(suiscore_r_beta_matrix[,3],
                         suiscore_r_beta_matrix[,6],
                         suiscore_r_beta_matrix[,9])

colnames(suiscore_p_only) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_p_only <- suiscore_p_only[1:30,]

sui_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                     c("both", "esc", "nor")))

sui_lowest_ps[which(suiscore_p_only<=0.05)] <- 
  suiscore_p_only[which(suiscore_p_only<=0.05)]

heatmap.2( suiscore_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(sui_lowest_ps,5),notecol="black",
           main="Suicidality score,\nrandomized-participants-only"
)

##

suiscore_betas_only_2 <- cbind(suiscore_r_beta_matrix_2[,1],
                             suiscore_r_beta_matrix_2[,4],
                             suiscore_r_beta_matrix_2[,7])

colnames(suiscore_betas_only_2) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_betas_only_2 <- suiscore_betas_only_2[1:30,]

suiscore_ses_only_2 <- cbind(suiscore_r_beta_matrix_2[,2],
                           suiscore_r_beta_matrix_2[,5],
                           suiscore_r_beta_matrix_2[,8])

colnames(suiscore_ses_only_2) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_ses_only_2 <- suiscore_ses_only_2[1:30,]

suiscore_p_only_2 <- cbind(suiscore_r_beta_matrix_2[,3],
                         suiscore_r_beta_matrix_2[,6],
                         suiscore_r_beta_matrix_2[,9])

colnames(suiscore_p_only_2) <- c("suiscore_both", "suiscore_esc","suiscore_nor")

suiscore_p_only_2 <- suiscore_p_only_2[1:30,]

sui_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                     c("both", "esc", "nor")))

sui_lowest_ps_2[which(suiscore_p_only_2<=0.05)] <- 
  suiscore_p_only_2[which(suiscore_p_only_2<=0.05)]

heatmap.2( suiscore_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(sui_lowest_ps_2,5),notecol="black",
           main="Suicidality score,\nrandomized-participants-only\ninf out rem"
)

##### 15) ANALYSE EFFECTS OF POLYGENIC SCORE ON TOTAL SIDE EFFECT BURDEN (TOTASEC) #####
##### 15a) Build models for totasec #####

# Generally all variables are included for theoretical reasons
# and because they have been included in previous studies.
# PCs 1-3 were important in screeplot when QC-ing SNP data.
# People report side effects more when zf1score (mood) is worse/higher
# - Uher et al (2009) "Adverse reactions to antidepressants".
# I did test using dose as a covariate but ultimately I 
# decided not to because, due to missing data, 66 people become excluded;
# also, Uher et al (2009) found it to be non-significant when time was
# taken into account (and weak when time wasn't)

# Below code verifies that week2 is not important in "both drugs" or
# "nortriptyline", but it is nominally significant in "escitalopram",
# and also that zf1score is significant

mod1 <- lmer( ztotasec ~ zf1score + week + week2 +
                sex + cage + drug + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

drop1(mod1) # week2 is not significant

mod2 <- lmer( ztotasec ~ zf1score + week + 
                sex + cage + drug + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP, REML=FALSE )

drop1(mod2)

hist(residuals(mod2)) 
predictmeans::residplot(mod2, newwd=FALSE)

mod1 <- lmer( ztotasec ~ zf1score + week + week2 +
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP_escit, REML=FALSE )

drop1(mod1) # week2 has p slightly >0.05; can we drop it?

mod2 <- lmer( ztotasec ~ zf1score + week + 
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP_escit, REML=FALSE )

drop1(mod2) # week becomes non-significant when we drop week2
# so go back to using week2

hist(residuals(mod1)) 
predictmeans::residplot(mod1, newwd=FALSE)

mod1 <- lmer( ztotasec ~ zf1score + week + week2 +
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP_nortrip, REML=FALSE )

drop1(mod1) # remove week2 and see if week becomes significant

mod2 <- lmer( ztotasec ~ zf1score + week + 
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GENDEP_nortrip, REML=FALSE )

drop1(mod2) # yes - remove it!

hist(residuals(mod2)) 
predictmeans::residplot(mod2, newwd=FALSE)

##### 15b) Fit mixed effects models for totasec #####

# Make lists to store models in

a_models_both <- vector(mode = "list", length = length(PRSs_31))
a_models_esc <- vector(mode = "list", length = length(PRSs_31))
a_models_nor <- vector(mode = "list", length = length(PRSs_31))

a_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))

# Make new data set without missing data in important variables

GENDEP_a <- subset( GENDEP, !is.na(totasec) )
GENDEP_a_e <- subset( GENDEP_a, drug==0 )
GENDEP_a_n <- subset( GENDEP_a, drug==1 )

# Make a vector to note models which failed to converge

ftc_totasec <- vector(mode = "list", length = 0)
ftc_totasec_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a <- lmer( ztotasec ~ zf1score + week + 
                           sex + cage + drug + GENDEP_a[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = GENDEP_a, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a)[-1]))==1){
    ftc_totasec[length(ftc_totasec)+1] <- lm_both_a
    names(ftc_totasec)[length(ftc_totasec)] <- PRSs[i]
  }
  a_models_both[[i]] <- lm_both_a
  
  removd <- romr.fnc(a_models_both[[i]], GENDEP_a, trim = 2.5)$data # remove influential outliers
  
  lm_both_a <- lmer( ztotasec ~ zf1score + week + 
                           sex + cage + drug + removd[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 +
                           (1|centreid) + (1|subjectid),
                         data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a)[-1]))==1){
    ftc_totasec_2[length(ftc_totasec_2)+1] <- lm_both_a
    names(ftc_totasec_2)[length(ftc_totasec_2)] <- PRSs[i]
  }
  a_models_both_2[[i]] <- lm_both_a
  

  lm_esc_a<- lmer( ztotasec ~ zf1score + week + week2 +
                          sex + cage + GENDEP_a_e[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 +
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_a_e, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a)[-1]))==1){
    ftc_totasec[length(ftc_totasec)+1] <- lm_esc_a
    names(ftc_totasec)[length(ftc_totasec)] <- PRSs[i]
  }
  a_models_esc[[i]] <- lm_esc_a
  
  removd <- romr.fnc(a_models_esc[[i]], GENDEP_a_e, trim = 2.5)$data # remove influential outliers
  
  lm_esc_a <- lmer( ztotasec ~ zf1score + week + week2 +
                          sex + cage + GENDEP_a_e[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 +
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_a_e, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a)[-1]))==1){
    ftc_totasec_2[length(ftc_totasec_2)+1] <- lm_esc_a
    names(ftc_totasec_2)[length(ftc_totasec_2)] <- PRSs[i]
  }
  a_models_esc_2[[i]] <- lm_esc_a
  

  lm_nor_a <- lmer( ztotasec ~ zf1score + week + 
                          sex + cage + GENDEP_a_n[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 +
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_a_n, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a)[-1]))==1){
    ftc_totasec[length(ftc_totasec)+1] <- lm_nor_a
    names(ftc_totasec)[length(ftc_totasec)] <- PRSs[i]
  }
  a_models_nor[[i]] <- lm_nor_a
  
  removd <- romr.fnc(a_models_nor[[i]], GENDEP_a_n, trim = 2.5)$data # remove influential outliers
  
  lm_nor_a <- lmer( ztotasec ~ zf1score + week + 
                      sex + cage + removd[,PRSs[i]] +
                      zPC1 + zPC2 + zPC3 +
                      (1|centreid) + (1|subjectid),
                    data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a)[-1]))==1){
    ftc_totasec_2[length(ftc_totasec_2)+1] <- lm_nor_a
    names(ftc_totasec_2)[length(ftc_totasec_2)] <- PRSs[i]
  }
  a_models_nor_2[[i]] <- lm_nor_a
  
}

names(a_models_both) <- PRSs_31
names(a_models_esc) <- PRSs_31
names(a_models_nor) <- PRSs_31

names(a_models_both_2) <- PRSs_31
names(a_models_esc_2) <- PRSs_31
names(a_models_nor_2) <- PRSs_31

# Look at zmadrs models without polygenic score

noPRS_both_a <- lmer( ztotasec ~ zf1score + week + 
                        sex + cage + drug + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_a, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a)[-1]))==1){
  ftc_totasec[length(ftc_totasec)+1] <- noPRS_both_a
  names(ftc_totasec)[length(ftc_totasec)] <- "No_PRS"
}
a_models_both[[31]] <- noPRS_both_a

removd <- romr.fnc(a_models_both[[31]], GENDEP_a, trim = 2.5)$data # remove influential outliers

noPRS_both_a <- lmer( ztotasec ~ zf1score + week + 
                        sex + cage + drug + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a)[-1]))==1){
  ftc_totasec_2[length(ftc_totasec_2)+1] <- noPRS_both_a
  names(ftc_totasec_2)[length(ftc_totasec_2)] <- "No_PRS"
}
a_models_both_2[[31]] <- noPRS_both_a

noPRS_esc_a <- lmer( ztotasec ~ zf1score + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_a_e, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a)[-1]))==1){
  ftc_totasec[length(ftc_totasec)+1] <- noPRS_esc_a
  names(ftc_totasec)[length(ftc_totasec)] <- "No_PRS"
}
a_models_esc[[31]] <- noPRS_esc_a

removd <- romr.fnc(a_models_esc[[31]], GENDEP_a_e, trim = 2.5)$data # remove influential outliers

noPRS_esc_a <- lmer( ztotasec ~ zf1score + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a)[-1]))==1){
  ftc_totasec_2[length(ftc_totasec_2)+1] <- noPRS_esc_a
  names(ftc_totasec_2)[length(ftc_totasec_2)] <- "No_PRS"
}
a_models_esc_2[[31]] <- noPRS_esc_a

noPRS_nor_a <- lmer( ztotasec ~ zf1score + week + 
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_a_n, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a)[-1]))==1){
  ftc_totasec[length(ftc_totasec)+1] <- noPRS_nor_a
  names(ftc_totasec)[length(ftc_totasec)] <- "No_PRS"
}
a_models_nor[[31]] <- noPRS_nor_a

removd <- romr.fnc(a_models_nor[[31]], GENDEP_a_n, trim = 2.5)$data # remove influential outliers

noPRS_nor_a <- lmer( ztotasec ~ zf1score + week + 
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a)[-1]))==1){
  ftc_totasec_2[length(ftc_totasec_2)+1] <- noPRS_nor_a
  names(ftc_totasec_2)[length(ftc_totasec_2)] <- "No_PRS"
}
a_models_nor_2[[31]] <- noPRS_nor_a

##### 15c) Get betas, standard errors and p-values for totasec #####

a_models_both_betas <- vector(mode = "list", length = length(PRSs_31))
a_models_esc_betas <- vector(mode = "list", length = length(PRSs_31))
a_models_nor_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  a_models_both_betas[i] <- list(summary(a_models_both[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_esc_betas[i] <- list(summary(a_models_esc[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_nor_betas[i] <- list(summary(a_models_nor[[i]])$coefficients[6,])
}

names(a_models_both_betas) <- PRSs_31
names(a_models_esc_betas) <- PRSs_31
names(a_models_nor_betas) <- PRSs_31

totasec_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a_models_both_betas))[,1],
    t(as.data.frame(a_models_both_betas))[,2],
    t(as.data.frame(a_models_both_betas))[,5],
    t(as.data.frame(a_models_esc_betas))[,1],
    t(as.data.frame(a_models_esc_betas))[,2],
    t(as.data.frame(a_models_esc_betas))[,5],
    t(as.data.frame(a_models_nor_betas))[,1],
    t(as.data.frame(a_models_nor_betas))[,2],
    t(as.data.frame(a_models_nor_betas))[,5]
  )
)

colnames( totasec_beta_matrix ) <- c(
  "totasec_both_B", "totasec_both_SE","totasec_both_p",
  "totasec_esc_B", "totasec_esc_SE", "totasec_esc_p",
  "totasec_nor_B", "totasec_nor_SE", "totasec_nor_p"
)

####

a_models_both_betas_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_esc_betas_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_nor_betas_2 <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  a_models_both_betas_2[i] <- list(summary(a_models_both_2[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_esc_betas_2[i] <- list(summary(a_models_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_nor_betas_2[i] <- list(summary(a_models_nor_2[[i]])$coefficients[6,])
}

names(a_models_both_betas_2) <- PRSs_31
names(a_models_esc_betas_2) <- PRSs_31
names(a_models_nor_betas_2) <- PRSs_31

totasec_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(a_models_both_betas_2))[,1],
    t(as.data.frame(a_models_both_betas_2))[,2],
    t(as.data.frame(a_models_both_betas_2))[,5],
    t(as.data.frame(a_models_esc_betas_2))[,1],
    t(as.data.frame(a_models_esc_betas_2))[,2],
    t(as.data.frame(a_models_esc_betas_2))[,5],
    t(as.data.frame(a_models_nor_betas_2))[,1],
    t(as.data.frame(a_models_nor_betas_2))[,2],
    t(as.data.frame(a_models_nor_betas_2))[,5]
  )
)

colnames( totasec_beta_matrix_2 ) <- c(
  "totasec_both_B", "totasec_both_SE","totasec_both_p",
  "totasec_esc_B", "totasec_esc_SE", "totasec_esc_p",
  "totasec_nor_B", "totasec_nor_SE", "totasec_nor_p"
)

##### 15d) Get PRSs' explained variances for totasec #####

a_matrix_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

a_matrix_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

a_matrix_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  a_matrix_both[i,] <- r.squaredGLMM(a_models_both[[i]])
  a_matrix_esc[i,] <- r.squaredGLMM(a_models_esc[[i]])
  a_matrix_nor[i,] <- r.squaredGLMM(a_models_nor[[i]])
}


####

a_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

a_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

a_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  a_matrix_both_2[i,] <- r.squaredGLMM(a_models_both_2[[i]])
  a_matrix_esc_2[i,] <- r.squaredGLMM(a_models_esc_2[[i]])
  a_matrix_nor_2[i,] <- r.squaredGLMM(a_models_nor_2[[i]])
}

##### 15e) Make output files for totasec #####

sink( "totasec.csv")

# Betas and p-values

cat("totasec_betas\n,")
write.table( totasec_beta_matrix[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(a_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a_matrix_nor,
            col.names=TRUE, sep=",")

sink()

####

sink( "totasec_inf_outliers_removed.csv")

# Betas and p-values

cat("totasec_betas\n,")
write.table( totasec_beta_matrix_2[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(a_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()

##### 15f) Make heatmap of effect sizes for totasec #####

a_betas_only <- cbind(
  totasec_beta_matrix[,1],
  totasec_beta_matrix[,4],
  totasec_beta_matrix[,7]
)

colnames(a_betas_only) <- c("a_both_B_r",
                            "a_esc_B_r",
                            "a_nor_B_r")

a_betas_only <- a_betas_only[1:30,]

a_SEs_only <- cbind(
  totasec_beta_matrix[,2],
  totasec_beta_matrix[,5],
  totasec_beta_matrix[,8]
)

colnames(a_SEs_only) <- c("a_both_SE_r",
                          "a_esc_SE_r",
                          "a_nor_SE_r")

a_SEs_only <- a_SEs_only[1:30,]

a_ps_only <- cbind(
  totasec_beta_matrix[,3],
  totasec_beta_matrix[,6],
  totasec_beta_matrix[,9]
)

colnames(a_ps_only) <- c("a_both_p_r",
                         "a_esc_p_r",
                         "a_nor_p_r")

a_ps_only <- a_ps_only[1:30,]

a_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                   c("both", "esc", "nor")))

a_lowest_ps[which(a_ps_only<=0.05)] <- 
  a_ps_only[which(a_ps_only<=0.05)]

heatmap.2( a_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(a_lowest_ps,5),notecol="black",
           main="Total ASEC"
)

####

a_betas_only_2 <- cbind(
  totasec_beta_matrix_2[,1],
  totasec_beta_matrix_2[,4],
  totasec_beta_matrix_2[,7]
)

colnames(a_betas_only_2) <- c("a_both_B_r",
                              "a_esc_B_r",
                              "a_nor_B_r")

a_betas_only_2 <- a_betas_only_2[1:30,]

a_SEs_only_2 <- cbind(
  totasec_beta_matrix_2[,2],
  totasec_beta_matrix_2[,5],
  totasec_beta_matrix_2[,8]
)

colnames(a_SEs_only_2) <- c("a_both_SE_r",
                            "a_esc_SE_r",
                            "a_nor_SE_r")

a_SEs_only_2 <- a_SEs_only_2[1:30,]

a_ps_only_2 <- cbind(
  totasec_beta_matrix_2[,3],
  totasec_beta_matrix_2[,6],
  totasec_beta_matrix_2[,9]
)

colnames(a_ps_only_2) <- c("a_both_p_r",
                           "a_esc_p_r",
                           "a_nor_p_r")

a_ps_only_2 <- a_ps_only_2[1:30,]

a_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                     c("both", "esc", "nor")))

a_lowest_ps_2[which(a_ps_only_2<=0.05)] <- 
  a_ps_only_2[which(a_ps_only_2<=0.05)]

heatmap.2( a_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(a_lowest_ps_2,5),notecol="black",
           main="Total ASEC, influential outliers removed"
)

##### 16) ANALYSE EFFECTS OF POLYGENIC SCORE ON TOTAL SIDE EFFECT BURDEN (TOTASEC): random only #####

##### 16a) Build models for totasec #####

# Generally all variables are included for theoretical reasons
# and because they have been included in previous studies.
# PCs 1-3 were important in screeplot when QC-ing SNP data.
# People report side effects more when zf1score (mood) is worse/higher
# - Uher et al (2009) "Adverse reactions to antidepressants".
# I did test using dose as a covariate but ultimately I 
# decided not to because, due to missing data, 66 people become excluded;
# also, Uher et al (2009) found it to be non-significant when time was
# taken into account (and weak when time wasn't)

# Below code verifies that zf1score is significant. I will keep the models
# the same as for the previous analysis (that includes non-randomized
# participants) for consistency, even though week and week2 perform diferently
# as covariates here

mod1 <- lmer( ztotasec ~ zf1score + week + week2 +
                sex + cage + drug + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

drop1(mod1)

mod2 <- lmer( ztotasec ~ zf1score + week + 
                sex + cage + drug + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand, REML=FALSE )

drop1(mod2) 

anova(mod1,mod2) # The models are not significantly different

hist(residuals(mod2)) 
predictmeans::residplot(mod2, newwd=FALSE)

mod1 <- lmer( ztotasec ~ zf1score + week + week2 +
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand_e, REML=FALSE )

drop1(mod1)

mod2 <- lmer( ztotasec ~ zf1score + week + 
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand_e, REML=FALSE )

drop1(mod2) 

hist(residuals(mod1)) 
predictmeans::residplot(mod1, newwd=FALSE)

mod1 <- lmer( ztotasec ~ zf1score + week + week2 +
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand_n, REML=FALSE )

drop1(mod1) 

mod2 <- lmer( ztotasec ~ zf1score + week + 
                sex + cage + SZ_0_0001 +
                zPC1 + zPC2 + zPC3 +
                (1|centreid) + (1|subjectid),
              data = GEND_rand_n, REML=FALSE )

drop1(mod2) 

hist(residuals(mod2)) 
predictmeans::residplot(mod2, newwd=FALSE)

##### 16b) Fit mixed effects models for totasec #####

# Make lists to store models in

a_models_r_both <- vector(mode = "list", length = length(PRSs_31))
a_models_r_esc <- vector(mode = "list", length = length(PRSs_31))
a_models_r_nor <- vector(mode = "list", length = length(PRSs_31))

a_models_r_both_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_r_esc_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_r_nor_2 <- vector(mode = "list", length = length(PRSs_31))

# Make new data set without missing data in important variables

GENDEP_a_r <- subset( GEND_rand, !is.na(totasec) )
GENDEP_a_r_e <- subset( GENDEP_a_r, drug==0 )
GENDEP_a_r_n <- subset( GENDEP_a_r, drug==1 )

# Make a vector to note models which failed to converge

ftc_r_totasec <- vector(mode = "list", length = 0)
ftc_r_totasec_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a <- lmer( ztotasec ~ zf1score + week + 
                       sex + cage + drug + GENDEP_a_r[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_a_r, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a)[-1]))==1){
    ftc_r_totasec[length(ftc_r_totasec)+1] <- lm_both_a
    names(ftc_r_totasec)[length(ftc_r_totasec)] <- PRSs[i]
  }
  a_models_r_both[[i]] <- lm_both_a
  
  removd <- romr.fnc(a_models_r_both[[i]], GENDEP_a_r, trim = 2.5)$data # remove influential outliers
  
  lm_both_a <- lmer( ztotasec ~ zf1score + week + 
                       sex + cage + drug + removd[,PRSs[i]] +
                       zPC1 + zPC2 + zPC3 +
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a)[-1]))==1){
    ftc_r_totasec_2[length(ftc_r_totasec_2)+1] <- lm_both_a
    names(ftc_r_totasec_2)[length(ftc_r_totasec_2)] <- PRSs[i]
  }
  a_models_r_both_2[[i]] <- lm_both_a
  
  
  lm_esc_a<- lmer( ztotasec ~ zf1score + week + week2 +
                     sex + cage + GENDEP_a_r_e[,PRSs[i]] +
                     zPC1 + zPC2 + zPC3 +
                     (1|centreid) + (1|subjectid),
                   data = GENDEP_a_r_e, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a)[-1]))==1){
    ftc_r_totasec[length(ftc_r_totasec)+1] <- lm_esc_a
    names(ftc_r_totasec)[length(ftc_r_totasec)] <- PRSs[i]
  }
  a_models_r_esc[[i]] <- lm_esc_a
  
  removd <- romr.fnc(a_models_r_esc[[i]], GENDEP_a_r_e, trim = 2.5)$data # remove influential outliers
  
  lm_esc_a <- lmer( ztotasec ~ zf1score + week + week2 +
                      sex + cage + GENDEP_a_r_e[,PRSs[i]] +
                      zPC1 + zPC2 + zPC3 +
                      (1|centreid) + (1|subjectid),
                    data = GENDEP_a_r_e, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a)[-1]))==1){
    ftc_r_totasec_2[length(ftc_r_totasec_2)+1] <- lm_esc_a
    names(ftc_r_totasec_2)[length(ftc_r_totasec_2)] <- PRSs[i]
  }
  a_models_r_esc_2[[i]] <- lm_esc_a
  
  
  lm_nor_a <- lmer( ztotasec ~ zf1score + week + 
                      sex + cage + GENDEP_a_r_n[,PRSs[i]] +
                      zPC1 + zPC2 + zPC3 +
                      (1|centreid) + (1|subjectid),
                    data = GENDEP_a_r_n, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a)[-1]))==1){
    ftc_r_totasec[length(ftc_r_totasec)+1] <- lm_nor_a
    names(ftc_r_totasec)[length(ftc_r_totasec)] <- PRSs[i]
  }
  a_models_r_nor[[i]] <- lm_nor_a
  
  removd <- romr.fnc(a_models_r_nor[[i]], GENDEP_a_r_n, trim = 2.5)$data # remove influential outliers
  
  lm_nor_a <- lmer( ztotasec ~ zf1score + week + 
                      sex + cage + removd[,PRSs[i]] +
                      zPC1 + zPC2 + zPC3 +
                      (1|centreid) + (1|subjectid),
                    data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a)[-1]))==1){
    ftc_r_totasec_2[length(ftc_r_totasec_2)+1] <- lm_nor_a
    names(ftc_r_totasec_2)[length(ftc_r_totasec_2)] <- PRSs[i]
  }
  a_models_r_nor_2[[i]] <- lm_nor_a
  
}

names(a_models_r_both) <- PRSs_31
names(a_models_r_esc) <- PRSs_31
names(a_models_r_nor) <- PRSs_31

names(a_models_r_both_2) <- PRSs_31
names(a_models_r_esc_2) <- PRSs_31
names(a_models_r_nor_2) <- PRSs_31

# Look at zmadrs models without polygenic score

noPRS_both_a <- lmer( ztotasec ~ zf1score + week + 
                        sex + cage + drug + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_a_r, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a)[-1]))==1){
  ftc_r_totasec[length(ftc_r_totasec)+1] <- noPRS_both_a
  names(ftc_r_totasec)[length(ftc_r_totasec)] <- "No_PRS"
}
a_models_r_both[[31]] <- noPRS_both_a

removd <- romr.fnc(a_models_r_both[[31]], GENDEP_a_r, trim = 2.5)$data # remove influential outliers

noPRS_both_a <- lmer( ztotasec ~ zf1score + week + 
                        sex + cage + drug + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a)[-1]))==1){
  ftc_r_totasec_2[length(ftc_r_totasec_2)+1] <- noPRS_both_a
  names(ftc_r_totasec_2)[length(ftc_r_totasec_2)] <- "No_PRS"
}
a_models_r_both_2[[31]] <- noPRS_both_a

noPRS_esc_a <- lmer( ztotasec ~ zf1score + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_a_r_e, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a)[-1]))==1){
  ftc_r_totasec[length(ftc_r_totasec)+1] <- noPRS_esc_a
  names(ftc_r_totasec)[length(ftc_r_totasec)] <- "No_PRS"
}
a_models_r_esc[[31]] <- noPRS_esc_a

removd <- romr.fnc(a_models_r_esc[[31]], GENDEP_a_r_e, trim = 2.5)$data # remove influential outliers

noPRS_esc_a <- lmer( ztotasec ~ zf1score + week + week2 +
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a)[-1]))==1){
  ftc_r_totasec_2[length(ftc_r_totasec_2)+1] <- noPRS_esc_a
  names(ftc_r_totasec_2)[length(ftc_r_totasec_2)] <- "No_PRS"
}
a_models_r_esc_2[[31]] <- noPRS_esc_a

noPRS_nor_a <- lmer( ztotasec ~ zf1score + week + 
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = GENDEP_a_r_n, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a)[-1]))==1){
  ftc_r_totasec[length(ftc_r_totasec)+1] <- noPRS_nor_a
  names(ftc_r_totasec)[length(ftc_r_totasec)] <- "No_PRS"
}
a_models_r_nor[[31]] <- noPRS_nor_a

removd <- romr.fnc(a_models_r_nor[[31]], GENDEP_a_r_n, trim = 2.5)$data # remove influential outliers

noPRS_nor_a <- lmer( ztotasec ~ zf1score + week + 
                       sex + cage + 
                       zPC1 + zPC2 + zPC3 + 
                       (1|centreid) + (1|subjectid),
                     data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a)[-1]))==1){
  ftc_r_totasec_2[length(ftc_r_totasec_2)+1] <- noPRS_nor_a
  names(ftc_r_totasec_2)[length(ftc_r_totasec_2)] <- "No_PRS"
}
a_models_r_nor_2[[31]] <- noPRS_nor_a

##### 16c) Get betas, standard errors and p-values for totasec #####

a_models_r_both_betas <- vector(mode = "list", length = length(PRSs_31))
a_models_r_esc_betas <- vector(mode = "list", length = length(PRSs_31))
a_models_r_nor_betas <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  a_models_r_both_betas[i] <- list(summary(a_models_r_both[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_r_esc_betas[i] <- list(summary(a_models_r_esc[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_r_nor_betas[i] <- list(summary(a_models_r_nor[[i]])$coefficients[6,])
}

names(a_models_r_both_betas) <- PRSs_31
names(a_models_r_esc_betas) <- PRSs_31
names(a_models_r_nor_betas) <- PRSs_31

totasec_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a_models_r_both_betas))[,1],
    t(as.data.frame(a_models_r_both_betas))[,2],
    t(as.data.frame(a_models_r_both_betas))[,5],
    t(as.data.frame(a_models_r_esc_betas))[,1],
    t(as.data.frame(a_models_r_esc_betas))[,2],
    t(as.data.frame(a_models_r_esc_betas))[,5],
    t(as.data.frame(a_models_r_nor_betas))[,1],
    t(as.data.frame(a_models_r_nor_betas))[,2],
    t(as.data.frame(a_models_r_nor_betas))[,5]
  )
)

colnames( totasec_beta_matrix ) <- c(
  "totasec_both_B", "totasec_both_SE","totasec_both_p",
  "totasec_esc_B", "totasec_esc_SE", "totasec_esc_p",
  "totasec_nor_B", "totasec_nor_SE", "totasec_nor_p"
)

####

a_models_r_both_betas_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_r_esc_betas_2 <- vector(mode = "list", length = length(PRSs_31))
a_models_r_nor_betas_2 <- vector(mode = "list", length = length(PRSs_31))

for (i in 1:length(PRSs_31)){
  a_models_r_both_betas_2[i] <- list(summary(a_models_r_both_2[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_r_esc_betas_2[i] <- list(summary(a_models_r_esc_2[[i]])$coefficients[7,])
}
for (i in 1:length(PRSs_31)){
  a_models_r_nor_betas_2[i] <- list(summary(a_models_r_nor_2[[i]])$coefficients[6,])
}

names(a_models_r_both_betas_2) <- PRSs_31
names(a_models_r_esc_betas_2) <- PRSs_31
names(a_models_r_nor_betas_2) <- PRSs_31

totasec_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(a_models_r_both_betas_2))[,1],
    t(as.data.frame(a_models_r_both_betas_2))[,2],
    t(as.data.frame(a_models_r_both_betas_2))[,5],
    t(as.data.frame(a_models_r_esc_betas_2))[,1],
    t(as.data.frame(a_models_r_esc_betas_2))[,2],
    t(as.data.frame(a_models_r_esc_betas_2))[,5],
    t(as.data.frame(a_models_r_nor_betas_2))[,1],
    t(as.data.frame(a_models_r_nor_betas_2))[,2],
    t(as.data.frame(a_models_r_nor_betas_2))[,5]
  )
)

colnames( totasec_beta_matrix_2 ) <- c(
  "totasec_both_B", "totasec_both_SE","totasec_both_p",
  "totasec_esc_B", "totasec_esc_SE", "totasec_esc_p",
  "totasec_nor_B", "totasec_nor_SE", "totasec_nor_p"
)

##### 16d) Get PRSs' explained variances for totasec #####

a_matrix_both <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

a_matrix_esc <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

a_matrix_nor <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  a_matrix_both[i,] <- r.squaredGLMM(a_models_r_both[[i]])
  a_matrix_esc[i,] <- r.squaredGLMM(a_models_r_esc[[i]])
  a_matrix_nor[i,] <- r.squaredGLMM(a_models_r_nor[[i]])
}


####

a_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),
  c("R2m","R2c")))

a_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

a_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,
  dimnames=list(c(
    PRSs_31
  ),c("R2m","R2c")))

for (i in 1:
     31){
  a_matrix_both_2[i,] <- r.squaredGLMM(a_models_r_both_2[[i]])
  a_matrix_esc_2[i,] <- r.squaredGLMM(a_models_r_esc_2[[i]])
  a_matrix_nor_2[i,] <- r.squaredGLMM(a_models_r_nor_2[[i]])
}

##### 16e) Make output files for totasec #####

sink( "totasec_rand_only.csv")

# Betas and p-values

cat("totasec_betas\n,")
write.table( totasec_beta_matrix[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(a_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a_matrix_nor,
            col.names=TRUE, sep=",")

sink()

####

sink( "totasec_rand_only_inf_outliers_removed.csv")

# Betas and p-values

cat("totasec_betas\n,")
write.table( totasec_beta_matrix_2[1:30,], col.names=TRUE, sep="," )

# Explained variance tables

cat("Both\n,")
write.table(a_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()

##### 16f) Make heatmap of effect sizes for totasec #####

a_betas_only <- cbind(
  totasec_beta_matrix[,1],
  totasec_beta_matrix[,4],
  totasec_beta_matrix[,7]
)

colnames(a_betas_only) <- c("a_both_B_r",
                            "a_esc_B_r",
                            "a_nor_B_r")

a_betas_only <- a_betas_only[1:30,]

a_SEs_only <- cbind(
  totasec_beta_matrix[,2],
  totasec_beta_matrix[,5],
  totasec_beta_matrix[,8]
)

colnames(a_SEs_only) <- c("a_both_SE_r",
                          "a_esc_SE_r",
                          "a_nor_SE_r")

a_SEs_only <- a_SEs_only[1:30,]

a_ps_only <- cbind(
  totasec_beta_matrix[,3],
  totasec_beta_matrix[,6],
  totasec_beta_matrix[,9]
)

colnames(a_ps_only) <- c("a_both_p_r",
                         "a_esc_p_r",
                         "a_nor_p_r")

a_ps_only <- a_ps_only[1:30,]

a_lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                   c("both", "esc", "nor")))

a_lowest_ps[which(a_ps_only<=0.05)] <- 
  a_ps_only[which(a_ps_only<=0.05)]

heatmap.2( a_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(a_lowest_ps,5),notecol="black",
           main="Total ASEC, randomized-only"
)

####

a_betas_only_2 <- cbind(
  totasec_beta_matrix_2[,1],
  totasec_beta_matrix_2[,4],
  totasec_beta_matrix_2[,7]
)

colnames(a_betas_only_2) <- c("a_both_B_r",
                              "a_esc_B_r",
                              "a_nor_B_r")

a_betas_only_2 <- a_betas_only_2[1:30,]

a_SEs_only_2 <- cbind(
  totasec_beta_matrix_2[,2],
  totasec_beta_matrix_2[,5],
  totasec_beta_matrix_2[,8]
)

colnames(a_SEs_only_2) <- c("a_both_SE_r",
                            "a_esc_SE_r",
                            "a_nor_SE_r")

a_SEs_only_2 <- a_SEs_only_2[1:30,]

a_ps_only_2 <- cbind(
  totasec_beta_matrix_2[,3],
  totasec_beta_matrix_2[,6],
  totasec_beta_matrix_2[,9]
)

colnames(a_ps_only_2) <- c("a_both_p_r",
                           "a_esc_p_r",
                           "a_nor_p_r")

a_ps_only_2 <- a_ps_only_2[1:30,]

a_lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,
                                                     c("both", "esc", "nor")))

a_lowest_ps_2[which(a_ps_only_2<=0.05)] <- 
  a_ps_only_2[which(a_ps_only_2<=0.05)]

heatmap.2( a_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(a_lowest_ps_2,5),notecol="black",
           main="Total ASEC, randomized-only, influential outliers removed"
)

##### ANALYSE EFFECTS OF POLYGENIC SCORE ON INDIVIDUAL SIDE EFFECTS (ASEC 1-21) #####

##### Optimise models for each side effect #####

## Start with week, week2, sex, cage, drug, 3x PCs, and
## random effects centre and subject

## I would like to exclude the F1 (mood) dimension, even though it
## is usually significant (Uher et al, 2009) because
#### a) some ASEC items are symptoms of anxiety, so correcting for
# mood (which includes anxiety items) seems wrong
#### b) some ASEC items are not significantly affected by mood
# (Uher et al, 2009)
#### c) I'm doing a separate analysis on the effects of PSs on mood
# so I can always re-run any significant findings for PSs that are
# correlated with mood, and include mood as a covariate
#### d) if a patient was more likely to feel that they have a 
# moderate or severe side effect, and this were to impair their experience
# of a drug, and their polygenic score is at play, we don't want to 
# miss this by "correcting for" mood - that was part of their bad day,
# and, in the real world, that effect can't be deducted

## Side effects may have different trajectories depending
## on drug, so optimise models separately

## Random effect for centre can be removed, as can other
## terms, if doing so helps prevent "fail to converge" or
## singular fit issues

## a1 #####

drop1(
  lmer( asec1wk ~ 
        CZ_0_0001 + 
        week + week2 +
        sex + cage + drug +
        zPC1 + zPC2 + zPC3 + 
        (1|centreid) + (1|subjectid),
      data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <- lmer( asec1wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1)

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec1wk ~ 
        CZ_0_0001 + 
        week + week2 +
        sex + cage + 
        zPC1 + zPC2 + zPC3 + 
        (1|centreid) + (1|subjectid),
      data = GENDEP_escit, REML=FALSE )
) # Remove week2 for escitalopram

drop1(
  lmer( asec1wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) 

drop1(
  lmer( asec1wk ~ 
        CZ_0_0001 + 
        week + week2 +
        sex + cage + 
        zPC1 + zPC2 + zPC3 + 
        (1|centreid) + (1|subjectid),
      data = GENDEP_nortrip, REML=FALSE )
) # Keep week2 for nortriptyline

## a2 #####

mod1 <-
  lmer( asec2wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE ) 

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec2wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec2wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) 

drop1(
  lmer( asec2wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Keep week2

## a3 #####

mod1 <-
  lmer( asec3wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec3wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec3wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) 

drop1(
  lmer( asec3wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec3wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
)

## a4 #####

drop1(
  lmer( asec4wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <-
  lmer( asec4wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) 

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec4wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Keep week2

drop1(
  lmer( asec4wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Keep week2

## a5 #####

mod1 <-
  lmer( asec5wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE ) 

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec5wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Keep week2

drop1(
  lmer( asec5wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Keep week2

## a6 #####

drop1(
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <- 
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1)

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Singular fit. Remove terms.

drop1(
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Singular fit. Remove terms.

drop1(
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week +
          sex + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Singular fit. Remove terms.

mod1 <-
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week +
          sex + 
          (1|subjectid),
        data = GENDEP_escit, REML=FALSE )

drop1(mod1) # Fine

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec6wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Keep week2

## a7 #####

mod1 <-
  lmer( asec7wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec7wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec7wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
)

drop1(
  lmer( asec7wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Keep week2

## a8 #####

drop1(
  lmer( asec8wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <-
  lmer( asec8wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) #  Week still not significant, but keep it in

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec8wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec8wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Week still not significant, but keep it in

drop1(
  lmer( asec8wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec8wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
)

## a9 #####

mod1 <-
  lmer( asec9wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec9wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Keep week2

drop1(
  lmer( asec9wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Keep week2

## a10 #####

mod1 <-
  lmer( asec10wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec10wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Keep week2

drop1(
  lmer( asec10wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec10wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
)

## a11 #####

drop1(
  lmer( asec11wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <-
  lmer( asec11wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Week still not significant, but keep in model

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec11wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec11wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Week still not significant, but keep in model

drop1(
  lmer( asec11wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec11wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
)  # Week still not significant, but keep in model

## a12 #####

drop1(
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <-
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Week still not significant, but keep in model

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE )
) # Week still not significant, but keep in model

drop1(
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
)  # Singular fit. Delete terms. PC1 is significant

drop1(
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 +  
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
)  # Singular fit. Delete terms. 

drop1(
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week +
          cage + 
          zPC1 +  
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Singular fit. Delete terms. 

mod1 <-
  lmer( asec12wk ~ 
          CZ_0_0001 + 
          week +
          cage + 
          zPC1 +  
          (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )

drop1(mod1) # Fine

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

## a13 #####

drop1(
  lmer( asec13wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <-
  lmer( asec13wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1)

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec13wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))
# Fine

drop1(
  lmer( asec13wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec13wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

## a14 #####

mod1 <-
  lmer( asec14wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec14wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))
# Fine

drop1(
  lmer( asec14wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec14wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

## a15 #####

mod1 <-
  lmer( asec15wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Keep week2

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec15wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

# Fine

drop1(
  lmer( asec15wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec15wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

## a16 #####

drop1(
  lmer( asec16wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <-
  lmer( asec16wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) 

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec16wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

# Remove week2

drop1(
  lmer( asec16wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

drop1(
  lmer( asec16wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec16wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

## a17 #####

drop1(
  lmer( asec17wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2, also it failed to converge

mod1 <-
  lmer( asec17wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) 

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec17wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

# Remove week2

drop1(
  lmer( asec17wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

drop1(
  lmer( asec17wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec17wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Week not significant but leave it in

## a18 #####

drop1(
  lmer( asec18wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2,

mod1 <-
  lmer( asec18wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) 

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec18wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

# Remove week2

drop1(
  lmer( asec18wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

drop1(
  lmer( asec18wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec18wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

## a19 #####

mod1 <-
  lmer( asec19wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1) # Fine

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec19wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

drop1(
  lmer( asec19wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec19wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

## a20 #####

drop1(
  lmer( asec20wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )
) # Remove week2

mod1 <-
  lmer( asec20wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1)

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec20wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

# Remove week2

drop1(
  lmer( asec20wk ~ 
          CZ_0_0001 + 
          week + 
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

drop1(
  lmer( asec20wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) # Remove week2

drop1(
  lmer( asec20wk ~ 
          CZ_0_0001 + 
          week +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

## a21 #####

mod1 <-
  lmer( asec21wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + drug +
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP, REML=FALSE )

drop1(mod1)

hist(residuals(mod1),breaks=20) 
predictmeans::residplot(mod1, newwd=FALSE)

drop1(
  lmer( asec21wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_escit, REML=FALSE ))

drop1(
  lmer( asec21wk ~ 
          CZ_0_0001 + 
          week + week2 +
          sex + cage + 
          zPC1 + zPC2 + zPC3 + 
          (1|centreid) + (1|subjectid),
        data = GENDEP_nortrip, REML=FALSE )
) 

##### Fit models for a1: dry mouth (cholinergic) #####

GENDEP_a1 <- subset( GENDEP, !is.na(asec1wk) )
GENDEP_a1_e <- subset( GENDEP_a1, drug==0 )
GENDEP_a1_n <- subset( GENDEP_a1, drug==1 )

a1_models_both <- vector(mode = "list", length = length(PRSs_31))
a1_models_esc <- vector(mode = "list", length = length(PRSs_31))
a1_models_nor <- vector(mode = "list", length = length(PRSs_31))

a1_models_both_2 <- vector(mode = "list", length = length(PRSs_31))
a1_models_esc_2 <- vector(mode = "list", length = length(PRSs_31))
a1_models_nor_2 <- vector(mode = "list", length = length(PRSs_31))

ftc_a1 <- vector(mode = "list", length = 0)
ftc_a1_2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a1_s <- lmer( asec1wk ~ week + 
                          sex + cage + drug + GENDEP_a1[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_a1, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a1_s)[-1]))==1){
    ftc_a1[length(ftc_a1)+1] <- lm_both_a1_s
    names(ftc_a1)[length(ftc_a1)] <- PRSs[i]
  }
  a1_models_both[[i]] <- lm_both_a1_s
  
  removd <- romr.fnc(a1_models_both[[i]], GENDEP_a1, trim = 2.5)$data # remove influential outliers
  
  lm_both_a1_s <- lmer( asec1wk ~ week + 
                          sex + cage + drug + removd[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a1_s)[-1]))==1){
    ftc_a1_2[length(ftc_a1_2)+1] <- lm_both_a1_s
    names(ftc_a1_2)[length(ftc_a1_2)] <- PRSs[i]
  }
  a1_models_both_2[[i]] <- lm_both_a1_s
  
  ##
  
  lm_esc_a1_s <- lmer( asec1wk ~ week + 
                         sex + cage + GENDEP_a1_e[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_a1_e, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a1_s)[-1]))==1){
    ftc_a1[length(ftc_a1)+1] <- lm_esc_a1_s
    names(ftc_a1)[length(ftc_a1)] <- PRSs[i]
  }
  a1_models_esc[[i]] <- lm_esc_a1_s
  
  removd <- romr.fnc(a1_models_esc[[i]], GENDEP_a1_e, trim = 2.5)$data # remove influential outliers
  
  lm_esc_a1_s <- lmer( asec1wk ~ week + 
                         sex + cage + removd[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a1_s)[-1]))==1){
    ftc_a1_2[length(ftc_a1_2)+1] <- lm_esc_a1_s
    names(ftc_a1_2)[length(ftc_a1_2)] <- PRSs[i]
  }
  a1_models_esc_2[[i]] <- lm_esc_a1_s
  
  ##
  
  lm_nor_a1_s <- lmer( asec1wk ~ week + week2 +
                         sex + cage + GENDEP_a1_n[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_a1_n, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a1_s)[-1]))==1){
    ftc_a1[length(ftc_a1)+1] <- lm_nor_a1_s
    names(ftc_a1)[length(ftc_a1)] <- PRSs[i]
  }
  a1_models_nor[[i]] <- lm_nor_a1_s
  
  removd <- romr.fnc(a1_models_nor[[i]], GENDEP_a1_n, trim = 2.5)$data # remove influential outliers
  
  lm_nor_a1_s <- lmer( asec1wk ~ week + week2 +
                         sex + cage + removd[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = removd, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a1_s)[-1]))==1){
    ftc_a1_2[length(ftc_a1_2)+1] <- lm_nor_a1_s
    names(ftc_a1_2)[length(ftc_a1_2)] <- PRSs[i]
  }
  a1_models_nor_2[[i]] <- lm_nor_a1_s
  
}  

noPRS_both_a1 <- lmer( asec1wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_a1, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a1)[-1]))==1){
  ftc_a1[length(ftc_a1)+1] <- noPRS_both_a1
  names(ftc_a1)[length(ftc_a1)] <- "No_PRS"
}
a1_models_both[[31]] <- noPRS_both_a1

removd <- romr.fnc(a1_models_both[[31]], GENDEP_a1, trim = 2.5)$data

noPRS_both_a1 <- lmer( asec1wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a1)[-1]))==1){
  ftc_a1_2[length(ftc_a1_2)+1] <- noPRS_both_a1
  names(ftc_a1_2)[length(ftc_a1_2)] <- "No_PRS"
}
a1_models_both_2[[31]] <- noPRS_both_a1

##

noPRS_esc_a1 <- lmer( asec1wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_a1_e, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a1)[-1]))==1){
  ftc_a1[length(ftc_a1)+1] <- noPRS_esc_a1
  names(ftc_a1)[length(ftc_a1)] <- "No_PRS"
}
a1_models_esc[[31]] <- noPRS_esc_a1

removd <- romr.fnc(a1_models_esc[[31]], GENDEP_a1_e, trim = 2.5)$data

noPRS_esc_a1 <- lmer( asec1wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a1)[-1]))==1){
  ftc_a1_2[length(ftc_a1_2)+1] <- noPRS_esc_a1
  names(ftc_a1_2)[length(ftc_a1_2)] <- "No_PRS"
}
a1_models_esc_2[[31]] <- noPRS_esc_a1

##

noPRS_nor_a1 <- lmer( asec1wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_a1_n, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a1)[-1]))==1){
  ftc_a1[length(ftc_a1)+1] <- noPRS_nor_a1
  names(ftc_a1)[length(ftc_a1)] <- "No_PRS"
}
a1_models_nor[[31]] <- noPRS_nor_a1

removd <- romr.fnc(a1_models_nor[[31]], GENDEP_a1_n, trim = 2.5)$data

noPRS_nor_a1 <- lmer( asec1wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = removd, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a1)[-1]))==1){
  ftc_a1_2[length(ftc_a1_2)+1] <- noPRS_nor_a1
  names(ftc_a1_2)[length(ftc_a1_2)] <- "No_PRS"
}
a1_models_nor_2[[31]] <- noPRS_nor_a1

names(a1_models_both) <- PRSs_31
names(a1_models_esc) <- PRSs_31
names(a1_models_nor) <- PRSs_31

names(a1_models_both_2) <- PRSs_31
names(a1_models_esc_2) <- PRSs_31
names(a1_models_nor_2) <- PRSs_31

##### Get betas, standard errors, and p-values for a1 #####

a1_models_both_betas <- vector(mode = "list", length = length(PRSs))
a1_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a1_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a1_models_both_betas[i] <- list(summary(a1_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a1_models_esc_betas[i] <- list(summary(a1_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a1_models_nor_betas[i] <- list(summary(a1_models_nor[[i]])$coefficients[6,])
}

names(a1_models_both_betas) <- PRSs
names(a1_models_esc_betas) <- PRSs
names(a1_models_nor_betas) <- PRSs

a1_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a1_models_both_betas))[,1],
    t(as.data.frame(a1_models_both_betas))[,2],
    t(as.data.frame(a1_models_both_betas))[,5],
    t(as.data.frame(a1_models_esc_betas))[,1],
    t(as.data.frame(a1_models_esc_betas))[,2],
    t(as.data.frame(a1_models_esc_betas))[,5],
    t(as.data.frame(a1_models_nor_betas))[,1],
    t(as.data.frame(a1_models_nor_betas))[,2],
    t(as.data.frame(a1_models_nor_betas))[,5]
  )
)

colnames( a1_beta_matrix ) <- c(
  "a1_both_B", "a1_both_se", "a1_both_p", 
  "a1_esc_B", "a1_esc_se", "a1_esc_p",
  "a1_nor_B", "a1_nor_se", "a1_nor_p"
)


a1_models_both_betas_2 <- vector(mode = "list", length = length(PRSs))
a1_models_esc_betas_2 <- vector(mode = "list", length = length(PRSs))
a1_models_nor_betas_2 <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a1_models_both_betas_2[i] <- list(summary(a1_models_both_2[[i]])$coefficients[6,])
}
for (i in 1:30){
  a1_models_esc_betas_2[i] <- list(summary(a1_models_esc_2[[i]])$coefficients[5,])
}
for (i in 1:30){
  a1_models_nor_betas_2[i] <- list(summary(a1_models_nor_2[[i]])$coefficients[6,])
}

names(a1_models_both_betas_2) <- PRSs
names(a1_models_esc_betas_2) <- PRSs
names(a1_models_nor_betas_2) <- PRSs

a1_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(a1_models_both_betas_2))[,1],
    t(as.data.frame(a1_models_both_betas_2))[,2],
    t(as.data.frame(a1_models_both_betas_2))[,5],
    t(as.data.frame(a1_models_esc_betas_2))[,1],
    t(as.data.frame(a1_models_esc_betas_2))[,2],
    t(as.data.frame(a1_models_esc_betas_2))[,5],
    t(as.data.frame(a1_models_nor_betas_2))[,1],
    t(as.data.frame(a1_models_nor_betas_2))[,2],
    t(as.data.frame(a1_models_nor_betas_2))[,5]
  )
)

colnames( a1_beta_matrix_2 ) <- c(
  "a1_both_B", "a1_both_se", "a1_both_p", 
  "a1_esc_B", "a1_esc_se", "a1_esc_p",
  "a1_nor_B", "a1_nor_se", "a1_nor_p"
)

##### Get PRSs' explained variance for a1 #####

a1_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a1_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a1_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a1_matrix_both[i,] <- r.squaredGLMM(a1_models_both[[i]])
  a1_matrix_esc[i,] <- r.squaredGLMM(a1_models_esc[[i]])
  a1_matrix_nor[i,] <- r.squaredGLMM(a1_models_nor[[i]])
}

####

a1_matrix_both_2 <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a1_matrix_esc_2 <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a1_matrix_nor_2 <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a1_matrix_both_2[i,] <- r.squaredGLMM(a1_models_both_2[[i]])
  a1_matrix_esc_2[i,] <- r.squaredGLMM(a1_models_esc_2[[i]])
  a1_matrix_nor_2[i,] <- r.squaredGLMM(a1_models_nor_2[[i]])
}

##### Make output files for a1 #####

sink( "a1.csv")

# Betas and p-values

cat("a1_betas\n,")
write.table( a1_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a1_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a1_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a1_matrix_nor,
            col.names=TRUE, sep=",")

sink()

####


sink( "a1_inf_out_removed.csv")

# Betas and p-values

cat("a1_betas\n,")
write.table( a1_beta_matrix_2, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a1_matrix_both_2,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a1_matrix_esc_2,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a1_matrix_nor_2,
            col.names=TRUE, sep=",")

sink()

##### Make a1 heatmap #####

a1_models_betas <- c("a1_both_B", "a1_esc_B", "a1_nor_B")
a1_models_ses <- c("a1_both_se", "a1_esc_se", "a1_nor_se")
a1_models_ps <- c("a1_both_p", "a1_esc_p", "a1_nor_p")

a1_models_betas_columns <- which( colnames(a1_beta_matrix)
                                  %in% a1_models_betas )
a1_models_ses_columns <- which( colnames(a1_beta_matrix)
                                %in% a1_models_ses )
a1_models_ps_columns <- which( colnames(a1_beta_matrix)
                               %in% a1_models_ps )

a1_betas_only <- cbind(a1_beta_matrix[,a1_models_betas_columns[1]],
                       a1_beta_matrix[,a1_models_betas_columns[2]],
                       a1_beta_matrix[,a1_models_betas_columns[3]])

colnames(a1_betas_only) <- c("a1_both", "a1_esc","a1_nor")

a1_ses_only <- cbind(a1_beta_matrix[,a1_models_ses_columns[1]],
                     a1_beta_matrix[,a1_models_ses_columns[2]],
                     a1_beta_matrix[,a1_models_ses_columns[3]])

colnames(a1_betas_only) <- c("a1_both", "a1_esc","a1_nor")

a1_p_only <- cbind(a1_beta_matrix[,a1_models_ps_columns[1]],
                   a1_beta_matrix[,a1_models_ps_columns[2]],
                   a1_beta_matrix[,a1_models_ps_columns[3]])

colnames(a1_p_only) <- c("a1_both", "a1_esc","a1_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a1_both",
  "a1_esc",
  "a1_nor")))

lowest_ps[which(a1_p_only<=0.05)]<-a1_p_only[which(a1_p_only<=0.05)]

lowest_ps_a1 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a1_both",
                                                            "a1_esc",
                                                            "a1_nor")))

lowest_ps_a1[which(a1_p_only<=0.05)]<-a1_p_only[which(a1_p_only<=0.05)]

heatmap.2( a1_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a1,5),notecol="black",
           main="Dry mouth"
)

####


a1_models_betas_2 <- c("a1_both_B", "a1_esc_B", "a1_nor_B")
a1_models_ses_2 <- c("a1_both_se", "a1_esc_se", "a1_nor_se")
a1_models_ps_2 <- c("a1_both_p", "a1_esc_p", "a1_nor_p")

a1_models_betas_2_columns <- which( colnames(a1_beta_matrix_2)
                                  %in% a1_models_betas_2 )
a1_models_ses_2_columns <- which( colnames(a1_beta_matrix_2)
                                %in% a1_models_ses_2 )
a1_models_ps_2_columns <- which( colnames(a1_beta_matrix_2)
                               %in% a1_models_ps_2 )

a1_betas_2_only_2 <- cbind(a1_beta_matrix_2[,a1_models_betas_2_columns[1]],
                       a1_beta_matrix_2[,a1_models_betas_2_columns[2]],
                       a1_beta_matrix_2[,a1_models_betas_2_columns[3]])

colnames(a1_betas_2_only_2) <- c("a1_both", "a1_esc","a1_nor")

a1_ses_2_only_2 <- cbind(a1_beta_matrix_2[,a1_models_ses_2_columns[1]],
                     a1_beta_matrix_2[,a1_models_ses_2_columns[2]],
                     a1_beta_matrix_2[,a1_models_ses_2_columns[3]])

colnames(a1_betas_2_only_2) <- c("a1_both", "a1_esc","a1_nor")

a1_p_only_2 <- cbind(a1_beta_matrix_2[,a1_models_ps_2_columns[1]],
                   a1_beta_matrix_2[,a1_models_ps_2_columns[2]],
                   a1_beta_matrix_2[,a1_models_ps_2_columns[3]])

colnames(a1_p_only_2) <- c("a1_both", "a1_esc","a1_nor")

lowest_ps_2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a1_both",
  "a1_esc",
  "a1_nor")))

lowest_ps_2[which(a1_p_only_2<=0.05)]<-a1_p_only_2[which(a1_p_only_2<=0.05)]

lowest_ps_2_a1 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a1_both",
                                                            "a1_esc",
                                                            "a1_nor")))

lowest_ps_2_a1[which(a1_p_only_2<=0.05)]<-a1_p_only_2[which(a1_p_only_2<=0.05)]

heatmap.2( a1_betas_2_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_2_a1,5),notecol="black",
           main="Dry mouth, inf out rem"
)

##### Fit models for a2: drowsiness (histaminergic) #####

a2_models_both <- vector(mode = "list", length = length(PRSs_31))
a2_models_esc <- vector(mode = "list", length = length(PRSs_31))
a2_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a2 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a2_s <- lmer( asec2wk ~ week + week2 +
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a2_s)[-1]))==1){
    ftc_a2[length(ftc_a2)+1] <- lm_both_a2_s
    names(ftc_a2)[length(ftc_a2)] <- PRSs[i]
  }
  a2_models_both[[i]] <- lm_both_a2_s

  lm_esc_a2_s <- lmer( asec2wk ~ week + 
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a2_s)[-1]))==1){
    ftc_a2[length(ftc_a2)+1] <- lm_esc_a2_s
    names(ftc_a2)[length(ftc_a2)] <- PRSs[i]
  }
  a2_models_esc[[i]] <- lm_esc_a2_s

  lm_nor_a2_s <- lmer( asec2wk ~ week + week2 +
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a2_s)[-1]))==1){
    ftc_a2[length(ftc_a2)+1] <- lm_nor_a2_s
    names(ftc_a2)[length(ftc_a2)] <- PRSs[i]
  }
  a2_models_nor[[i]] <- lm_nor_a2_s
}  

noPRS_both_a2 <- lmer( asec2wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a2)[-1]))==1){
  ftc_a2[length(ftc_a2)+1] <- noPRS_both_a2
  names(ftc_a2)[length(ftc_a2)] <- "No_PRS"
}
a2_models_both[[31]] <- noPRS_both_a2

noPRS_esc_a2 <- lmer( asec2wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a2)[-1]))==1){
  ftc_a2[length(ftc_a2)+1] <- noPRS_esc_a2
  names(ftc_a2)[length(ftc_a2)] <- "No_PRS"
}
a2_models_esc[[31]] <- noPRS_esc_a2

noPRS_nor_a2 <- lmer( asec2wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a2)[-1]))==1){
  ftc_a2[length(ftc_a2)+1] <- noPRS_nor_a2
  names(ftc_a2)[length(ftc_a2)] <- "No_PRS"
}
a2_models_nor[[31]] <- noPRS_nor_a2

names(a2_models_both) <- PRSs_31
names(a2_models_esc) <- PRSs_31
names(a2_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a2 #####

a2_models_both_betas <- vector(mode = "list", length = length(PRSs))
a2_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a2_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a2_models_both_betas[i] <- list(summary(a2_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a2_models_esc_betas[i] <- list(summary(a2_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a2_models_nor_betas[i] <- list(summary(a2_models_nor[[i]])$coefficients[6,])
}

names(a2_models_both_betas) <- PRSs
names(a2_models_esc_betas) <- PRSs
names(a2_models_nor_betas) <- PRSs

a2_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a2_models_both_betas))[,1],
    t(as.data.frame(a2_models_both_betas))[,2],
    t(as.data.frame(a2_models_both_betas))[,5],
    t(as.data.frame(a2_models_esc_betas))[,1],
    t(as.data.frame(a2_models_esc_betas))[,2],
    t(as.data.frame(a2_models_esc_betas))[,5],
    t(as.data.frame(a2_models_nor_betas))[,1],
    t(as.data.frame(a2_models_nor_betas))[,2],
    t(as.data.frame(a2_models_nor_betas))[,5]
  )
)

colnames( a2_beta_matrix ) <- c(
  "a2_both_B", "a2_both_se", "a2_both_p", 
  "a2_esc_B", "a2_esc_se", "a2_esc_p",
  "a2_nor_B", "a2_nor_se", "a2_nor_p"
)

##### Get PRSs' explained variance for a2 #####

a2_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a2_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a2_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a2_matrix_both[i,] <- r.squaredGLMM(a2_models_both[[i]])
  a2_matrix_esc[i,] <- r.squaredGLMM(a2_models_esc[[i]])
  a2_matrix_nor[i,] <- r.squaredGLMM(a2_models_nor[[i]])
}

##### Make output files for a2 #####

sink( "a2.csv")

# Betas and p-values

cat("a2_betas\n,")
write.table( a2_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a2_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a2_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a2_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a2 heatmap #####

a2_models_betas <- c("a2_both_B", "a2_esc_B", "a2_nor_B")
a2_models_ses <- c("a2_both_se", "a2_esc_se", "a2_nor_se")
a2_models_ps <- c("a2_both_p", "a2_esc_p", "a2_nor_p")

a2_models_betas_columns <- which( colnames(a2_beta_matrix)
                                  %in% a2_models_betas )
a2_models_ses_columns <- which( colnames(a2_beta_matrix)
                                %in% a2_models_ses )
a2_models_ps_columns <- which( colnames(a2_beta_matrix)
                               %in% a2_models_ps )

a2_betas_only <- cbind(a2_beta_matrix[,a2_models_betas_columns[1]],
                       a2_beta_matrix[,a2_models_betas_columns[2]],
                       a2_beta_matrix[,a2_models_betas_columns[3]])

colnames(a2_betas_only) <- c("a2_both", "a2_esc","a2_nor")

a2_ses_only <- cbind(a2_beta_matrix[,a2_models_ses_columns[1]],
                     a2_beta_matrix[,a2_models_ses_columns[2]],
                     a2_beta_matrix[,a2_models_ses_columns[3]])

colnames(a2_betas_only) <- c("a2_both", "a2_esc","a2_nor")

a2_p_only <- cbind(a2_beta_matrix[,a2_models_ps_columns[1]],
                   a2_beta_matrix[,a2_models_ps_columns[2]],
                   a2_beta_matrix[,a2_models_ps_columns[3]])

colnames(a2_p_only) <- c("a2_both", "a2_esc","a2_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a2_both",
  "a2_esc",
  "a2_nor")))

lowest_ps[which(a2_p_only<=0.05)]<-a2_p_only[which(a2_p_only<=0.05)]

lowest_ps_a2 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a2_both",
                                                            "a2_esc",
                                                            "a2_nor")))

lowest_ps_a2[which(a2_p_only<=0.05)]<-a2_p_only[which(a2_p_only<=0.05)]

heatmap.2( a2_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a2,5),notecol="black",
           main="Drowsiness"
)

##### Fit models for a3: insomnia (serotonergic) #####

a3_models_both <- vector(mode = "list", length = length(PRSs_31))
a3_models_esc <- vector(mode = "list", length = length(PRSs_31))
a3_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a3 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a3_s <- lmer( asec3wk ~ week + week2 +
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a3_s)[-1]))==1){
    ftc_a3[length(ftc_a3)+1] <- lm_both_a3_s
    names(ftc_a3)[length(ftc_a3)] <- PRSs[i]
  }
  a3_models_both[[i]] <- lm_both_a3_s
  
  lm_esc_a3_s <- lmer( asec3wk ~ week + 
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a3_s)[-1]))==1){
    ftc_a3[length(ftc_a3)+1] <- lm_esc_a3_s
    names(ftc_a3)[length(ftc_a3)] <- PRSs[i]
  }
  a3_models_esc[[i]] <- lm_esc_a3_s
  
  lm_nor_a3_s <- lmer( asec3wk ~ week + 
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a3_s)[-1]))==1){
    ftc_a3[length(ftc_a3)+1] <- lm_nor_a3_s
    names(ftc_a3)[length(ftc_a3)] <- PRSs[i]
  }
  a3_models_nor[[i]] <- lm_nor_a3_s
}  

noPRS_both_a3 <- lmer( asec3wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a3)[-1]))==1){
  ftc_a3[length(ftc_a3)+1] <- noPRS_both_a3
  names(ftc_a3)[length(ftc_a3)] <- "No_PRS"
}
a3_models_both[[31]] <- noPRS_both_a3

noPRS_esc_a3 <- lmer( asec3wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a3)[-1]))==1){
  ftc_a3[length(ftc_a3)+1] <- noPRS_esc_a3
  names(ftc_a3)[length(ftc_a3)] <- "No_PRS"
}
a3_models_esc[[31]] <- noPRS_esc_a3


noPRS_nor_a3 <- lmer( asec3wk ~ week + 
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a3)[-1]))==1){
  ftc_a3[length(ftc_a3)+1] <- noPRS_nor_a3
  names(ftc_a3)[length(ftc_a3)] <- "No_PRS"
}
a3_models_nor[[31]] <- noPRS_nor_a3

names(a3_models_both) <- PRSs_31
names(a3_models_esc) <- PRSs_31
names(a3_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a3 #####

a3_models_both_betas <- vector(mode = "list", length = length(PRSs))
a3_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a3_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a3_models_both_betas[i] <- list(summary(a3_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a3_models_esc_betas[i] <- list(summary(a3_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a3_models_nor_betas[i] <- list(summary(a3_models_nor[[i]])$coefficients[5,])
}

names(a3_models_both_betas) <- PRSs
names(a3_models_esc_betas) <- PRSs
names(a3_models_nor_betas) <- PRSs

a3_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a3_models_both_betas))[,1],
    t(as.data.frame(a3_models_both_betas))[,2],
    t(as.data.frame(a3_models_both_betas))[,5],
    t(as.data.frame(a3_models_esc_betas))[,1],
    t(as.data.frame(a3_models_esc_betas))[,2],
    t(as.data.frame(a3_models_esc_betas))[,5],
    t(as.data.frame(a3_models_nor_betas))[,1],
    t(as.data.frame(a3_models_nor_betas))[,2],
    t(as.data.frame(a3_models_nor_betas))[,5]
  )
)

colnames( a3_beta_matrix ) <- c(
  "a3_both_B", "a3_both_se", "a3_both_p", 
  "a3_esc_B", "a3_esc_se", "a3_esc_p",
  "a3_nor_B", "a3_nor_se", "a3_nor_p"
)

##### Get PRSs' explained variance for a3 #####

a3_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a3_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a3_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a3_matrix_both[i,] <- r.squaredGLMM(a3_models_both[[i]])
  a3_matrix_esc[i,] <- r.squaredGLMM(a3_models_esc[[i]])
  a3_matrix_nor[i,] <- r.squaredGLMM(a3_models_nor[[i]])
}

##### Make output files for a3 #####

sink( "a3.csv")

# Betas and p-values

cat("a3_betas\n,")
write.table( a3_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a3_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a3_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a3_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a3 heatmap #####

a3_models_betas <- c("a3_both_B", "a3_esc_B", "a3_nor_B")
a3_models_ses <- c("a3_both_se", "a3_esc_se", "a3_nor_se")
a3_models_ps <- c("a3_both_p", "a3_esc_p", "a3_nor_p")

a3_models_betas_columns <- which( colnames(a3_beta_matrix)
                                  %in% a3_models_betas )
a3_models_ses_columns <- which( colnames(a3_beta_matrix)
                                %in% a3_models_ses )
a3_models_ps_columns <- which( colnames(a3_beta_matrix)
                               %in% a3_models_ps )

a3_betas_only <- cbind(a3_beta_matrix[,a3_models_betas_columns[1]],
                       a3_beta_matrix[,a3_models_betas_columns[2]],
                       a3_beta_matrix[,a3_models_betas_columns[3]])

colnames(a3_betas_only) <- c("a3_both", "a3_esc","a3_nor")

a3_ses_only <- cbind(a3_beta_matrix[,a3_models_ses_columns[1]],
                     a3_beta_matrix[,a3_models_ses_columns[2]],
                     a3_beta_matrix[,a3_models_ses_columns[3]])

colnames(a3_betas_only) <- c("a3_both", "a3_esc","a3_nor")

a3_p_only <- cbind(a3_beta_matrix[,a3_models_ps_columns[1]],
                   a3_beta_matrix[,a3_models_ps_columns[2]],
                   a3_beta_matrix[,a3_models_ps_columns[3]])

colnames(a3_p_only) <- c("a3_both", "a3_esc","a3_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a3_both",
  "a3_esc",
  "a3_nor")))

lowest_ps[which(a3_p_only<=0.05)]<-a3_p_only[which(a3_p_only<=0.05)]

lowest_ps_a3 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a3_both",
                                                            "a3_esc",
                                                            "a3_nor")))

lowest_ps_a3[which(a3_p_only<=0.05)]<-a3_p_only[which(a3_p_only<=0.05)]

heatmap.2( a3_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a3,5),notecol="black",
           main="Insomnia"
)

##### Fit models for a4: blurred vision (cholinergic) #####

a4_models_both <- vector(mode = "list", length = length(PRSs_31))
a4_models_esc <- vector(mode = "list", length = length(PRSs_31))
a4_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a4 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a4_s <- lmer( asec4wk ~ week + 
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a4_s)[-1]))==1){
    ftc_a4[length(ftc_a4)+1] <- lm_both_a4_s
    names(ftc_a4)[length(ftc_a4)] <- PRSs[i]
  }
  a4_models_both[[i]] <- lm_both_a4_s
  
  lm_esc_a4_s <- lmer( asec4wk ~ week + week2 +
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a4_s)[-1]))==1){
    ftc_a4[length(ftc_a4)+1] <- lm_esc_a4_s
    names(ftc_a4)[length(ftc_a4)] <- PRSs[i]
  }
  a4_models_esc[[i]] <- lm_esc_a4_s
  
  lm_nor_a4_s <- lmer( asec4wk ~ week + week2 +
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a4_s)[-1]))==1){
    ftc_a4[length(ftc_a4)+1] <- lm_nor_a4_s
    names(ftc_a4)[length(ftc_a4)] <- PRSs[i]
  }
  a4_models_nor[[i]] <- lm_nor_a4_s
}  

noPRS_both_a4 <- lmer( asec4wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a4)[-1]))==1){
  ftc_a4[length(ftc_a4)+1] <- noPRS_both_a4
  names(ftc_a4)[length(ftc_a4)] <- "No_PRS"
}
a4_models_both[[31]] <- noPRS_both_a4

noPRS_esc_a4 <- lmer( asec4wk ~ week + week2 +
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a4)[-1]))==1){
  ftc_a4[length(ftc_a4)+1] <- noPRS_esc_a4
  names(ftc_a4)[length(ftc_a4)] <- "No_PRS"
}
a4_models_esc[[31]] <- noPRS_esc_a4

noPRS_nor_a4 <- lmer( asec4wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a4)[-1]))==1){
  ftc_a4[length(ftc_a4)+1] <- noPRS_nor_a4
  names(ftc_a4)[length(ftc_a4)] <- "No_PRS"
}
a4_models_nor[[31]] <- noPRS_nor_a4

names(a4_models_both) <- PRSs_31
names(a4_models_esc) <- PRSs_31
names(a4_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a4 #####

a4_models_both_betas <- vector(mode = "list", length = length(PRSs))
a4_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a4_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a4_models_both_betas[i] <- list(summary(a4_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a4_models_esc_betas[i] <- list(summary(a4_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a4_models_nor_betas[i] <- list(summary(a4_models_nor[[i]])$coefficients[6,])
}

names(a4_models_both_betas) <- PRSs
names(a4_models_esc_betas) <- PRSs
names(a4_models_nor_betas) <- PRSs

a4_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a4_models_both_betas))[,1],
    t(as.data.frame(a4_models_both_betas))[,2],
    t(as.data.frame(a4_models_both_betas))[,5],
    t(as.data.frame(a4_models_esc_betas))[,1],
    t(as.data.frame(a4_models_esc_betas))[,2],
    t(as.data.frame(a4_models_esc_betas))[,5],
    t(as.data.frame(a4_models_nor_betas))[,1],
    t(as.data.frame(a4_models_nor_betas))[,2],
    t(as.data.frame(a4_models_nor_betas))[,5]
  )
)

colnames( a4_beta_matrix ) <- c(
  "a4_both_B", "a4_both_se", "a4_both_p", 
  "a4_esc_B", "a4_esc_se", "a4_esc_p",
  "a4_nor_B", "a4_nor_se", "a4_nor_p"
)

##### Get PRSs' explained variance for a4 #####

a4_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a4_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a4_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a4_matrix_both[i,] <- r.squaredGLMM(a4_models_both[[i]])
  a4_matrix_esc[i,] <- r.squaredGLMM(a4_models_esc[[i]])
  a4_matrix_nor[i,] <- r.squaredGLMM(a4_models_nor[[i]])
}

##### Make output files for a4 #####

sink( "a4.csv")

# Betas and p-values

cat("a4_betas\n,")
write.table( a4_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a4_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a4_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a4_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a4 heatmap #####

a4_models_betas <- c("a4_both_B", "a4_esc_B", "a4_nor_B")
a4_models_ses <- c("a4_both_se", "a4_esc_se", "a4_nor_se")
a4_models_ps <- c("a4_both_p", "a4_esc_p", "a4_nor_p")

a4_models_betas_columns <- which( colnames(a4_beta_matrix)
                                  %in% a4_models_betas )
a4_models_ses_columns <- which( colnames(a4_beta_matrix)
                                %in% a4_models_ses )
a4_models_ps_columns <- which( colnames(a4_beta_matrix)
                               %in% a4_models_ps )

a4_betas_only <- cbind(a4_beta_matrix[,a4_models_betas_columns[1]],
                       a4_beta_matrix[,a4_models_betas_columns[2]],
                       a4_beta_matrix[,a4_models_betas_columns[3]])

colnames(a4_betas_only) <- c("a4_both", "a4_esc","a4_nor")

a4_ses_only <- cbind(a4_beta_matrix[,a4_models_ses_columns[1]],
                     a4_beta_matrix[,a4_models_ses_columns[2]],
                     a4_beta_matrix[,a4_models_ses_columns[3]])

colnames(a4_betas_only) <- c("a4_both", "a4_esc","a4_nor")

a4_p_only <- cbind(a4_beta_matrix[,a4_models_ps_columns[1]],
                   a4_beta_matrix[,a4_models_ps_columns[2]],
                   a4_beta_matrix[,a4_models_ps_columns[3]])

colnames(a4_p_only) <- c("a4_both", "a4_esc","a4_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a4_both",
  "a4_esc",
  "a4_nor")))

lowest_ps[which(a4_p_only<=0.05)]<-a4_p_only[which(a4_p_only<=0.05)]

lowest_ps_a4 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a4_both",
                                                            "a4_esc",
                                                            "a4_nor")))

lowest_ps_a4[which(a4_p_only<=0.05)]<-a4_p_only[which(a4_p_only<=0.05)]

heatmap.2( a4_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a4,5),notecol="black",
           main="Blurred vision"
)

##### Fit models for a5: headache #####

a5_models_both <- vector(mode = "list", length = length(PRSs_31))
a5_models_esc <- vector(mode = "list", length = length(PRSs_31))
a5_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a5 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a5_s <- lmer( asec5wk ~ week + week2 +
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a5_s)[-1]))==1){
    ftc_a5[length(ftc_a5)+1] <- lm_both_a5_s
    names(ftc_a5)[length(ftc_a5)] <- PRSs[i]
  }
  a5_models_both[[i]] <- lm_both_a5_s
  
  lm_esc_a5_s <- lmer( asec5wk ~ week + week2 +
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a5_s)[-1]))==1){
    ftc_a5[length(ftc_a5)+1] <- lm_esc_a5_s
    names(ftc_a5)[length(ftc_a5)] <- PRSs[i]
  }
  a5_models_esc[[i]] <- lm_esc_a5_s
  
  lm_nor_a5_s <- lmer( asec5wk ~ week + week2 +
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a5_s)[-1]))==1){
    ftc_a5[length(ftc_a5)+1] <- lm_nor_a5_s
    names(ftc_a5)[length(ftc_a5)] <- PRSs[i]
  }
  a5_models_nor[[i]] <- lm_nor_a5_s
}  

noPRS_both_a5 <- lmer( asec5wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a5)[-1]))==1){
  ftc_a5[length(ftc_a5)+1] <- noPRS_both_a5
  names(ftc_a5)[length(ftc_a5)] <- "No_PRS"
}
a5_models_both[[31]] <- noPRS_both_a5

noPRS_esc_a5 <- lmer( asec5wk ~ week + week2 +
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a5)[-1]))==1){
  ftc_a5[length(ftc_a5)+1] <- noPRS_esc_a5
  names(ftc_a5)[length(ftc_a5)] <- "No_PRS"
}
a5_models_esc[[31]] <- noPRS_esc_a5

noPRS_nor_a5 <- lmer( asec5wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a5)[-1]))==1){
  ftc_a5[length(ftc_a5)+1] <- noPRS_nor_a5
  names(ftc_a5)[length(ftc_a5)] <- "No_PRS"
}
a5_models_nor[[31]] <- noPRS_nor_a5

names(a5_models_both) <- PRSs_31
names(a5_models_esc) <- PRSs_31
names(a5_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a5 #####

a5_models_both_betas <- vector(mode = "list", length = length(PRSs))
a5_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a5_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a5_models_both_betas[i] <- list(summary(a5_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a5_models_esc_betas[i] <- list(summary(a5_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a5_models_nor_betas[i] <- list(summary(a5_models_nor[[i]])$coefficients[6,])
}

names(a5_models_both_betas) <- PRSs
names(a5_models_esc_betas) <- PRSs
names(a5_models_nor_betas) <- PRSs

a5_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a5_models_both_betas))[,1],
    t(as.data.frame(a5_models_both_betas))[,2],
    t(as.data.frame(a5_models_both_betas))[,5],
    t(as.data.frame(a5_models_esc_betas))[,1],
    t(as.data.frame(a5_models_esc_betas))[,2],
    t(as.data.frame(a5_models_esc_betas))[,5],
    t(as.data.frame(a5_models_nor_betas))[,1],
    t(as.data.frame(a5_models_nor_betas))[,2],
    t(as.data.frame(a5_models_nor_betas))[,5]
  )
)

colnames( a5_beta_matrix ) <- c(
  "a5_both_B", "a5_both_se", "a5_both_p", 
  "a5_esc_B", "a5_esc_se", "a5_esc_p",
  "a5_nor_B", "a5_nor_se", "a5_nor_p"
)

##### Get PRSs' explained variance for a5 #####

a5_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a5_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a5_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a5_matrix_both[i,] <- r.squaredGLMM(a5_models_both[[i]])
  a5_matrix_esc[i,] <- r.squaredGLMM(a5_models_esc[[i]])
  a5_matrix_nor[i,] <- r.squaredGLMM(a5_models_nor[[i]])
}

##### Make output files for a5 #####

sink( "a5.csv")

# Betas and p-values

cat("a5_betas\n,")
write.table( a5_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a5_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a5_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a5_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a5 heatmap #####

a5_models_betas <- c("a5_both_B", "a5_esc_B", "a5_nor_B")
a5_models_ses <- c("a5_both_se", "a5_esc_se", "a5_nor_se")
a5_models_ps <- c("a5_both_p", "a5_esc_p", "a5_nor_p")

a5_models_betas_columns <- which( colnames(a5_beta_matrix)
                                  %in% a5_models_betas )
a5_models_ses_columns <- which( colnames(a5_beta_matrix)
                                %in% a5_models_ses )
a5_models_ps_columns <- which( colnames(a5_beta_matrix)
                               %in% a5_models_ps )

a5_betas_only <- cbind(a5_beta_matrix[,a5_models_betas_columns[1]],
                       a5_beta_matrix[,a5_models_betas_columns[2]],
                       a5_beta_matrix[,a5_models_betas_columns[3]])

colnames(a5_betas_only) <- c("a5_both", "a5_esc","a5_nor")

a5_ses_only <- cbind(a5_beta_matrix[,a5_models_ses_columns[1]],
                     a5_beta_matrix[,a5_models_ses_columns[2]],
                     a5_beta_matrix[,a5_models_ses_columns[3]])

colnames(a5_betas_only) <- c("a5_both", "a5_esc","a5_nor")

a5_p_only <- cbind(a5_beta_matrix[,a5_models_ps_columns[1]],
                   a5_beta_matrix[,a5_models_ps_columns[2]],
                   a5_beta_matrix[,a5_models_ps_columns[3]])

colnames(a5_p_only) <- c("a5_both", "a5_esc","a5_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a5_both",
  "a5_esc",
  "a5_nor")))

lowest_ps[which(a5_p_only<=0.05)]<-a5_p_only[which(a5_p_only<=0.05)]

lowest_ps_a5 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a5_both",
                                                            "a5_esc",
                                                            "a5_nor")))

lowest_ps_a5[which(a5_p_only<=0.05)]<-a5_p_only[which(a5_p_only<=0.05)]

heatmap.2( a5_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a5,5),notecol="black",
           main="Headache"
)

##### Fit models for a6: constipation (cholinergic) #####

a6_models_both <- vector(mode = "list", length = length(PRSs_31))
a6_models_esc <- vector(mode = "list", length = length(PRSs_31))
a6_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a6 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a6_s <- lmer( asec6wk ~ week + 
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a6_s)[-1]))==1){
    ftc_a6[length(ftc_a6)+1] <- lm_both_a6_s
    names(ftc_a6)[length(ftc_a6)] <- PRSs[i]
  }
  a6_models_both[[i]] <- lm_both_a6_s
  
  lm_esc_a6_s <- lmer( asec6wk ~ week + 
                         sex + GENDEP_escit[,PRSs[i]] +
                         (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a6_s)[-1]))==1){
    ftc_a6[length(ftc_a6)+1] <- lm_esc_a6_s
    names(ftc_a6)[length(ftc_a6)] <- PRSs[i]
  }
  a6_models_esc[[i]] <- lm_esc_a6_s
  
  lm_nor_a6_s <- lmer( asec6wk ~ week + week2 +
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a6_s)[-1]))==1){
    ftc_a6[length(ftc_a6)+1] <- lm_nor_a6_s
    names(ftc_a6)[length(ftc_a6)] <- PRSs[i]
  }
  a6_models_nor[[i]] <- lm_nor_a6_s
}  

noPRS_both_a6 <- lmer( asec6wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a6)[-1]))==1){
  ftc_a6[length(ftc_a6)+1] <- noPRS_both_a6
  names(ftc_a6)[length(ftc_a6)] <- "No_PRS"
}
a6_models_both[[31]] <- noPRS_both_a6

noPRS_esc_a6 <- lmer( asec6wk ~ week + 
                        sex +
                        (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a6)[-1]))==1){
  ftc_a6[length(ftc_a6)+1] <- noPRS_esc_a6
  names(ftc_a6)[length(ftc_a6)] <- "No_PRS"
}
a6_models_esc[[31]] <- noPRS_esc_a6

noPRS_nor_a6 <- lmer( asec6wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a6)[-1]))==1){
  ftc_a6[length(ftc_a6)+1] <- noPRS_nor_a6
  names(ftc_a6)[length(ftc_a6)] <- "No_PRS"
}
a6_models_nor[[31]] <- noPRS_nor_a6

names(a6_models_both) <- PRSs_31
names(a6_models_esc) <- PRSs_31
names(a6_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a6 #####

a6_models_both_betas <- vector(mode = "list", length = length(PRSs))
a6_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a6_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a6_models_both_betas[i] <- list(summary(a6_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a6_models_esc_betas[i] <- list(summary(a6_models_esc[[i]])$coefficients[4,])
}
for (i in 1:30){
  a6_models_nor_betas[i] <- list(summary(a6_models_nor[[i]])$coefficients[6,])
}

names(a6_models_both_betas) <- PRSs
names(a6_models_esc_betas) <- PRSs
names(a6_models_nor_betas) <- PRSs

a6_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a6_models_both_betas))[,1],
    t(as.data.frame(a6_models_both_betas))[,2],
    t(as.data.frame(a6_models_both_betas))[,5],
    t(as.data.frame(a6_models_esc_betas))[,1],
    t(as.data.frame(a6_models_esc_betas))[,2],
    t(as.data.frame(a6_models_esc_betas))[,5],
    t(as.data.frame(a6_models_nor_betas))[,1],
    t(as.data.frame(a6_models_nor_betas))[,2],
    t(as.data.frame(a6_models_nor_betas))[,5]
  )
)

colnames( a6_beta_matrix ) <- c(
  "a6_both_B", "a6_both_se", "a6_both_p", 
  "a6_esc_B", "a6_esc_se", "a6_esc_p",
  "a6_nor_B", "a6_nor_se", "a6_nor_p"
)

##### Get PRSs' explained variance for a6 #####

a6_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a6_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a6_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a6_matrix_both[i,] <- r.squaredGLMM(a6_models_both[[i]])
  a6_matrix_esc[i,] <- r.squaredGLMM(a6_models_esc[[i]])
  a6_matrix_nor[i,] <- r.squaredGLMM(a6_models_nor[[i]])
}

##### Make output files for a6 #####

sink( "a6.csv")

# Betas and p-values

cat("a6_betas\n,")
write.table( a6_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a6_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a6_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a6_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a6 heatmap #####

a6_models_betas <- c("a6_both_B", "a6_esc_B", "a6_nor_B")
a6_models_ses <- c("a6_both_se", "a6_esc_se", "a6_nor_se")
a6_models_ps <- c("a6_both_p", "a6_esc_p", "a6_nor_p")

a6_models_betas_columns <- which( colnames(a6_beta_matrix)
                                  %in% a6_models_betas )
a6_models_ses_columns <- which( colnames(a6_beta_matrix)
                                %in% a6_models_ses )
a6_models_ps_columns <- which( colnames(a6_beta_matrix)
                               %in% a6_models_ps )

a6_betas_only <- cbind(a6_beta_matrix[,a6_models_betas_columns[1]],
                       a6_beta_matrix[,a6_models_betas_columns[2]],
                       a6_beta_matrix[,a6_models_betas_columns[3]])

colnames(a6_betas_only) <- c("a6_both", "a6_esc","a6_nor")

a6_ses_only <- cbind(a6_beta_matrix[,a6_models_ses_columns[1]],
                     a6_beta_matrix[,a6_models_ses_columns[2]],
                     a6_beta_matrix[,a6_models_ses_columns[3]])

colnames(a6_betas_only) <- c("a6_both", "a6_esc","a6_nor")

a6_p_only <- cbind(a6_beta_matrix[,a6_models_ps_columns[1]],
                   a6_beta_matrix[,a6_models_ps_columns[2]],
                   a6_beta_matrix[,a6_models_ps_columns[3]])

colnames(a6_p_only) <- c("a6_both", "a6_esc","a6_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a6_both",
  "a6_esc",
  "a6_nor")))

lowest_ps[which(a6_p_only<=0.05)]<-a6_p_only[which(a6_p_only<=0.05)]

lowest_ps_a6 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a6_both",
                                                            "a6_esc",
                                                            "a6_nor")))

lowest_ps_a6[which(a6_p_only<=0.05)]<-a6_p_only[which(a6_p_only<=0.05)]

heatmap.2( a6_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a6,5),notecol="black",
           main="Constipation"
)


##### Fit models for a7: diarrhoea (serotonergic) #####

a7_models_both <- vector(mode = "list", length = length(PRSs_31))
a7_models_esc <- vector(mode = "list", length = length(PRSs_31))
a7_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a7 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a7_s <- lmer( asec7wk ~ week + week2 +
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a7_s)[-1]))==1){
    ftc_a7[length(ftc_a7)+1] <- lm_both_a7_s
    names(ftc_a7)[length(ftc_a7)] <- PRSs[i]
  }
  a7_models_both[[i]] <- lm_both_a7_s
  
  lm_esc_a7_s <- lmer( asec7wk ~ week + 
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a7_s)[-1]))==1){
    ftc_a7[length(ftc_a7)+1] <- lm_esc_a7_s
    names(ftc_a7)[length(ftc_a7)] <- PRSs[i]
  }
  a7_models_esc[[i]] <- lm_esc_a7_s
  
  lm_nor_a7_s <- lmer( asec7wk ~ week + week2 +
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a7_s)[-1]))==1){
    ftc_a7[length(ftc_a7)+1] <- lm_nor_a7_s
    names(ftc_a7)[length(ftc_a7)] <- PRSs[i]
  }
  a7_models_nor[[i]] <- lm_nor_a7_s
}  

noPRS_both_a7 <- lmer( asec7wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a7)[-1]))==1){
  ftc_a7[length(ftc_a7)+1] <- noPRS_both_a7
  names(ftc_a7)[length(ftc_a7)] <- "No_PRS"
}
a7_models_both[[31]] <- noPRS_both_a7

noPRS_esc_a7 <- lmer( asec7wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a7)[-1]))==1){
  ftc_a7[length(ftc_a7)+1] <- noPRS_esc_a7
  names(ftc_a7)[length(ftc_a7)] <- "No_PRS"
}
a7_models_esc[[31]] <- noPRS_esc_a7

noPRS_nor_a7 <- lmer( asec7wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a7)[-1]))==1){
  ftc_a7[length(ftc_a7)+1] <- noPRS_nor_a7
  names(ftc_a7)[length(ftc_a7)] <- "No_PRS"
}
a7_models_nor[[31]] <- noPRS_nor_a7

names(a7_models_both) <- PRSs_31
names(a7_models_esc) <- PRSs_31
names(a7_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a7 #####

a7_models_both_betas <- vector(mode = "list", length = length(PRSs))
a7_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a7_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a7_models_both_betas[i] <- list(summary(a7_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a7_models_esc_betas[i] <- list(summary(a7_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a7_models_nor_betas[i] <- list(summary(a7_models_nor[[i]])$coefficients[6,])
}

names(a7_models_both_betas) <- PRSs
names(a7_models_esc_betas) <- PRSs
names(a7_models_nor_betas) <- PRSs

a7_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a7_models_both_betas))[,1],
    t(as.data.frame(a7_models_both_betas))[,2],
    t(as.data.frame(a7_models_both_betas))[,5],
    t(as.data.frame(a7_models_esc_betas))[,1],
    t(as.data.frame(a7_models_esc_betas))[,2],
    t(as.data.frame(a7_models_esc_betas))[,5],
    t(as.data.frame(a7_models_nor_betas))[,1],
    t(as.data.frame(a7_models_nor_betas))[,2],
    t(as.data.frame(a7_models_nor_betas))[,5]
  )
)

colnames( a7_beta_matrix ) <- c(
  "a7_both_B", "a7_both_se", "a7_both_p", 
  "a7_esc_B", "a7_esc_se", "a7_esc_p",
  "a7_nor_B", "a7_nor_se", "a7_nor_p"
)

##### Get PRSs' explained variance for a7 #####

a7_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a7_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a7_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a7_matrix_both[i,] <- r.squaredGLMM(a7_models_both[[i]])
  a7_matrix_esc[i,] <- r.squaredGLMM(a7_models_esc[[i]])
  a7_matrix_nor[i,] <- r.squaredGLMM(a7_models_nor[[i]])
}

##### Make output files for a7 #####

sink( "a7.csv")

# Betas and p-values

cat("a7_betas\n,")
write.table( a7_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a7_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a7_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a7_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a7 heatmap #####

a7_models_betas <- c("a7_both_B", "a7_esc_B", "a7_nor_B")
a7_models_ses <- c("a7_both_se", "a7_esc_se", "a7_nor_se")
a7_models_ps <- c("a7_both_p", "a7_esc_p", "a7_nor_p")

a7_models_betas_columns <- which( colnames(a7_beta_matrix)
                                  %in% a7_models_betas )
a7_models_ses_columns <- which( colnames(a7_beta_matrix)
                                %in% a7_models_ses )
a7_models_ps_columns <- which( colnames(a7_beta_matrix)
                               %in% a7_models_ps )

a7_betas_only <- cbind(a7_beta_matrix[,a7_models_betas_columns[1]],
                       a7_beta_matrix[,a7_models_betas_columns[2]],
                       a7_beta_matrix[,a7_models_betas_columns[3]])

colnames(a7_betas_only) <- c("a7_both", "a7_esc","a7_nor")

a7_ses_only <- cbind(a7_beta_matrix[,a7_models_ses_columns[1]],
                     a7_beta_matrix[,a7_models_ses_columns[2]],
                     a7_beta_matrix[,a7_models_ses_columns[3]])

colnames(a7_betas_only) <- c("a7_both", "a7_esc","a7_nor")

a7_p_only <- cbind(a7_beta_matrix[,a7_models_ps_columns[1]],
                   a7_beta_matrix[,a7_models_ps_columns[2]],
                   a7_beta_matrix[,a7_models_ps_columns[3]])

colnames(a7_p_only) <- c("a7_both", "a7_esc","a7_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a7_both",
  "a7_esc",
  "a7_nor")))

lowest_ps[which(a7_p_only<=0.05)]<-a7_p_only[which(a7_p_only<=0.05)]

lowest_ps_a7 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a7_both",
                                                            "a7_esc",
                                                            "a7_nor")))

lowest_ps_a7[which(a7_p_only<=0.05)]<-a7_p_only[which(a7_p_only<=0.05)]

heatmap.2( a7_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a7,5),notecol="black",
           main="Diarrhoea"
)

##### Fit models for a8: increased appetite (histaminergic) #####

a8_models_both <- vector(mode = "list", length = length(PRSs_31))
a8_models_esc <- vector(mode = "list", length = length(PRSs_31))
a8_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a8 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a8_s <- lmer( asec8wk ~ week + 
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a8_s)[-1]))==1){
    ftc_a8[length(ftc_a8)+1] <- lm_both_a8_s
    names(ftc_a8)[length(ftc_a8)] <- PRSs[i]
  }
  a8_models_both[[i]] <- lm_both_a8_s
  
  lm_esc_a8_s <- lmer( asec8wk ~ week + 
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a8_s)[-1]))==1){
    ftc_a8[length(ftc_a8)+1] <- lm_esc_a8_s
    names(ftc_a8)[length(ftc_a8)] <- PRSs[i]
  }
  a8_models_esc[[i]] <- lm_esc_a8_s
  
  lm_nor_a8_s <- lmer( asec8wk ~ week + 
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a8_s)[-1]))==1){
    ftc_a8[length(ftc_a8)+1] <- lm_nor_a8_s
    names(ftc_a8)[length(ftc_a8)] <- PRSs[i]
  }
  a8_models_nor[[i]] <- lm_nor_a8_s
}  

noPRS_both_a8 <- lmer( asec8wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a8)[-1]))==1){
  ftc_a8[length(ftc_a8)+1] <- noPRS_both_a8
  names(ftc_a8)[length(ftc_a8)] <- "No_PRS"
}
a8_models_both[[31]] <- noPRS_both_a8

noPRS_esc_a8 <- lmer( asec8wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a8)[-1]))==1){
  ftc_a8[length(ftc_a8)+1] <- noPRS_esc_a8
  names(ftc_a8)[length(ftc_a8)] <- "No_PRS"
}
a8_models_esc[[31]] <- noPRS_esc_a8

noPRS_nor_a8 <- lmer( asec8wk ~ week + 
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a8)[-1]))==1){
  ftc_a8[length(ftc_a8)+1] <- noPRS_nor_a8
  names(ftc_a8)[length(ftc_a8)] <- "No_PRS"
}
a8_models_nor[[31]] <- noPRS_nor_a8

names(a8_models_both) <- PRSs_31
names(a8_models_esc) <- PRSs_31
names(a8_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a8 #####

a8_models_both_betas <- vector(mode = "list", length = length(PRSs))
a8_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a8_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a8_models_both_betas[i] <- list(summary(a8_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a8_models_esc_betas[i] <- list(summary(a8_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a8_models_nor_betas[i] <- list(summary(a8_models_nor[[i]])$coefficients[5,])
}

names(a8_models_both_betas) <- PRSs
names(a8_models_esc_betas) <- PRSs
names(a8_models_nor_betas) <- PRSs

a8_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a8_models_both_betas))[,1],
    t(as.data.frame(a8_models_both_betas))[,2],
    t(as.data.frame(a8_models_both_betas))[,5],
    t(as.data.frame(a8_models_esc_betas))[,1],
    t(as.data.frame(a8_models_esc_betas))[,2],
    t(as.data.frame(a8_models_esc_betas))[,5],
    t(as.data.frame(a8_models_nor_betas))[,1],
    t(as.data.frame(a8_models_nor_betas))[,2],
    t(as.data.frame(a8_models_nor_betas))[,5]
  )
)

colnames( a8_beta_matrix ) <- c(
  "a8_both_B", "a8_both_se", "a8_both_p", 
  "a8_esc_B", "a8_esc_se", "a8_esc_p",
  "a8_nor_B", "a8_nor_se", "a8_nor_p"
)

##### Get PRSs' explained variance for a8 #####

a8_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a8_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a8_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a8_matrix_both[i,] <- r.squaredGLMM(a8_models_both[[i]])
  a8_matrix_esc[i,] <- r.squaredGLMM(a8_models_esc[[i]])
  a8_matrix_nor[i,] <- r.squaredGLMM(a8_models_nor[[i]])
}

##### Make output files for a8 #####

sink( "a8.csv")

# Betas and p-values

cat("a8_betas\n,")
write.table( a8_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a8_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a8_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a8_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a8 heatmap #####

a8_models_betas <- c("a8_both_B", "a8_esc_B", "a8_nor_B")
a8_models_ses <- c("a8_both_se", "a8_esc_se", "a8_nor_se")
a8_models_ps <- c("a8_both_p", "a8_esc_p", "a8_nor_p")

a8_models_betas_columns <- which( colnames(a8_beta_matrix)
                                  %in% a8_models_betas )
a8_models_ses_columns <- which( colnames(a8_beta_matrix)
                                %in% a8_models_ses )
a8_models_ps_columns <- which( colnames(a8_beta_matrix)
                               %in% a8_models_ps )

a8_betas_only <- cbind(a8_beta_matrix[,a8_models_betas_columns[1]],
                       a8_beta_matrix[,a8_models_betas_columns[2]],
                       a8_beta_matrix[,a8_models_betas_columns[3]])

colnames(a8_betas_only) <- c("a8_both", "a8_esc","a8_nor")

a8_ses_only <- cbind(a8_beta_matrix[,a8_models_ses_columns[1]],
                     a8_beta_matrix[,a8_models_ses_columns[2]],
                     a8_beta_matrix[,a8_models_ses_columns[3]])

colnames(a8_betas_only) <- c("a8_both", "a8_esc","a8_nor")

a8_p_only <- cbind(a8_beta_matrix[,a8_models_ps_columns[1]],
                   a8_beta_matrix[,a8_models_ps_columns[2]],
                   a8_beta_matrix[,a8_models_ps_columns[3]])

colnames(a8_p_only) <- c("a8_both", "a8_esc","a8_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a8_both",
  "a8_esc",
  "a8_nor")))

lowest_ps[which(a8_p_only<=0.05)]<-a8_p_only[which(a8_p_only<=0.05)]

lowest_ps_a8 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a8_both",
                                                            "a8_esc",
                                                            "a8_nor")))

lowest_ps_a8[which(a8_p_only<=0.05)]<-a8_p_only[which(a8_p_only<=0.05)]

heatmap.2( a8_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a8,5),notecol="black",
           main="Increased appetite"
)


##### Fit models for a9: decreased appetite (serotonergic) #####

a9_models_both <- vector(mode = "list", length = length(PRSs_31))
a9_models_esc <- vector(mode = "list", length = length(PRSs_31))
a9_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a9 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a9_s <- lmer( asec9wk ~ week + week2 +
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a9_s)[-1]))==1){
    ftc_a9[length(ftc_a9)+1] <- lm_both_a9_s
    names(ftc_a9)[length(ftc_a9)] <- PRSs[i]
  }
  a9_models_both[[i]] <- lm_both_a9_s

  lm_esc_a9_s <- lmer( asec9wk ~ week + week2 +
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a9_s)[-1]))==1){
    ftc_a9[length(ftc_a9)+1] <- lm_esc_a9_s
    names(ftc_a9)[length(ftc_a9)] <- PRSs[i]
  }
  a9_models_esc[[i]] <- lm_esc_a9_s

  lm_nor_a9_s <- lmer( asec9wk ~ week + week2 +
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a9_s)[-1]))==1){
    ftc_a9[length(ftc_a9)+1] <- lm_nor_a9_s
    names(ftc_a9)[length(ftc_a9)] <- PRSs[i]
  }
  a9_models_nor[[i]] <- lm_nor_a9_s
}  

noPRS_both_a9 <- lmer( asec9wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a9)[-1]))==1){
  ftc_a9[length(ftc_a9)+1] <- noPRS_both_a9
  names(ftc_a9)[length(ftc_a9)] <- "No_PRS"
}
a9_models_both[[31]] <- noPRS_both_a9

noPRS_esc_a9 <- lmer( asec9wk ~ week + week2 +
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a9)[-1]))==1){
  ftc_a9[length(ftc_a9)+1] <- noPRS_esc_a9
  names(ftc_a9)[length(ftc_a9)] <- "No_PRS"
}
a9_models_esc[[31]] <- noPRS_esc_a9

noPRS_nor_a9 <- lmer( asec9wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a9)[-1]))==1){
  ftc_a9[length(ftc_a9)+1] <- noPRS_nor_a9
  names(ftc_a9)[length(ftc_a9)] <- "No_PRS"
}
a9_models_nor[[31]] <- noPRS_nor_a9

names(a9_models_both) <- PRSs_31
names(a9_models_esc) <- PRSs_31
names(a9_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a9 #####

a9_models_both_betas <- vector(mode = "list", length = length(PRSs))
a9_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a9_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a9_models_both_betas[i] <- list(summary(a9_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a9_models_esc_betas[i] <- list(summary(a9_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a9_models_nor_betas[i] <- list(summary(a9_models_nor[[i]])$coefficients[6,])
}

names(a9_models_both_betas) <- PRSs
names(a9_models_esc_betas) <- PRSs
names(a9_models_nor_betas) <- PRSs

a9_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a9_models_both_betas))[,1],
    t(as.data.frame(a9_models_both_betas))[,2],
    t(as.data.frame(a9_models_both_betas))[,5],
    t(as.data.frame(a9_models_esc_betas))[,1],
    t(as.data.frame(a9_models_esc_betas))[,2],
    t(as.data.frame(a9_models_esc_betas))[,5],
    t(as.data.frame(a9_models_nor_betas))[,1],
    t(as.data.frame(a9_models_nor_betas))[,2],
    t(as.data.frame(a9_models_nor_betas))[,5]
  )
)

colnames( a9_beta_matrix ) <- c(
  "a9_both_B", "a9_both_se", "a9_both_p", 
  "a9_esc_B", "a9_esc_se", "a9_esc_p",
  "a9_nor_B", "a9_nor_se", "a9_nor_p"
)

##### Get PRSs' explained variance for a9 #####

a9_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a9_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a9_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a9_matrix_both[i,] <- r.squaredGLMM(a9_models_both[[i]])
  a9_matrix_esc[i,] <- r.squaredGLMM(a9_models_esc[[i]])
  a9_matrix_nor[i,] <- r.squaredGLMM(a9_models_nor[[i]])
}

##### Make output files for a9 #####

sink( "a9.csv")

# Betas and p-values

cat("a9_betas\n,")
write.table( a9_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a9_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a9_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a9_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a9 heatmap #####

a9_models_betas <- c("a9_both_B", "a9_esc_B", "a9_nor_B")
a9_models_ses <- c("a9_both_se", "a9_esc_se", "a9_nor_se")
a9_models_ps <- c("a9_both_p", "a9_esc_p", "a9_nor_p")

a9_models_betas_columns <- which( colnames(a9_beta_matrix)
                                  %in% a9_models_betas )
a9_models_ses_columns <- which( colnames(a9_beta_matrix)
                                %in% a9_models_ses )
a9_models_ps_columns <- which( colnames(a9_beta_matrix)
                               %in% a9_models_ps )

a9_betas_only <- cbind(a9_beta_matrix[,a9_models_betas_columns[1]],
                       a9_beta_matrix[,a9_models_betas_columns[2]],
                       a9_beta_matrix[,a9_models_betas_columns[3]])

colnames(a9_betas_only) <- c("a9_both", "a9_esc","a9_nor")

a9_ses_only <- cbind(a9_beta_matrix[,a9_models_ses_columns[1]],
                     a9_beta_matrix[,a9_models_ses_columns[2]],
                     a9_beta_matrix[,a9_models_ses_columns[3]])

colnames(a9_betas_only) <- c("a9_both", "a9_esc","a9_nor")

a9_p_only <- cbind(a9_beta_matrix[,a9_models_ps_columns[1]],
                   a9_beta_matrix[,a9_models_ps_columns[2]],
                   a9_beta_matrix[,a9_models_ps_columns[3]])

colnames(a9_p_only) <- c("a9_both", "a9_esc","a9_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a9_both",
  "a9_esc",
  "a9_nor")))

lowest_ps[which(a9_p_only<=0.05)]<-a9_p_only[which(a9_p_only<=0.05)]

lowest_ps_a9 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a9_both",
                                                            "a9_esc",
                                                            "a9_nor")))

lowest_ps_a9[which(a9_p_only<=0.05)]<-a9_p_only[which(a9_p_only<=0.05)]

heatmap.2( a9_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a9,5),notecol="black",
           main="Decreased appetite"
)

##### Fit models for a10: nausea or vomiting (serotonergic) #####

a10_models_both <- vector(mode = "list", length = length(PRSs_31))
a10_models_esc <- vector(mode = "list", length = length(PRSs_31))
a10_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a10 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a10_s <- lmer( asec10wk ~ week + week2 +
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a10_s)[-1]))==1){
    ftc_a10[length(ftc_a10)+1] <- lm_both_a10_s
    names(ftc_a10)[length(ftc_a10)] <- PRSs[i]
  }
  a10_models_both[[i]] <- lm_both_a10_s
  
  lm_esc_a10_s <- lmer( asec10wk ~ week + week2 +
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a10_s)[-1]))==1){
    ftc_a10[length(ftc_a10)+1] <- lm_esc_a10_s
    names(ftc_a10)[length(ftc_a10)] <- PRSs[i]
  }
  a10_models_esc[[i]] <- lm_esc_a10_s
  
  lm_nor_a10_s <- lmer( asec10wk ~ week + 
                         sex + cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a10_s)[-1]))==1){
    ftc_a10[length(ftc_a10)+1] <- lm_nor_a10_s
    names(ftc_a10)[length(ftc_a10)] <- PRSs[i]
  }
  a10_models_nor[[i]] <- lm_nor_a10_s
}  

noPRS_both_a10 <- lmer( asec10wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a10)[-1]))==1){
  ftc_a10[length(ftc_a10)+1] <- noPRS_both_a10
  names(ftc_a10)[length(ftc_a10)] <- "No_PRS"
}
a10_models_both[[31]] <- noPRS_both_a10

noPRS_esc_a10 <- lmer( asec10wk ~ week + week2 +
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a10)[-1]))==1){
  ftc_a10[length(ftc_a10)+1] <- noPRS_esc_a10
  names(ftc_a10)[length(ftc_a10)] <- "No_PRS"
}
a10_models_esc[[31]] <- noPRS_esc_a10

noPRS_nor_a10 <- lmer( asec10wk ~ week + 
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a10)[-1]))==1){
  ftc_a10[length(ftc_a10)+1] <- noPRS_nor_a10
  names(ftc_a10)[length(ftc_a10)] <- "No_PRS"
}
a10_models_nor[[31]] <- noPRS_nor_a10

names(a10_models_both) <- PRSs_31
names(a10_models_esc) <- PRSs_31
names(a10_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a10 #####

a10_models_both_betas <- vector(mode = "list", length = length(PRSs))
a10_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a10_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a10_models_both_betas[i] <- list(summary(a10_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a10_models_esc_betas[i] <- list(summary(a10_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a10_models_nor_betas[i] <- list(summary(a10_models_nor[[i]])$coefficients[5,])
}

names(a10_models_both_betas) <- PRSs
names(a10_models_esc_betas) <- PRSs
names(a10_models_nor_betas) <- PRSs

a10_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a10_models_both_betas))[,1],
    t(as.data.frame(a10_models_both_betas))[,2],
    t(as.data.frame(a10_models_both_betas))[,5],
    t(as.data.frame(a10_models_esc_betas))[,1],
    t(as.data.frame(a10_models_esc_betas))[,2],
    t(as.data.frame(a10_models_esc_betas))[,5],
    t(as.data.frame(a10_models_nor_betas))[,1],
    t(as.data.frame(a10_models_nor_betas))[,2],
    t(as.data.frame(a10_models_nor_betas))[,5]
  )
)

colnames( a10_beta_matrix ) <- c(
  "a10_both_B", "a10_both_se", "a10_both_p", 
  "a10_esc_B", "a10_esc_se", "a10_esc_p",
  "a10_nor_B", "a10_nor_se", "a10_nor_p"
)

##### Get PRSs' explained variance for a10 #####

a10_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a10_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a10_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a10_matrix_both[i,] <- r.squaredGLMM(a10_models_both[[i]])
  a10_matrix_esc[i,] <- r.squaredGLMM(a10_models_esc[[i]])
  a10_matrix_nor[i,] <- r.squaredGLMM(a10_models_nor[[i]])
}

##### Make output files for a10 #####

sink( "a10.csv")

# Betas and p-values

cat("a10_betas\n,")
write.table( a10_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a10_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a10_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a10_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a10 heatmap #####

a10_models_betas <- c("a10_both_B", "a10_esc_B", "a10_nor_B")
a10_models_ses <- c("a10_both_se", "a10_esc_se", "a10_nor_se")
a10_models_ps <- c("a10_both_p", "a10_esc_p", "a10_nor_p")

a10_models_betas_columns <- which( colnames(a10_beta_matrix)
                                  %in% a10_models_betas )
a10_models_ses_columns <- which( colnames(a10_beta_matrix)
                                %in% a10_models_ses )
a10_models_ps_columns <- which( colnames(a10_beta_matrix)
                               %in% a10_models_ps )

a10_betas_only <- cbind(a10_beta_matrix[,a10_models_betas_columns[1]],
                       a10_beta_matrix[,a10_models_betas_columns[2]],
                       a10_beta_matrix[,a10_models_betas_columns[3]])

colnames(a10_betas_only) <- c("a10_both", "a10_esc","a10_nor")

a10_ses_only <- cbind(a10_beta_matrix[,a10_models_ses_columns[1]],
                     a10_beta_matrix[,a10_models_ses_columns[2]],
                     a10_beta_matrix[,a10_models_ses_columns[3]])

colnames(a10_betas_only) <- c("a10_both", "a10_esc","a10_nor")

a10_p_only <- cbind(a10_beta_matrix[,a10_models_ps_columns[1]],
                   a10_beta_matrix[,a10_models_ps_columns[2]],
                   a10_beta_matrix[,a10_models_ps_columns[3]])

colnames(a10_p_only) <- c("a10_both", "a10_esc","a10_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a10_both",
  "a10_esc",
  "a10_nor")))

lowest_ps[which(a10_p_only<=0.05)]<-a10_p_only[which(a10_p_only<=0.05)]

lowest_ps_a10 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a10_both",
                                                            "a10_esc",
                                                            "a10_nor")))

lowest_ps_a10[which(a10_p_only<=0.05)]<-a10_p_only[which(a10_p_only<=0.05)]

heatmap.2( a10_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a10,5),notecol="black",
           main="Nausea or vomiting"
)

##### Fit models for a11: problems with urination (cholinergic) #####

a11_models_both <- vector(mode = "list", length = length(PRSs_31))
a11_models_esc <- vector(mode = "list", length = length(PRSs_31))
a11_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a11 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a11_s <- lmer( asec11wk ~ week + 
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a11_s)[-1]))==1){
    ftc_a11[length(ftc_a11)+1] <- lm_both_a11_s
    names(ftc_a11)[length(ftc_a11)] <- PRSs[i]
  }
  a11_models_both[[i]] <- lm_both_a11_s

  lm_esc_a11_s <- lmer( asec11wk ~ week + 
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a11_s)[-1]))==1){
    ftc_a11[length(ftc_a11)+1] <- lm_esc_a11_s
    names(ftc_a11)[length(ftc_a11)] <- PRSs[i]
  }
  a11_models_esc[[i]] <- lm_esc_a11_s
  
  lm_nor_a11_s <- lmer( asec11wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a11_s)[-1]))==1){
    ftc_a11[length(ftc_a11)+1] <- lm_nor_a11_s
    names(ftc_a11)[length(ftc_a11)] <- PRSs[i]
  }
  a11_models_nor[[i]] <- lm_nor_a11_s
}  

noPRS_both_a11 <- lmer( asec11wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a11)[-1]))==1){
  ftc_a11[length(ftc_a11)+1] <- noPRS_both_a11
  names(ftc_a11)[length(ftc_a11)] <- "No_PRS"
}
a11_models_both[[31]] <- noPRS_both_a11

noPRS_esc_a11 <- lmer( asec11wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a11)[-1]))==1){
  ftc_a11[length(ftc_a11)+1] <- noPRS_esc_a11
  names(ftc_a11)[length(ftc_a11)] <- "No_PRS"
}
a11_models_esc[[31]] <- noPRS_esc_a11

noPRS_nor_a11 <- lmer( asec11wk ~ week + 
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a11)[-1]))==1){
  ftc_a11[length(ftc_a11)+1] <- noPRS_nor_a11
  names(ftc_a11)[length(ftc_a11)] <- "No_PRS"
}
a11_models_nor[[31]] <- noPRS_nor_a11

names(a11_models_both) <- PRSs_31
names(a11_models_esc) <- PRSs_31
names(a11_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a11 #####

a11_models_both_betas <- vector(mode = "list", length = length(PRSs))
a11_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a11_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a11_models_both_betas[i] <- list(summary(a11_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a11_models_esc_betas[i] <- list(summary(a11_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a11_models_nor_betas[i] <- list(summary(a11_models_nor[[i]])$coefficients[5,])
}

names(a11_models_both_betas) <- PRSs
names(a11_models_esc_betas) <- PRSs
names(a11_models_nor_betas) <- PRSs

a11_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a11_models_both_betas))[,1],
    t(as.data.frame(a11_models_both_betas))[,2],
    t(as.data.frame(a11_models_both_betas))[,5],
    t(as.data.frame(a11_models_esc_betas))[,1],
    t(as.data.frame(a11_models_esc_betas))[,2],
    t(as.data.frame(a11_models_esc_betas))[,5],
    t(as.data.frame(a11_models_nor_betas))[,1],
    t(as.data.frame(a11_models_nor_betas))[,2],
    t(as.data.frame(a11_models_nor_betas))[,5]
  )
)

colnames( a11_beta_matrix ) <- c(
  "a11_both_B", "a11_both_se", "a11_both_p", 
  "a11_esc_B", "a11_esc_se", "a11_esc_p",
  "a11_nor_B", "a11_nor_se", "a11_nor_p"
)

##### Get PRSs' explained variance for a11 #####

a11_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a11_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a11_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a11_matrix_both[i,] <- r.squaredGLMM(a11_models_both[[i]])
  a11_matrix_esc[i,] <- r.squaredGLMM(a11_models_esc[[i]])
  a11_matrix_nor[i,] <- r.squaredGLMM(a11_models_nor[[i]])
}

##### Make output files for a11 #####

sink( "a11.csv")

# Betas and p-values

cat("a11_betas\n,")
write.table( a11_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a11_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a11_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a11_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a11 heatmap #####

a11_models_betas <- c("a11_both_B", "a11_esc_B", "a11_nor_B")
a11_models_ses <- c("a11_both_se", "a11_esc_se", "a11_nor_se")
a11_models_ps <- c("a11_both_p", "a11_esc_p", "a11_nor_p")

a11_models_betas_columns <- which( colnames(a11_beta_matrix)
                                  %in% a11_models_betas )
a11_models_ses_columns <- which( colnames(a11_beta_matrix)
                                %in% a11_models_ses )
a11_models_ps_columns <- which( colnames(a11_beta_matrix)
                               %in% a11_models_ps )

a11_betas_only <- cbind(a11_beta_matrix[,a11_models_betas_columns[1]],
                       a11_beta_matrix[,a11_models_betas_columns[2]],
                       a11_beta_matrix[,a11_models_betas_columns[3]])

colnames(a11_betas_only) <- c("a11_both", "a11_esc","a11_nor")

a11_ses_only <- cbind(a11_beta_matrix[,a11_models_ses_columns[1]],
                     a11_beta_matrix[,a11_models_ses_columns[2]],
                     a11_beta_matrix[,a11_models_ses_columns[3]])

colnames(a11_betas_only) <- c("a11_both", "a11_esc","a11_nor")

a11_p_only <- cbind(a11_beta_matrix[,a11_models_ps_columns[1]],
                   a11_beta_matrix[,a11_models_ps_columns[2]],
                   a11_beta_matrix[,a11_models_ps_columns[3]])

colnames(a11_p_only) <- c("a11_both", "a11_esc","a11_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a11_both",
  "a11_esc",
  "a11_nor")))

lowest_ps[which(a11_p_only<=0.05)]<-a11_p_only[which(a11_p_only<=0.05)]

lowest_ps_a11 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a11_both",
                                                            "a11_esc",
                                                            "a11_nor")))

lowest_ps_a11[which(a11_p_only<=0.05)]<-a11_p_only[which(a11_p_only<=0.05)]

heatmap.2( a11_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a11,5),notecol="black",
           main="Urination problems"
)



##### Fit models for a12: problems with sexual function (serotonergic) #####

a12_models_both <- vector(mode = "list", length = length(PRSs_31))
a12_models_esc <- vector(mode = "list", length = length(PRSs_31))
a12_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a12 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a12_s <- lmer( asec12wk ~ week + 
                          sex + cage + drug + GENDEP[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a12_s)[-1]))==1){
    ftc_a12[length(ftc_a12)+1] <- lm_both_a12_s
    names(ftc_a12)[length(ftc_a12)] <- PRSs[i]
  }
  a12_models_both[[i]] <- lm_both_a12_s
  
  lm_esc_a12_s <- lmer( asec12wk ~ week + 
                         sex + cage + GENDEP_escit[,PRSs[i]] +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a12_s)[-1]))==1){
    ftc_a12[length(ftc_a12)+1] <- lm_esc_a12_s
    names(ftc_a12)[length(ftc_a12)] <- PRSs[i]
  }
  a12_models_esc[[i]] <- lm_esc_a12_s
  
  lm_nor_a12_s <- lmer( asec12wk ~ week + 
                         cage + GENDEP_nortrip[,PRSs[i]] +
                         zPC1 +
                          (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a12_s)[-1]))==1){
    ftc_a12[length(ftc_a12)+1] <- lm_nor_a12_s
    names(ftc_a12)[length(ftc_a12)] <- PRSs[i]
  }
  a12_models_nor[[i]] <- lm_nor_a12_s
}  

noPRS_both_a12 <- lmer( asec12wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a12)[-1]))==1){
  ftc_a12[length(ftc_a12)+1] <- noPRS_both_a12
  names(ftc_a12)[length(ftc_a12)] <- "No_PRS"
}
a12_models_both[[31]] <- noPRS_both_a12

noPRS_esc_a12 <- lmer( asec12wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a12)[-1]))==1){
  ftc_a12[length(ftc_a12)+1] <- noPRS_esc_a12
  names(ftc_a12)[length(ftc_a12)] <- "No_PRS"
}
a12_models_esc[[31]] <- noPRS_esc_a12

noPRS_nor_a12 <- lmer( asec12wk ~ week + 
                        cage +
                        zPC1 + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a12)[-1]))==1){
  ftc_a12[length(ftc_a12)+1] <- noPRS_nor_a12
  names(ftc_a12)[length(ftc_a12)] <- "No_PRS"
}
a12_models_nor[[31]] <- noPRS_nor_a12

names(a12_models_both) <- PRSs_31
names(a12_models_esc) <- PRSs_31
names(a12_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a12 #####

a12_models_both_betas <- vector(mode = "list", length = length(PRSs))
a12_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a12_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a12_models_both_betas[i] <- list(summary(a12_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a12_models_esc_betas[i] <- list(summary(a12_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a12_models_nor_betas[i] <- list(summary(a12_models_nor[[i]])$coefficients[4,])
}

names(a12_models_both_betas) <- PRSs
names(a12_models_esc_betas) <- PRSs
names(a12_models_nor_betas) <- PRSs

a12_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a12_models_both_betas))[,1],
    t(as.data.frame(a12_models_both_betas))[,2],
    t(as.data.frame(a12_models_both_betas))[,5],
    t(as.data.frame(a12_models_esc_betas))[,1],
    t(as.data.frame(a12_models_esc_betas))[,2],
    t(as.data.frame(a12_models_esc_betas))[,5],
    t(as.data.frame(a12_models_nor_betas))[,1],
    t(as.data.frame(a12_models_nor_betas))[,2],
    t(as.data.frame(a12_models_nor_betas))[,5]
  )
)

colnames( a12_beta_matrix ) <- c(
  "a12_both_B", "a12_both_se", "a12_both_p", 
  "a12_esc_B", "a12_esc_se", "a12_esc_p",
  "a12_nor_B", "a12_nor_se", "a12_nor_p"
)

##### Get PRSs' explained variance for a12 #####

a12_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a12_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a12_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a12_matrix_both[i,] <- r.squaredGLMM(a12_models_both[[i]])
  a12_matrix_esc[i,] <- r.squaredGLMM(a12_models_esc[[i]])
  a12_matrix_nor[i,] <- r.squaredGLMM(a12_models_nor[[i]])
}

##### Make output files for a12 #####

sink( "a12.csv")

# Betas and p-values

cat("a12_betas\n,")
write.table( a12_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a12_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a12_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a12_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a12 heatmap #####

a12_models_betas <- c("a12_both_B", "a12_esc_B", "a12_nor_B")
a12_models_ses <- c("a12_both_se", "a12_esc_se", "a12_nor_se")
a12_models_ps <- c("a12_both_p", "a12_esc_p", "a12_nor_p")

a12_models_betas_columns <- which( colnames(a12_beta_matrix)
                                  %in% a12_models_betas )
a12_models_ses_columns <- which( colnames(a12_beta_matrix)
                                %in% a12_models_ses )
a12_models_ps_columns <- which( colnames(a12_beta_matrix)
                               %in% a12_models_ps )

a12_betas_only <- cbind(a12_beta_matrix[,a12_models_betas_columns[1]],
                       a12_beta_matrix[,a12_models_betas_columns[2]],
                       a12_beta_matrix[,a12_models_betas_columns[3]])

colnames(a12_betas_only) <- c("a12_both", "a12_esc","a12_nor")

a12_ses_only <- cbind(a12_beta_matrix[,a12_models_ses_columns[1]],
                     a12_beta_matrix[,a12_models_ses_columns[2]],
                     a12_beta_matrix[,a12_models_ses_columns[3]])

colnames(a12_betas_only) <- c("a12_both", "a12_esc","a12_nor")

a12_p_only <- cbind(a12_beta_matrix[,a12_models_ps_columns[1]],
                   a12_beta_matrix[,a12_models_ps_columns[2]],
                   a12_beta_matrix[,a12_models_ps_columns[3]])

colnames(a12_p_only) <- c("a12_both", "a12_esc","a12_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a12_both",
  "a12_esc",
  "a12_nor")))

lowest_ps[which(a12_p_only<=0.05)]<-a12_p_only[which(a12_p_only<=0.05)]

lowest_ps_a12 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a12_both",
                                                            "a12_esc",
                                                            "a12_nor")))

lowest_ps_a12[which(a12_p_only<=0.05)]<-a12_p_only[which(a12_p_only<=0.05)]

heatmap.2( a12_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a12,5),notecol="black",
           main="Sexual dysfunction"
)

##### Fit models for a13: paliptations #####

a13_models_both <- vector(mode = "list", length = length(PRSs_31))
a13_models_esc <- vector(mode = "list", length = length(PRSs_31))
a13_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a13 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a13_s <- lmer( asec13wk ~ week + 
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a13_s)[-1]))==1){
    ftc_a13[length(ftc_a13)+1] <- lm_both_a13_s
    names(ftc_a13)[length(ftc_a13)] <- PRSs[i]
  }
  a13_models_both[[i]] <- lm_both_a13_s
  
  lm_esc_a13_s <- lmer( asec13wk ~ week + week2 +
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a13_s)[-1]))==1){
    ftc_a13[length(ftc_a13)+1] <- lm_esc_a13_s
    names(ftc_a13)[length(ftc_a13)] <- PRSs[i]
  }
  a13_models_esc[[i]] <- lm_esc_a13_s
  
  lm_nor_a13_s <- lmer( asec13wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a13_s)[-1]))==1){
    ftc_a13[length(ftc_a13)+1] <- lm_nor_a13_s
    names(ftc_a13)[length(ftc_a13)] <- PRSs[i]
  }
  a13_models_nor[[i]] <- lm_nor_a13_s
}  

noPRS_both_a13 <- lmer( asec13wk ~ week + 
                          sex + cage + drug + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a13)[-1]))==1){
  ftc_a13[length(ftc_a13)+1] <- noPRS_both_a13
  names(ftc_a13)[length(ftc_a13)] <- "No_PRS"
}
a13_models_both[[31]] <- noPRS_both_a13

noPRS_esc_a13 <- lmer( asec13wk ~ week + week2 +
                         sex + cage + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a13)[-1]))==1){
  ftc_a13[length(ftc_a13)+1] <- noPRS_esc_a13
  names(ftc_a13)[length(ftc_a13)] <- "No_PRS"
}
a13_models_esc[[31]] <- noPRS_esc_a13

noPRS_nor_a13 <- lmer( asec13wk ~ week + 
                         sex + cage +
                         zPC1 + zPC2 + zPC3 +
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a13)[-1]))==1){
  ftc_a13[length(ftc_a13)+1] <- noPRS_nor_a13
  names(ftc_a13)[length(ftc_a13)] <- "No_PRS"
}
a13_models_nor[[31]] <- noPRS_nor_a13

names(a13_models_both) <- PRSs_31
names(a13_models_esc) <- PRSs_31
names(a13_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a13 #####

a13_models_both_betas <- vector(mode = "list", length = length(PRSs))
a13_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a13_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a13_models_both_betas[i] <- list(summary(a13_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a13_models_esc_betas[i] <- list(summary(a13_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a13_models_nor_betas[i] <- list(summary(a13_models_nor[[i]])$coefficients[5,])
}

names(a13_models_both_betas) <- PRSs
names(a13_models_esc_betas) <- PRSs
names(a13_models_nor_betas) <- PRSs

a13_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a13_models_both_betas))[,1],
    t(as.data.frame(a13_models_both_betas))[,2],
    t(as.data.frame(a13_models_both_betas))[,5],
    t(as.data.frame(a13_models_esc_betas))[,1],
    t(as.data.frame(a13_models_esc_betas))[,2],
    t(as.data.frame(a13_models_esc_betas))[,5],
    t(as.data.frame(a13_models_nor_betas))[,1],
    t(as.data.frame(a13_models_nor_betas))[,2],
    t(as.data.frame(a13_models_nor_betas))[,5]
  )
)

colnames( a13_beta_matrix ) <- c(
  "a13_both_B", "a13_both_se", "a13_both_p", 
  "a13_esc_B", "a13_esc_se", "a13_esc_p",
  "a13_nor_B", "a13_nor_se", "a13_nor_p"
)

##### Get PRSs' explained variance for a13 #####

a13_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a13_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a13_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a13_matrix_both[i,] <- r.squaredGLMM(a13_models_both[[i]])
  a13_matrix_esc[i,] <- r.squaredGLMM(a13_models_esc[[i]])
  a13_matrix_nor[i,] <- r.squaredGLMM(a13_models_nor[[i]])
}

##### Make output files for a13 #####

sink( "a13.csv")

# Betas and p-values

cat("a13_betas\n,")
write.table( a13_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a13_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a13_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a13_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a13 heatmap #####

a13_models_betas <- c("a13_both_B", "a13_esc_B", "a13_nor_B")
a13_models_ses <- c("a13_both_se", "a13_esc_se", "a13_nor_se")
a13_models_ps <- c("a13_both_p", "a13_esc_p", "a13_nor_p")

a13_models_betas_columns <- which( colnames(a13_beta_matrix)
                                   %in% a13_models_betas )
a13_models_ses_columns <- which( colnames(a13_beta_matrix)
                                 %in% a13_models_ses )
a13_models_ps_columns <- which( colnames(a13_beta_matrix)
                                %in% a13_models_ps )

a13_betas_only <- cbind(a13_beta_matrix[,a13_models_betas_columns[1]],
                        a13_beta_matrix[,a13_models_betas_columns[2]],
                        a13_beta_matrix[,a13_models_betas_columns[3]])

colnames(a13_betas_only) <- c("a13_both", "a13_esc","a13_nor")

a13_ses_only <- cbind(a13_beta_matrix[,a13_models_ses_columns[1]],
                      a13_beta_matrix[,a13_models_ses_columns[2]],
                      a13_beta_matrix[,a13_models_ses_columns[3]])

colnames(a13_betas_only) <- c("a13_both", "a13_esc","a13_nor")

a13_p_only <- cbind(a13_beta_matrix[,a13_models_ps_columns[1]],
                    a13_beta_matrix[,a13_models_ps_columns[2]],
                    a13_beta_matrix[,a13_models_ps_columns[3]])

colnames(a13_p_only) <- c("a13_both", "a13_esc","a13_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a13_both",
  "a13_esc",
  "a13_nor")))

lowest_ps[which(a13_p_only<=0.05)]<-a13_p_only[which(a13_p_only<=0.05)]

lowest_ps_a13 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a13_both",
                                                             "a13_esc",
                                                             "a13_nor")))

lowest_ps_a13[which(a13_p_only<=0.05)]<-a13_p_only[which(a13_p_only<=0.05)]

heatmap.2( a13_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a13,5),notecol="black",
           main="Palpitations"
)

##### Fit models for a14: light-headedness upon standing (adrenergic) #####

a14_models_both <- vector(mode = "list", length = length(PRSs_31))
a14_models_esc <- vector(mode = "list", length = length(PRSs_31))
a14_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a14 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a14_s <- lmer( asec14wk ~ week + week2 +
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a14_s)[-1]))==1){
    ftc_a14[length(ftc_a14)+1] <- lm_both_a14_s
    names(ftc_a14)[length(ftc_a14)] <- PRSs[i]
  }
  a14_models_both[[i]] <- lm_both_a14_s

  lm_esc_a14_s <- lmer( asec14wk ~ week + week2 +
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a14_s)[-1]))==1){
    ftc_a14[length(ftc_a14)+1] <- lm_esc_a14_s
    names(ftc_a14)[length(ftc_a14)] <- PRSs[i]
  }
  a14_models_esc[[i]] <- lm_esc_a14_s

  lm_nor_a14_s <- lmer( asec14wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a14_s)[-1]))==1){
    ftc_a14[length(ftc_a14)+1] <- lm_nor_a14_s
    names(ftc_a14)[length(ftc_a14)] <- PRSs[i]
  }
  a14_models_nor[[i]] <- lm_nor_a14_s
}  

noPRS_both_a14 <- lmer( asec14wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a14)[-1]))==1){
  ftc_a14[length(ftc_a14)+1] <- noPRS_both_a14
  names(ftc_a14)[length(ftc_a14)] <- "No_PRS"
}
a14_models_both[[31]] <- noPRS_both_a14

noPRS_esc_a14 <- lmer( asec14wk ~ week + week2 +
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a14)[-1]))==1){
  ftc_a14[length(ftc_a14)+1] <- noPRS_esc_a14
  names(ftc_a14)[length(ftc_a14)] <- "No_PRS"
}
a14_models_esc[[31]] <- noPRS_esc_a14

noPRS_nor_a14 <- lmer( asec14wk ~ week + 
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a14)[-1]))==1){
  ftc_a14[length(ftc_a14)+1] <- noPRS_nor_a14
  names(ftc_a14)[length(ftc_a14)] <- "No_PRS"
}
a14_models_nor[[31]] <- noPRS_nor_a14

names(a14_models_both) <- PRSs_31
names(a14_models_esc) <- PRSs_31
names(a14_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a14 #####

a14_models_both_betas <- vector(mode = "list", length = length(PRSs))
a14_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a14_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a14_models_both_betas[i] <- list(summary(a14_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a14_models_esc_betas[i] <- list(summary(a14_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a14_models_nor_betas[i] <- list(summary(a14_models_nor[[i]])$coefficients[5,])
}

names(a14_models_both_betas) <- PRSs
names(a14_models_esc_betas) <- PRSs
names(a14_models_nor_betas) <- PRSs

a14_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a14_models_both_betas))[,1],
    t(as.data.frame(a14_models_both_betas))[,2],
    t(as.data.frame(a14_models_both_betas))[,5],
    t(as.data.frame(a14_models_esc_betas))[,1],
    t(as.data.frame(a14_models_esc_betas))[,2],
    t(as.data.frame(a14_models_esc_betas))[,5],
    t(as.data.frame(a14_models_nor_betas))[,1],
    t(as.data.frame(a14_models_nor_betas))[,2],
    t(as.data.frame(a14_models_nor_betas))[,5]
  )
)

colnames( a14_beta_matrix ) <- c(
  "a14_both_B", "a14_both_se", "a14_both_p", 
  "a14_esc_B", "a14_esc_se", "a14_esc_p",
  "a14_nor_B", "a14_nor_se", "a14_nor_p"
)

##### Get PRSs' explained variance for a14 #####

a14_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a14_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a14_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a14_matrix_both[i,] <- r.squaredGLMM(a14_models_both[[i]])
  a14_matrix_esc[i,] <- r.squaredGLMM(a14_models_esc[[i]])
  a14_matrix_nor[i,] <- r.squaredGLMM(a14_models_nor[[i]])
}

##### Make output files for a14 #####

sink( "a14.csv")

# Betas and p-values

cat("a14_betas\n,")
write.table( a14_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a14_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a14_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a14_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a14 heatmap #####

a14_models_betas <- c("a14_both_B", "a14_esc_B", "a14_nor_B")
a14_models_ses <- c("a14_both_se", "a14_esc_se", "a14_nor_se")
a14_models_ps <- c("a14_both_p", "a14_esc_p", "a14_nor_p")

a14_models_betas_columns <- which( colnames(a14_beta_matrix)
                                  %in% a14_models_betas )
a14_models_ses_columns <- which( colnames(a14_beta_matrix)
                                %in% a14_models_ses )
a14_models_ps_columns <- which( colnames(a14_beta_matrix)
                               %in% a14_models_ps )

a14_betas_only <- cbind(a14_beta_matrix[,a14_models_betas_columns[1]],
                       a14_beta_matrix[,a14_models_betas_columns[2]],
                       a14_beta_matrix[,a14_models_betas_columns[3]])

colnames(a14_betas_only) <- c("a14_both", "a14_esc","a14_nor")

a14_ses_only <- cbind(a14_beta_matrix[,a14_models_ses_columns[1]],
                     a14_beta_matrix[,a14_models_ses_columns[2]],
                     a14_beta_matrix[,a14_models_ses_columns[3]])

colnames(a14_betas_only) <- c("a14_both", "a14_esc","a14_nor")

a14_p_only <- cbind(a14_beta_matrix[,a14_models_ps_columns[1]],
                   a14_beta_matrix[,a14_models_ps_columns[2]],
                   a14_beta_matrix[,a14_models_ps_columns[3]])

colnames(a14_p_only) <- c("a14_both", "a14_esc","a14_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a14_both",
  "a14_esc",
  "a14_nor")))

lowest_ps[which(a14_p_only<=0.05)]<-a14_p_only[which(a14_p_only<=0.05)]

lowest_ps_a14 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a14_both",
                                                            "a14_esc",
                                                            "a14_nor")))

lowest_ps_a14[which(a14_p_only<=0.05)]<-a14_p_only[which(a14_p_only<=0.05)]

heatmap.2( a14_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a14,5),notecol="black",
           main="Light-headedness"
)

##### Fit models for a15: room-spinning (adrenergic) #####

a15_models_both <- vector(mode = "list", length = length(PRSs_31))
a15_models_esc <- vector(mode = "list", length = length(PRSs_31))
a15_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a15 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a15_s <- lmer( asec15wk ~ week + week2 +
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a15_s)[-1]))==1){
    ftc_a15[length(ftc_a15)+1] <- lm_both_a15_s
    names(ftc_a15)[length(ftc_a15)] <- PRSs[i]
  }
  a15_models_both[[i]] <- lm_both_a15_s

  lm_esc_a15_s <- lmer( asec15wk ~ week + week2 +
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a15_s)[-1]))==1){
    ftc_a15[length(ftc_a15)+1] <- lm_esc_a15_s
    names(ftc_a15)[length(ftc_a15)] <- PRSs[i]
  }
  a15_models_esc[[i]] <- lm_esc_a15_s

  lm_nor_a15_s <- lmer( asec15wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a15_s)[-1]))==1){
    ftc_a15[length(ftc_a15)+1] <- lm_nor_a15_s
    names(ftc_a15)[length(ftc_a15)] <- PRSs[i]
  }
  a15_models_nor[[i]] <- lm_nor_a15_s
}  

noPRS_both_a15 <- lmer( asec15wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a15)[-1]))==1){
  ftc_a15[length(ftc_a15)+1] <- noPRS_both_a15
  names(ftc_a15)[length(ftc_a15)] <- "No_PRS"
}
a15_models_both[[31]] <- noPRS_both_a15

noPRS_esc_a15 <- lmer( asec15wk ~ week + week2 +
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a15)[-1]))==1){
  ftc_a15[length(ftc_a15)+1] <- noPRS_esc_a15
  names(ftc_a15)[length(ftc_a15)] <- "No_PRS"
}
a15_models_esc[[31]] <- noPRS_esc_a15

noPRS_nor_a15 <- lmer( asec15wk ~ week + 
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a15)[-1]))==1){
  ftc_a15[length(ftc_a15)+1] <- noPRS_nor_a15
  names(ftc_a15)[length(ftc_a15)] <- "No_PRS"
}
a15_models_nor[[31]] <- noPRS_nor_a15

names(a15_models_both) <- PRSs_31
names(a15_models_esc) <- PRSs_31
names(a15_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a15 #####

a15_models_both_betas <- vector(mode = "list", length = length(PRSs))
a15_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a15_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a15_models_both_betas[i] <- list(summary(a15_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a15_models_esc_betas[i] <- list(summary(a15_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a15_models_nor_betas[i] <- list(summary(a15_models_nor[[i]])$coefficients[5,])
}

names(a15_models_both_betas) <- PRSs
names(a15_models_esc_betas) <- PRSs
names(a15_models_nor_betas) <- PRSs

a15_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a15_models_both_betas))[,1],
    t(as.data.frame(a15_models_both_betas))[,2],
    t(as.data.frame(a15_models_both_betas))[,5],
    t(as.data.frame(a15_models_esc_betas))[,1],
    t(as.data.frame(a15_models_esc_betas))[,2],
    t(as.data.frame(a15_models_esc_betas))[,5],
    t(as.data.frame(a15_models_nor_betas))[,1],
    t(as.data.frame(a15_models_nor_betas))[,2],
    t(as.data.frame(a15_models_nor_betas))[,5]
  )
)

colnames( a15_beta_matrix ) <- c(
  "a15_both_B", "a15_both_se", "a15_both_p", 
  "a15_esc_B", "a15_esc_se", "a15_esc_p",
  "a15_nor_B", "a15_nor_se", "a15_nor_p"
)

##### Get PRSs' explained variance for a15 #####

a15_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a15_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a15_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a15_matrix_both[i,] <- r.squaredGLMM(a15_models_both[[i]])
  a15_matrix_esc[i,] <- r.squaredGLMM(a15_models_esc[[i]])
  a15_matrix_nor[i,] <- r.squaredGLMM(a15_models_nor[[i]])
}

##### Make output files for a15 #####

sink( "a15.csv")

# Betas and p-values

cat("a15_betas\n,")
write.table( a15_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a15_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a15_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a15_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a15 heatmap #####

a15_models_betas <- c("a15_both_B", "a15_esc_B", "a15_nor_B")
a15_models_ses <- c("a15_both_se", "a15_esc_se", "a15_nor_se")
a15_models_ps <- c("a15_both_p", "a15_esc_p", "a15_nor_p")

a15_models_betas_columns <- which( colnames(a15_beta_matrix)
                                  %in% a15_models_betas )
a15_models_ses_columns <- which( colnames(a15_beta_matrix)
                                %in% a15_models_ses )
a15_models_ps_columns <- which( colnames(a15_beta_matrix)
                               %in% a15_models_ps )

a15_betas_only <- cbind(a15_beta_matrix[,a15_models_betas_columns[1]],
                       a15_beta_matrix[,a15_models_betas_columns[2]],
                       a15_beta_matrix[,a15_models_betas_columns[3]])

colnames(a15_betas_only) <- c("a15_both", "a15_esc","a15_nor")

a15_ses_only <- cbind(a15_beta_matrix[,a15_models_ses_columns[1]],
                     a15_beta_matrix[,a15_models_ses_columns[2]],
                     a15_beta_matrix[,a15_models_ses_columns[3]])

colnames(a15_betas_only) <- c("a15_both", "a15_esc","a15_nor")

a15_p_only <- cbind(a15_beta_matrix[,a15_models_ps_columns[1]],
                   a15_beta_matrix[,a15_models_ps_columns[2]],
                   a15_beta_matrix[,a15_models_ps_columns[3]])

colnames(a15_p_only) <- c("a15_both", "a15_esc","a15_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a15_both",
  "a15_esc",
  "a15_nor")))

lowest_ps[which(a15_p_only<=0.05)]<-a15_p_only[which(a15_p_only<=0.05)]

lowest_ps_a15 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a15_both",
                                                            "a15_esc",
                                                            "a15_nor")))

lowest_ps_a15[which(a15_p_only<=0.05)]<-a15_p_only[which(a15_p_only<=0.05)]

heatmap.2( a15_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a15,5),notecol="black",
           main="Room-spinning"
)

##### Fit models for a16: sweating #####

a16_models_both <- vector(mode = "list", length = length(PRSs_31))
a16_models_esc <- vector(mode = "list", length = length(PRSs_31))
a16_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a16 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a16_s <- lmer( asec16wk ~ week + 
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a16_s)[-1]))==1){
    ftc_a16[length(ftc_a16)+1] <- lm_both_a16_s
    names(ftc_a16)[length(ftc_a16)] <- PRSs[i]
  }
  a16_models_both[[i]] <- lm_both_a16_s
  
  lm_esc_a16_s <- lmer( asec16wk ~ week + 
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a16_s)[-1]))==1){
    ftc_a16[length(ftc_a16)+1] <- lm_esc_a16_s
    names(ftc_a16)[length(ftc_a16)] <- PRSs[i]
  }
  a16_models_esc[[i]] <- lm_esc_a16_s
  
  lm_nor_a16_s <- lmer( asec16wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a16_s)[-1]))==1){
    ftc_a16[length(ftc_a16)+1] <- lm_nor_a16_s
    names(ftc_a16)[length(ftc_a16)] <- PRSs[i]
  }
  a16_models_nor[[i]] <- lm_nor_a16_s
}  

noPRS_both_a16 <- lmer( asec16wk ~ week + 
                          sex + cage + drug + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a16)[-1]))==1){
  ftc_a16[length(ftc_a16)+1] <- noPRS_both_a16
  names(ftc_a16)[length(ftc_a16)] <- "No_PRS"
}
a16_models_both[[31]] <- noPRS_both_a16

noPRS_esc_a16 <- lmer( asec16wk ~ week +
                         sex + cage + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a16)[-1]))==1){
  ftc_a16[length(ftc_a16)+1] <- noPRS_esc_a16
  names(ftc_a16)[length(ftc_a16)] <- "No_PRS"
}
a16_models_esc[[31]] <- noPRS_esc_a16

noPRS_nor_a16 <- lmer( asec16wk ~ week + 
                         sex + cage +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a16)[-1]))==1){
  ftc_a16[length(ftc_a16)+1] <- noPRS_nor_a16
  names(ftc_a16)[length(ftc_a16)] <- "No_PRS"
}
a16_models_nor[[31]] <- noPRS_nor_a16

names(a16_models_both) <- PRSs_31
names(a16_models_esc) <- PRSs_31
names(a16_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a16 #####

a16_models_both_betas <- vector(mode = "list", length = length(PRSs))
a16_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a16_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a16_models_both_betas[i] <- list(summary(a16_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a16_models_esc_betas[i] <- list(summary(a16_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a16_models_nor_betas[i] <- list(summary(a16_models_nor[[i]])$coefficients[5,])
}

names(a16_models_both_betas) <- PRSs
names(a16_models_esc_betas) <- PRSs
names(a16_models_nor_betas) <- PRSs

a16_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a16_models_both_betas))[,1],
    t(as.data.frame(a16_models_both_betas))[,2],
    t(as.data.frame(a16_models_both_betas))[,5],
    t(as.data.frame(a16_models_esc_betas))[,1],
    t(as.data.frame(a16_models_esc_betas))[,2],
    t(as.data.frame(a16_models_esc_betas))[,5],
    t(as.data.frame(a16_models_nor_betas))[,1],
    t(as.data.frame(a16_models_nor_betas))[,2],
    t(as.data.frame(a16_models_nor_betas))[,5]
  )
)

colnames( a16_beta_matrix ) <- c(
  "a16_both_B", "a16_both_se", "a16_both_p", 
  "a16_esc_B", "a16_esc_se", "a16_esc_p",
  "a16_nor_B", "a16_nor_se", "a16_nor_p"
)

##### Get PRSs' explained variance for a16 #####

a16_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a16_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a16_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a16_matrix_both[i,] <- r.squaredGLMM(a16_models_both[[i]])
  a16_matrix_esc[i,] <- r.squaredGLMM(a16_models_esc[[i]])
  a16_matrix_nor[i,] <- r.squaredGLMM(a16_models_nor[[i]])
}

##### Make output files for a16 #####

sink( "a16.csv")

# Betas and p-values

cat("a16_betas\n,")
write.table( a16_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a16_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a16_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a16_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a16 heatmap #####

a16_models_betas <- c("a16_both_B", "a16_esc_B", "a16_nor_B")
a16_models_ses <- c("a16_both_se", "a16_esc_se", "a16_nor_se")
a16_models_ps <- c("a16_both_p", "a16_esc_p", "a16_nor_p")

a16_models_betas_columns <- which( colnames(a16_beta_matrix)
                                   %in% a16_models_betas )
a16_models_ses_columns <- which( colnames(a16_beta_matrix)
                                 %in% a16_models_ses )
a16_models_ps_columns <- which( colnames(a16_beta_matrix)
                                %in% a16_models_ps )

a16_betas_only <- cbind(a16_beta_matrix[,a16_models_betas_columns[1]],
                        a16_beta_matrix[,a16_models_betas_columns[2]],
                        a16_beta_matrix[,a16_models_betas_columns[3]])

colnames(a16_betas_only) <- c("a16_both", "a16_esc","a16_nor")

a16_ses_only <- cbind(a16_beta_matrix[,a16_models_ses_columns[1]],
                      a16_beta_matrix[,a16_models_ses_columns[2]],
                      a16_beta_matrix[,a16_models_ses_columns[3]])

colnames(a16_betas_only) <- c("a16_both", "a16_esc","a16_nor")

a16_p_only <- cbind(a16_beta_matrix[,a16_models_ps_columns[1]],
                    a16_beta_matrix[,a16_models_ps_columns[2]],
                    a16_beta_matrix[,a16_models_ps_columns[3]])

colnames(a16_p_only) <- c("a16_both", "a16_esc","a16_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a16_both",
  "a16_esc",
  "a16_nor")))

lowest_ps[which(a16_p_only<=0.05)]<-a16_p_only[which(a16_p_only<=0.05)]

lowest_ps_a16 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a16_both",
                                                             "a16_esc",
                                                             "a16_nor")))

lowest_ps_a16[which(a16_p_only<=0.05)]<-a16_p_only[which(a16_p_only<=0.05)]

heatmap.2( a16_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a16,5),notecol="black",
           main="Sweating"
)

##### Fit models for a17: increased body temperature #####

a17_models_both <- vector(mode = "list", length = length(PRSs_31))
a17_models_esc <- vector(mode = "list", length = length(PRSs_31))
a17_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a17 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a17_s <- lmer( asec17wk ~ week + 
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a17_s)[-1]))==1){
    ftc_a17[length(ftc_a17)+1] <- lm_both_a17_s
    names(ftc_a17)[length(ftc_a17)] <- PRSs[i]
  }
  a17_models_both[[i]] <- lm_both_a17_s
  
  lm_esc_a17_s <- lmer( asec17wk ~ week + 
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a17_s)[-1]))==1){
    ftc_a17[length(ftc_a17)+1] <- lm_esc_a17_s
    names(ftc_a17)[length(ftc_a17)] <- PRSs[i]
  }
  a17_models_esc[[i]] <- lm_esc_a17_s
  
  lm_nor_a17_s <- lmer( asec17wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a17_s)[-1]))==1){
    ftc_a17[length(ftc_a17)+1] <- lm_nor_a17_s
    names(ftc_a17)[length(ftc_a17)] <- PRSs[i]
  }
  a17_models_nor[[i]] <- lm_nor_a17_s
}  

noPRS_both_a17 <- lmer( asec17wk ~ week + 
                          sex + cage + drug + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a17)[-1]))==1){
  ftc_a17[length(ftc_a17)+1] <- noPRS_both_a17
  names(ftc_a17)[length(ftc_a17)] <- "No_PRS"
}
a17_models_both[[31]] <- noPRS_both_a17

noPRS_esc_a17 <- lmer( asec17wk ~ week +
                         sex + cage + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a17)[-1]))==1){
  ftc_a17[length(ftc_a17)+1] <- noPRS_esc_a17
  names(ftc_a17)[length(ftc_a17)] <- "No_PRS"
}
a17_models_esc[[31]] <- noPRS_esc_a17

noPRS_nor_a17 <- lmer( asec17wk ~ week + 
                         sex + cage +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a17)[-1]))==1){
  ftc_a17[length(ftc_a17)+1] <- noPRS_nor_a17
  names(ftc_a17)[length(ftc_a17)] <- "No_PRS"
}
a17_models_nor[[31]] <- noPRS_nor_a17

names(a17_models_both) <- PRSs_31
names(a17_models_esc) <- PRSs_31
names(a17_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a17 #####

a17_models_both_betas <- vector(mode = "list", length = length(PRSs))
a17_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a17_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a17_models_both_betas[i] <- list(summary(a17_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a17_models_esc_betas[i] <- list(summary(a17_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a17_models_nor_betas[i] <- list(summary(a17_models_nor[[i]])$coefficients[5,])
}

names(a17_models_both_betas) <- PRSs
names(a17_models_esc_betas) <- PRSs
names(a17_models_nor_betas) <- PRSs

a17_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a17_models_both_betas))[,1],
    t(as.data.frame(a17_models_both_betas))[,2],
    t(as.data.frame(a17_models_both_betas))[,5],
    t(as.data.frame(a17_models_esc_betas))[,1],
    t(as.data.frame(a17_models_esc_betas))[,2],
    t(as.data.frame(a17_models_esc_betas))[,5],
    t(as.data.frame(a17_models_nor_betas))[,1],
    t(as.data.frame(a17_models_nor_betas))[,2],
    t(as.data.frame(a17_models_nor_betas))[,5]
  )
)

colnames( a17_beta_matrix ) <- c(
  "a17_both_B", "a17_both_se", "a17_both_p", 
  "a17_esc_B", "a17_esc_se", "a17_esc_p",
  "a17_nor_B", "a17_nor_se", "a17_nor_p"
)

##### Get PRSs' explained variance for a17 #####

a17_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a17_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a17_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a17_matrix_both[i,] <- r.squaredGLMM(a17_models_both[[i]])
  a17_matrix_esc[i,] <- r.squaredGLMM(a17_models_esc[[i]])
  a17_matrix_nor[i,] <- r.squaredGLMM(a17_models_nor[[i]])
}

##### Make output files for a17 #####

sink( "a17.csv")

# Betas and p-values

cat("a17_betas\n,")
write.table( a17_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a17_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a17_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a17_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a17 heatmap #####

a17_models_betas <- c("a17_both_B", "a17_esc_B", "a17_nor_B")
a17_models_ses <- c("a17_both_se", "a17_esc_se", "a17_nor_se")
a17_models_ps <- c("a17_both_p", "a17_esc_p", "a17_nor_p")

a17_models_betas_columns <- which( colnames(a17_beta_matrix)
                                   %in% a17_models_betas )
a17_models_ses_columns <- which( colnames(a17_beta_matrix)
                                 %in% a17_models_ses )
a17_models_ps_columns <- which( colnames(a17_beta_matrix)
                                %in% a17_models_ps )

a17_betas_only <- cbind(a17_beta_matrix[,a17_models_betas_columns[1]],
                        a17_beta_matrix[,a17_models_betas_columns[2]],
                        a17_beta_matrix[,a17_models_betas_columns[3]])

colnames(a17_betas_only) <- c("a17_both", "a17_esc","a17_nor")

a17_ses_only <- cbind(a17_beta_matrix[,a17_models_ses_columns[1]],
                      a17_beta_matrix[,a17_models_ses_columns[2]],
                      a17_beta_matrix[,a17_models_ses_columns[3]])

colnames(a17_betas_only) <- c("a17_both", "a17_esc","a17_nor")

a17_p_only <- cbind(a17_beta_matrix[,a17_models_ps_columns[1]],
                    a17_beta_matrix[,a17_models_ps_columns[2]],
                    a17_beta_matrix[,a17_models_ps_columns[3]])

colnames(a17_p_only) <- c("a17_both", "a17_esc","a17_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a17_both",
  "a17_esc",
  "a17_nor")))

lowest_ps[which(a17_p_only<=0.05)]<-a17_p_only[which(a17_p_only<=0.05)]

lowest_ps_a17 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a17_both",
                                                             "a17_esc",
                                                             "a17_nor")))

lowest_ps_a17[which(a17_p_only<=0.05)]<-a17_p_only[which(a17_p_only<=0.05)]

heatmap.2( a17_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a17,5),notecol="black",
           main="Increased body temperature"
)

##### Fit models for a18: tremor #####

a18_models_both <- vector(mode = "list", length = length(PRSs_31))
a18_models_esc <- vector(mode = "list", length = length(PRSs_31))
a18_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a18 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a18_s <- lmer( asec18wk ~ week + 
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a18_s)[-1]))==1){
    ftc_a18[length(ftc_a18)+1] <- lm_both_a18_s
    names(ftc_a18)[length(ftc_a18)] <- PRSs[i]
  }
  a18_models_both[[i]] <- lm_both_a18_s
  
  lm_esc_a18_s <- lmer( asec18wk ~ week + 
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a18_s)[-1]))==1){
    ftc_a18[length(ftc_a18)+1] <- lm_esc_a18_s
    names(ftc_a18)[length(ftc_a18)] <- PRSs[i]
  }
  a18_models_esc[[i]] <- lm_esc_a18_s
  
  lm_nor_a18_s <- lmer( asec18wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a18_s)[-1]))==1){
    ftc_a18[length(ftc_a18)+1] <- lm_nor_a18_s
    names(ftc_a18)[length(ftc_a18)] <- PRSs[i]
  }
  a18_models_nor[[i]] <- lm_nor_a18_s
}  

noPRS_both_a18 <- lmer( asec18wk ~ week + 
                          sex + cage + drug + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a18)[-1]))==1){
  ftc_a18[length(ftc_a18)+1] <- noPRS_both_a18
  names(ftc_a18)[length(ftc_a18)] <- "No_PRS"
}
a18_models_both[[31]] <- noPRS_both_a18

noPRS_esc_a18 <- lmer( asec18wk ~ week +
                         sex + cage + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a18)[-1]))==1){
  ftc_a18[length(ftc_a18)+1] <- noPRS_esc_a18
  names(ftc_a18)[length(ftc_a18)] <- "No_PRS"
}
a18_models_esc[[31]] <- noPRS_esc_a18

noPRS_nor_a18 <- lmer( asec18wk ~ week + 
                         sex + cage +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a18)[-1]))==1){
  ftc_a18[length(ftc_a18)+1] <- noPRS_nor_a18
  names(ftc_a18)[length(ftc_a18)] <- "No_PRS"
}
a18_models_nor[[31]] <- noPRS_nor_a18

names(a18_models_both) <- PRSs_31
names(a18_models_esc) <- PRSs_31
names(a18_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a18 #####

a18_models_both_betas <- vector(mode = "list", length = length(PRSs))
a18_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a18_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a18_models_both_betas[i] <- list(summary(a18_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a18_models_esc_betas[i] <- list(summary(a18_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a18_models_nor_betas[i] <- list(summary(a18_models_nor[[i]])$coefficients[5,])
}

names(a18_models_both_betas) <- PRSs
names(a18_models_esc_betas) <- PRSs
names(a18_models_nor_betas) <- PRSs

a18_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a18_models_both_betas))[,1],
    t(as.data.frame(a18_models_both_betas))[,2],
    t(as.data.frame(a18_models_both_betas))[,5],
    t(as.data.frame(a18_models_esc_betas))[,1],
    t(as.data.frame(a18_models_esc_betas))[,2],
    t(as.data.frame(a18_models_esc_betas))[,5],
    t(as.data.frame(a18_models_nor_betas))[,1],
    t(as.data.frame(a18_models_nor_betas))[,2],
    t(as.data.frame(a18_models_nor_betas))[,5]
  )
)

colnames( a18_beta_matrix ) <- c(
  "a18_both_B", "a18_both_se", "a18_both_p", 
  "a18_esc_B", "a18_esc_se", "a18_esc_p",
  "a18_nor_B", "a18_nor_se", "a18_nor_p"
)

##### Get PRSs' explained variance for a18 #####

a18_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a18_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a18_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a18_matrix_both[i,] <- r.squaredGLMM(a18_models_both[[i]])
  a18_matrix_esc[i,] <- r.squaredGLMM(a18_models_esc[[i]])
  a18_matrix_nor[i,] <- r.squaredGLMM(a18_models_nor[[i]])
}

##### Make output files for a18 #####

sink( "a18.csv")

# Betas and p-values

cat("a18_betas\n,")
write.table( a18_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a18_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a18_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a18_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a18 heatmap #####

a18_models_betas <- c("a18_both_B", "a18_esc_B", "a18_nor_B")
a18_models_ses <- c("a18_both_se", "a18_esc_se", "a18_nor_se")
a18_models_ps <- c("a18_both_p", "a18_esc_p", "a18_nor_p")

a18_models_betas_columns <- which( colnames(a18_beta_matrix)
                                   %in% a18_models_betas )
a18_models_ses_columns <- which( colnames(a18_beta_matrix)
                                 %in% a18_models_ses )
a18_models_ps_columns <- which( colnames(a18_beta_matrix)
                                %in% a18_models_ps )

a18_betas_only <- cbind(a18_beta_matrix[,a18_models_betas_columns[1]],
                        a18_beta_matrix[,a18_models_betas_columns[2]],
                        a18_beta_matrix[,a18_models_betas_columns[3]])

colnames(a18_betas_only) <- c("a18_both", "a18_esc","a18_nor")

a18_ses_only <- cbind(a18_beta_matrix[,a18_models_ses_columns[1]],
                      a18_beta_matrix[,a18_models_ses_columns[2]],
                      a18_beta_matrix[,a18_models_ses_columns[3]])

colnames(a18_betas_only) <- c("a18_both", "a18_esc","a18_nor")

a18_p_only <- cbind(a18_beta_matrix[,a18_models_ps_columns[1]],
                    a18_beta_matrix[,a18_models_ps_columns[2]],
                    a18_beta_matrix[,a18_models_ps_columns[3]])

colnames(a18_p_only) <- c("a18_both", "a18_esc","a18_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a18_both",
  "a18_esc",
  "a18_nor")))

lowest_ps[which(a18_p_only<=0.05)]<-a18_p_only[which(a18_p_only<=0.05)]

lowest_ps_a18 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a18_both",
                                                             "a18_esc",
                                                             "a18_nor")))

lowest_ps_a18[which(a18_p_only<=0.05)]<-a18_p_only[which(a18_p_only<=0.05)]

heatmap.2( a18_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a18,5),notecol="black",
           main="Tremor"
)
##### Fit models for a19: disorientation #####

a19_models_both <- vector(mode = "list", length = length(PRSs_31))
a19_models_esc <- vector(mode = "list", length = length(PRSs_31))
a19_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a19 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a19_s <- lmer( asec19wk ~ week + week2 +
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a19_s)[-1]))==1){
    ftc_a19[length(ftc_a19)+1] <- lm_both_a19_s
    names(ftc_a19)[length(ftc_a19)] <- PRSs[i]
  }
  a19_models_both[[i]] <- lm_both_a19_s

  lm_esc_a19_s <- lmer( asec19wk ~ week + week2 +
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a19_s)[-1]))==1){
    ftc_a19[length(ftc_a19)+1] <- lm_esc_a19_s
    names(ftc_a19)[length(ftc_a19)] <- PRSs[i]
  }
  a19_models_esc[[i]] <- lm_esc_a19_s

  lm_nor_a19_s <- lmer( asec19wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a19_s)[-1]))==1){
    ftc_a19[length(ftc_a19)+1] <- lm_nor_a19_s
    names(ftc_a19)[length(ftc_a19)] <- PRSs[i]
  }
  a19_models_nor[[i]] <- lm_nor_a19_s
}  

noPRS_both_a19 <- lmer( asec19wk ~ week + week2 +
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a19)[-1]))==1){
  ftc_a19[length(ftc_a19)+1] <- noPRS_both_a19
  names(ftc_a19)[length(ftc_a19)] <- "No_PRS"
}
a19_models_both[[31]] <- noPRS_both_a19

noPRS_esc_a19 <- lmer( asec19wk ~ week + week2 +
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a19)[-1]))==1){
  ftc_a19[length(ftc_a19)+1] <- noPRS_esc_a19
  names(ftc_a19)[length(ftc_a19)] <- "No_PRS"
}
a19_models_esc[[31]] <- noPRS_esc_a19

noPRS_nor_a19 <- lmer( asec19wk ~ week + week2 +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a19)[-1]))==1){
  ftc_a19[length(ftc_a19)+1] <- noPRS_nor_a19
  names(ftc_a19)[length(ftc_a19)] <- "No_PRS"
}
a19_models_nor[[31]] <- noPRS_nor_a19

names(a19_models_both) <- PRSs_31
names(a19_models_esc) <- PRSs_31
names(a19_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a19 #####

a19_models_both_betas <- vector(mode = "list", length = length(PRSs))
a19_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a19_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a19_models_both_betas[i] <- list(summary(a19_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a19_models_esc_betas[i] <- list(summary(a19_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a19_models_nor_betas[i] <- list(summary(a19_models_nor[[i]])$coefficients[5,])
}

names(a19_models_both_betas) <- PRSs
names(a19_models_esc_betas) <- PRSs
names(a19_models_nor_betas) <- PRSs

a19_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a19_models_both_betas))[,1],
    t(as.data.frame(a19_models_both_betas))[,2],
    t(as.data.frame(a19_models_both_betas))[,5],
    t(as.data.frame(a19_models_esc_betas))[,1],
    t(as.data.frame(a19_models_esc_betas))[,2],
    t(as.data.frame(a19_models_esc_betas))[,5],
    t(as.data.frame(a19_models_nor_betas))[,1],
    t(as.data.frame(a19_models_nor_betas))[,2],
    t(as.data.frame(a19_models_nor_betas))[,5]
  )
)

colnames( a19_beta_matrix ) <- c(
  "a19_both_B", "a19_both_se", "a19_both_p", 
  "a19_esc_B", "a19_esc_se", "a19_esc_p",
  "a19_nor_B", "a19_nor_se", "a19_nor_p"
)

##### Get PRSs' explained variance for a19 #####

a19_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a19_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a19_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a19_matrix_both[i,] <- r.squaredGLMM(a19_models_both[[i]])
  a19_matrix_esc[i,] <- r.squaredGLMM(a19_models_esc[[i]])
  a19_matrix_nor[i,] <- r.squaredGLMM(a19_models_nor[[i]])
}

##### Make output files for a19 #####

sink( "a19.csv")

# Betas and p-values

cat("a19_betas\n,")
write.table( a19_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a19_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a19_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a19_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a19 heatmap #####

a19_models_betas <- c("a19_both_B", "a19_esc_B", "a19_nor_B")
a19_models_ses <- c("a19_both_se", "a19_esc_se", "a19_nor_se")
a19_models_ps <- c("a19_both_p", "a19_esc_p", "a19_nor_p")

a19_models_betas_columns <- which( colnames(a19_beta_matrix)
                                  %in% a19_models_betas )
a19_models_ses_columns <- which( colnames(a19_beta_matrix)
                                %in% a19_models_ses )
a19_models_ps_columns <- which( colnames(a19_beta_matrix)
                               %in% a19_models_ps )

a19_betas_only <- cbind(a19_beta_matrix[,a19_models_betas_columns[1]],
                       a19_beta_matrix[,a19_models_betas_columns[2]],
                       a19_beta_matrix[,a19_models_betas_columns[3]])

colnames(a19_betas_only) <- c("a19_both", "a19_esc","a19_nor")

a19_ses_only <- cbind(a19_beta_matrix[,a19_models_ses_columns[1]],
                     a19_beta_matrix[,a19_models_ses_columns[2]],
                     a19_beta_matrix[,a19_models_ses_columns[3]])

colnames(a19_betas_only) <- c("a19_both", "a19_esc","a19_nor")

a19_p_only <- cbind(a19_beta_matrix[,a19_models_ps_columns[1]],
                   a19_beta_matrix[,a19_models_ps_columns[2]],
                   a19_beta_matrix[,a19_models_ps_columns[3]])

colnames(a19_p_only) <- c("a19_both", "a19_esc","a19_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a19_both",
  "a19_esc",
  "a19_nor")))

lowest_ps[which(a19_p_only<=0.05)]<-a19_p_only[which(a19_p_only<=0.05)]

lowest_ps_a19 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a19_both",
                                                            "a19_esc",
                                                            "a19_nor")))

lowest_ps_a19[which(a19_p_only<=0.05)]<-a19_p_only[which(a19_p_only<=0.05)]

heatmap.2( a19_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a19,5),notecol="black",
           main="Disorientation"
)

##### Fit models for a20: yawning #####

a20_models_both <- vector(mode = "list", length = length(PRSs_31))
a20_models_esc <- vector(mode = "list", length = length(PRSs_31))
a20_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a20 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a20_s <- lmer( asec20wk ~ week + 
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a20_s)[-1]))==1){
    ftc_a20[length(ftc_a20)+1] <- lm_both_a20_s
    names(ftc_a20)[length(ftc_a20)] <- PRSs[i]
  }
  a20_models_both[[i]] <- lm_both_a20_s

  lm_esc_a20_s <- lmer( asec20wk ~ week + 
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a20_s)[-1]))==1){
    ftc_a20[length(ftc_a20)+1] <- lm_esc_a20_s
    names(ftc_a20)[length(ftc_a20)] <- PRSs[i]
  }
  a20_models_esc[[i]] <- lm_esc_a20_s

  lm_nor_a20_s <- lmer( asec20wk ~ week + 
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a20_s)[-1]))==1){
    ftc_a20[length(ftc_a20)+1] <- lm_nor_a20_s
    names(ftc_a20)[length(ftc_a20)] <- PRSs[i]
  }
  a20_models_nor[[i]] <- lm_nor_a20_s
  
}  

noPRS_both_a20 <- lmer( asec20wk ~ week + 
                         sex + cage + drug + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a20)[-1]))==1){
  ftc_a20[length(ftc_a20)+1] <- noPRS_both_a20
  names(ftc_a20)[length(ftc_a20)] <- "No_PRS"
}
a20_models_both[[31]] <- noPRS_both_a20

noPRS_esc_a20 <- lmer( asec20wk ~ week + 
                        sex + cage + 
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a20)[-1]))==1){
  ftc_a20[length(ftc_a20)+1] <- noPRS_esc_a20
  names(ftc_a20)[length(ftc_a20)] <- "No_PRS"
}
a20_models_esc[[31]] <- noPRS_esc_a20

noPRS_nor_a20 <- lmer( asec20wk ~ week +
                        sex + cage +
                        zPC1 + zPC2 + zPC3 + 
                        (1|centreid) + (1|subjectid),
                      data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a20)[-1]))==1){
  ftc_a20[length(ftc_a20)+1] <- noPRS_nor_a20
  names(ftc_a20)[length(ftc_a20)] <- "No_PRS"
}
a20_models_nor[[31]] <- noPRS_nor_a20

names(a20_models_both) <- PRSs_31
names(a20_models_esc) <- PRSs_31
names(a20_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a20 #####

a20_models_both_betas <- vector(mode = "list", length = length(PRSs))
a20_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a20_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a20_models_both_betas[i] <- list(summary(a20_models_both[[i]])$coefficients[6,])
}
for (i in 1:30){
  a20_models_esc_betas[i] <- list(summary(a20_models_esc[[i]])$coefficients[5,])
}
for (i in 1:30){
  a20_models_nor_betas[i] <- list(summary(a20_models_nor[[i]])$coefficients[5,])
}

names(a20_models_both_betas) <- PRSs
names(a20_models_esc_betas) <- PRSs
names(a20_models_nor_betas) <- PRSs

a20_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a20_models_both_betas))[,1],
    t(as.data.frame(a20_models_both_betas))[,2],
    t(as.data.frame(a20_models_both_betas))[,5],
    t(as.data.frame(a20_models_esc_betas))[,1],
    t(as.data.frame(a20_models_esc_betas))[,2],
    t(as.data.frame(a20_models_esc_betas))[,5],
    t(as.data.frame(a20_models_nor_betas))[,1],
    t(as.data.frame(a20_models_nor_betas))[,2],
    t(as.data.frame(a20_models_nor_betas))[,5]
  )
)

colnames( a20_beta_matrix ) <- c(
  "a20_both_B", "a20_both_se", "a20_both_p", 
  "a20_esc_B", "a20_esc_se", "a20_esc_p",
  "a20_nor_B", "a20_nor_se", "a20_nor_p"
)

##### Get PRSs' explained variance for a20 #####

a20_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a20_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a20_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a20_matrix_both[i,] <- r.squaredGLMM(a20_models_both[[i]])
  a20_matrix_esc[i,] <- r.squaredGLMM(a20_models_esc[[i]])
  a20_matrix_nor[i,] <- r.squaredGLMM(a20_models_nor[[i]])
}

##### Make output files for a20 #####

sink( "a20.csv")

# Betas and p-values

cat("a20_betas\n,")
write.table( a20_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a20_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a20_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a20_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a20 heatmap #####

a20_models_betas <- c("a20_both_B", "a20_esc_B", "a20_nor_B")
a20_models_ses <- c("a20_both_se", "a20_esc_se", "a20_nor_se")
a20_models_ps <- c("a20_both_p", "a20_esc_p", "a20_nor_p")

a20_models_betas_columns <- which( colnames(a20_beta_matrix)
                                  %in% a20_models_betas )
a20_models_ses_columns <- which( colnames(a20_beta_matrix)
                                %in% a20_models_ses )
a20_models_ps_columns <- which( colnames(a20_beta_matrix)
                               %in% a20_models_ps )

a20_betas_only <- cbind(a20_beta_matrix[,a20_models_betas_columns[1]],
                       a20_beta_matrix[,a20_models_betas_columns[2]],
                       a20_beta_matrix[,a20_models_betas_columns[3]])

colnames(a20_betas_only) <- c("a20_both", "a20_esc","a20_nor")

a20_ses_only <- cbind(a20_beta_matrix[,a20_models_ses_columns[1]],
                     a20_beta_matrix[,a20_models_ses_columns[2]],
                     a20_beta_matrix[,a20_models_ses_columns[3]])

colnames(a20_betas_only) <- c("a20_both", "a20_esc","a20_nor")

a20_p_only <- cbind(a20_beta_matrix[,a20_models_ps_columns[1]],
                   a20_beta_matrix[,a20_models_ps_columns[2]],
                   a20_beta_matrix[,a20_models_ps_columns[3]])

colnames(a20_p_only) <- c("a20_both", "a20_esc","a20_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a20_both",
  "a20_esc",
  "a20_nor")))

lowest_ps[which(a20_p_only<=0.05)]<-a20_p_only[which(a20_p_only<=0.05)]

lowest_ps_a20 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a20_both",
                                                            "a20_esc",
                                                            "a20_nor")))

lowest_ps_a20[which(a20_p_only<=0.05)]<-a20_p_only[which(a20_p_only<=0.05)]

heatmap.2( a20_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a20,5),notecol="black",
           main="Yawning"
)

##### Fit models for a21: weight gain (histaminergic) #####

a21_models_both <- vector(mode = "list", length = length(PRSs_31))
a21_models_esc <- vector(mode = "list", length = length(PRSs_31))
a21_models_nor <- vector(mode = "list", length = length(PRSs_31))

ftc_a21 <- vector(mode = "list", length = 0)

for ( i in 1:length( PRSs ) ){
  
  lm_both_a21_s <- lmer( asec21wk ~ week + week2 +
                           sex + cage + drug + GENDEP[,PRSs[i]] +
                           zPC1 + zPC2 + zPC3 + 
                           (1|centreid) + (1|subjectid),
                         data = GENDEP, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_both_a21_s)[-1]))==1){
    ftc_a21[length(ftc_a21)+1] <- lm_both_a21_s
    names(ftc_a21)[length(ftc_a21)] <- PRSs[i]
  }
  a21_models_both[[i]] <- lm_both_a21_s
  
  lm_esc_a21_s <- lmer( asec21wk ~ week + week2 +
                          sex + cage + GENDEP_escit[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_escit, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_esc_a21_s)[-1]))==1){
    ftc_a21[length(ftc_a21)+1] <- lm_esc_a21_s
    names(ftc_a21)[length(ftc_a21)] <- PRSs[i]
  }
  a21_models_esc[[i]] <- lm_esc_a21_s
  
  lm_nor_a21_s <- lmer( asec21wk ~ week + week2 +
                          sex + cage + GENDEP_nortrip[,PRSs[i]] +
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP_nortrip, REML=FALSE )
  if (length(grep("warnings",capture.output(lm_nor_a21_s)[-1]))==1){
    ftc_a21[length(ftc_a21)+1] <- lm_nor_a21_s
    names(ftc_a21)[length(ftc_a21)] <- PRSs[i]
  }
  a21_models_nor[[i]] <- lm_nor_a21_s
  
}  

noPRS_both_a21 <- lmer( asec21wk ~ week + week2 +
                          sex + cage + drug + 
                          zPC1 + zPC2 + zPC3 + 
                          (1|centreid) + (1|subjectid),
                        data = GENDEP, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_both_a21)[-1]))==1){
  ftc_a21[length(ftc_a21)+1] <- noPRS_both_a21
  names(ftc_a21)[length(ftc_a21)] <- "No_PRS"
}
a21_models_both[[31]] <- noPRS_both_a21

noPRS_esc_a21 <- lmer( asec21wk ~ week + week2 +
                         sex + cage + 
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_escit, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_esc_a21)[-1]))==1){
  ftc_a21[length(ftc_a21)+1] <- noPRS_esc_a21
  names(ftc_a21)[length(ftc_a21)] <- "No_PRS"
}
a21_models_esc[[31]] <- noPRS_esc_a21

noPRS_nor_a21 <- lmer( asec21wk ~ week + week2 +
                         sex + cage +
                         zPC1 + zPC2 + zPC3 + 
                         (1|centreid) + (1|subjectid),
                       data = GENDEP_nortrip, REML=FALSE )
if (length(grep("warnings",capture.output(noPRS_nor_a21)[-1]))==1){
  ftc_a21[length(ftc_a21)+1] <- noPRS_nor_a21
  names(ftc_a21)[length(ftc_a21)] <- "No_PRS"
}
a21_models_nor[[31]] <- noPRS_nor_a21

names(a21_models_both) <- PRSs_31
names(a21_models_esc) <- PRSs_31
names(a21_models_nor) <- PRSs_31

##### Get betas, standard errors, and p-values for a21 #####

a21_models_both_betas <- vector(mode = "list", length = length(PRSs))
a21_models_esc_betas <- vector(mode = "list", length = length(PRSs))
a21_models_nor_betas <- vector(mode = "list", length = length(PRSs))

for (i in 1:30){
  a21_models_both_betas[i] <- list(summary(a21_models_both[[i]])$coefficients[7,])
}
for (i in 1:30){
  a21_models_esc_betas[i] <- list(summary(a21_models_esc[[i]])$coefficients[6,])
}
for (i in 1:30){
  a21_models_nor_betas[i] <- list(summary(a21_models_nor[[i]])$coefficients[6,])
}

names(a21_models_both_betas) <- PRSs
names(a21_models_esc_betas) <- PRSs
names(a21_models_nor_betas) <- PRSs

a21_beta_matrix <- as.matrix(
  cbind(
    t(as.data.frame(a21_models_both_betas))[,1],
    t(as.data.frame(a21_models_both_betas))[,2],
    t(as.data.frame(a21_models_both_betas))[,5],
    t(as.data.frame(a21_models_esc_betas))[,1],
    t(as.data.frame(a21_models_esc_betas))[,2],
    t(as.data.frame(a21_models_esc_betas))[,5],
    t(as.data.frame(a21_models_nor_betas))[,1],
    t(as.data.frame(a21_models_nor_betas))[,2],
    t(as.data.frame(a21_models_nor_betas))[,5]
  )
)

colnames( a21_beta_matrix ) <- c(
  "a21_both_B", "a21_both_se", "a21_both_p", 
  "a21_esc_B", "a21_esc_se", "a21_esc_p",
  "a21_nor_B", "a21_nor_se", "a21_nor_p"
)

##### Get PRSs' explained variance for a21 #####

a21_matrix_both <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a21_matrix_esc <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

a21_matrix_nor <- matrix(
  nrow=31,
  ncol=2,dimnames=list(c(
    PRSs,
    "No_PRS"
  ),c("R2m","R2c")))

for (i in 1:31){
  a21_matrix_both[i,] <- r.squaredGLMM(a21_models_both[[i]])
  a21_matrix_esc[i,] <- r.squaredGLMM(a21_models_esc[[i]])
  a21_matrix_nor[i,] <- r.squaredGLMM(a21_models_nor[[i]])
}

##### Make output files for a21 #####

sink( "a21.csv")

# Betas and p-values

cat("a21_betas\n,")
write.table( a21_beta_matrix, col.names=TRUE, sep="," )
cat( "\n" )

# Explained variance tables

cat("Both\n,")
write.table(a21_matrix_both,
            col.names=TRUE, sep=",")
cat("\nEscitalopram\n,")
write.table(a21_matrix_esc,
            col.names=TRUE, sep=",")
cat("\nNortriptyline\n,")
write.table(a21_matrix_nor,
            col.names=TRUE, sep=",")

sink()

##### Make a21 heatmap #####

a21_models_betas <- c("a21_both_B", "a21_esc_B", "a21_nor_B")
a21_models_ses <- c("a21_both_se", "a21_esc_se", "a21_nor_se")
a21_models_ps <- c("a21_both_p", "a21_esc_p", "a21_nor_p")

a21_models_betas_columns <- which( colnames(a21_beta_matrix)
                                   %in% a21_models_betas )
a21_models_ses_columns <- which( colnames(a21_beta_matrix)
                                 %in% a21_models_ses )
a21_models_ps_columns <- which( colnames(a21_beta_matrix)
                                %in% a21_models_ps )

a21_betas_only <- cbind(a21_beta_matrix[,a21_models_betas_columns[1]],
                        a21_beta_matrix[,a21_models_betas_columns[2]],
                        a21_beta_matrix[,a21_models_betas_columns[3]])

colnames(a21_betas_only) <- c("a21_both", "a21_esc","a21_nor")

a21_ses_only <- cbind(a21_beta_matrix[,a21_models_ses_columns[1]],
                      a21_beta_matrix[,a21_models_ses_columns[2]],
                      a21_beta_matrix[,a21_models_ses_columns[3]])

colnames(a21_betas_only) <- c("a21_both", "a21_esc","a21_nor")

a21_p_only <- cbind(a21_beta_matrix[,a21_models_ps_columns[1]],
                    a21_beta_matrix[,a21_models_ps_columns[2]],
                    a21_beta_matrix[,a21_models_ps_columns[3]])

colnames(a21_p_only) <- c("a21_both", "a21_esc","a21_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c(
  "a21_both",
  "a21_esc",
  "a21_nor")))

lowest_ps[which(a21_p_only<=0.05)]<-a21_p_only[which(a21_p_only<=0.05)]

lowest_ps_a21 <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "a21_both",
                                                             "a21_esc",
                                                             "a21_nor")))

lowest_ps_a21[which(a21_p_only<=0.05)]<-a21_p_only[which(a21_p_only<=0.05)]

heatmap.2( a21_betas_only, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps_a21,5),notecol="black",
           main="Weight gain"
)


##### ANALYSE EFFECTS OF POLYGENIC SCORE ON INDIVIDUAL SIDE EFFECTS (ASEC 1-21): random only #####

# This was achieved by changing the GENDEP dataset to include just the randomized participants 
# and re-running the code, with differently named csv files, rather than by copying and pasting the
# code in again with different datasets

##### Make deciles and quintiles for each PRS #####

# Sort #####

GEND_A1 <- unique(GENDEP[order(GENDEP$AZ_0_0001),]$subjectid)
GEND_A2 <- unique(GENDEP[order(GENDEP$AZ_0_01),]$subjectid)
GEND_A3 <- unique(GENDEP[order(GENDEP$AZ_0_05),]$subjectid)
GEND_A4 <- unique(GENDEP[order(GENDEP$AZ_0_1),]$subjectid)
GEND_A5 <- unique(GENDEP[order(GENDEP$AZ_0_5),]$subjectid)
GEND_A6 <- unique(GENDEP[order(GENDEP$AZ_1),]$subjectid)
GEND_C1 <- unique(GENDEP[order(GENDEP$CZ_0_0001),]$subjectid)
GEND_C2 <- unique(GENDEP[order(GENDEP$CZ_0_01),]$subjectid)
GEND_C3 <- unique(GENDEP[order(GENDEP$CZ_0_05),]$subjectid)
GEND_C4 <- unique(GENDEP[order(GENDEP$CZ_0_1),]$subjectid)
GEND_C5 <- unique(GENDEP[order(GENDEP$CZ_0_5),]$subjectid)
GEND_C6 <- unique(GENDEP[order(GENDEP$CZ_1),]$subjectid)
GEND_N1 <- unique(GENDEP[order(GENDEP$NZ_0_0001),]$subjectid)
GEND_N2 <- unique(GENDEP[order(GENDEP$NZ_0_01),]$subjectid)
GEND_N3 <- unique(GENDEP[order(GENDEP$NZ_0_05),]$subjectid)
GEND_N4 <- unique(GENDEP[order(GENDEP$NZ_0_1),]$subjectid)
GEND_N5 <- unique(GENDEP[order(GENDEP$NZ_0_5),]$subjectid)
GEND_N6 <- unique(GENDEP[order(GENDEP$NZ_1),]$subjectid)
GEND_S1 <- unique(GENDEP[order(GENDEP$SZ_0_0001),]$subjectid)
GEND_S2 <- unique(GENDEP[order(GENDEP$SZ_0_01),]$subjectid)
GEND_S3 <- unique(GENDEP[order(GENDEP$SZ_0_05),]$subjectid)
GEND_S4 <- unique(GENDEP[order(GENDEP$SZ_0_1),]$subjectid)
GEND_S5 <- unique(GENDEP[order(GENDEP$SZ_0_5),]$subjectid)
GEND_S6 <- unique(GENDEP[order(GENDEP$SZ_1),]$subjectid)
GEND_H1 <- unique(GENDEP[order(GENDEP$HZ_0_0001),]$subjectid)
GEND_H2 <- unique(GENDEP[order(GENDEP$HZ_0_01),]$subjectid)
GEND_H3 <- unique(GENDEP[order(GENDEP$HZ_0_05),]$subjectid)
GEND_H4 <- unique(GENDEP[order(GENDEP$HZ_0_1),]$subjectid)
GEND_H5 <- unique(GENDEP[order(GENDEP$HZ_0_5),]$subjectid)
GEND_H6 <- unique(GENDEP[order(GENDEP$HZ_1),]$subjectid)

# Make decile and quintile vectors #####

# deciles <- c( rep( "01", 65), rep( "02", 65), rep( "03", 65), rep( "04", 65),
#               rep( "05", 64), rep( "06", 64), rep( "07", 64), rep( "08", 65),
#               rep( "09", 65), rep( "10", 65)
#               )

quintiles <- c( rep( "01", 129), rep( "02", 130), rep( "03", 129), rep( "04", 130),
              rep( "05", 129)
)

# Make deciles for each PRS #####
                                           
# GEND_A1 <- cbind( as.character( GEND_A1), deciles )
# GEND_A2 <- cbind( as.character( GEND_A2), deciles )
# GEND_A3 <- cbind( as.character( GEND_A3), deciles )
# GEND_A4 <- cbind( as.character( GEND_A4), deciles )
# GEND_A5 <- cbind( as.character( GEND_A5), deciles )
# GEND_A6 <- cbind( as.character( GEND_A6), deciles )
# GEND_C1 <- cbind( as.character( GEND_C1), deciles )
# GEND_C2 <- cbind( as.character( GEND_C2), deciles )
# GEND_C3 <- cbind( as.character( GEND_C3), deciles )
# GEND_C4 <- cbind( as.character( GEND_C4), deciles )
# GEND_C5 <- cbind( as.character( GEND_C5), deciles )
# GEND_C6 <- cbind( as.character( GEND_C6), deciles )
# GEND_N1 <- cbind( as.character( GEND_N1), deciles )
# GEND_N2 <- cbind( as.character( GEND_N2), deciles )
# GEND_N3 <- cbind( as.character( GEND_N3), deciles )
# GEND_N4 <- cbind( as.character( GEND_N4), deciles )
# GEND_N5 <- cbind( as.character( GEND_N5), deciles )
# GEND_N6 <- cbind( as.character( GEND_N6), deciles )
# GEND_S1 <- cbind( as.character( GEND_S1), deciles )
# GEND_S2 <- cbind( as.character( GEND_S2), deciles )
# GEND_S3 <- cbind( as.character( GEND_S3), deciles )
# GEND_S4 <- cbind( as.character( GEND_S4), deciles )
# GEND_S5 <- cbind( as.character( GEND_S5), deciles )
# GEND_S6 <- cbind( as.character( GEND_S6), deciles )
# GEND_H1 <- cbind( as.character( GEND_H1), deciles )
# GEND_H2 <- cbind( as.character( GEND_H2), deciles )
# GEND_H3 <- cbind( as.character( GEND_H3), deciles )
# GEND_H4 <- cbind( as.character( GEND_H4), deciles )
# GEND_H5 <- cbind( as.character( GEND_H5), deciles )
# GEND_H6 <- cbind( as.character( GEND_H6), deciles )
# 
# colnames( GEND_A1 ) <- c( "subjectid", "A1_decile" )
# colnames( GEND_A2 ) <- c( "subjectid", "A2_decile" )
# colnames( GEND_A3 ) <- c( "subjectid", "A3_decile" )
# colnames( GEND_A4 ) <- c( "subjectid", "A4_decile" )
# colnames( GEND_A5 ) <- c( "subjectid", "A5_decile" )
# colnames( GEND_A6 ) <- c( "subjectid", "A6_decile" )
# 
# colnames( GEND_C1 ) <- c( "subjectid", "C1_decile" )
# colnames( GEND_C2 ) <- c( "subjectid", "C2_decile" )
# colnames( GEND_C3 ) <- c( "subjectid", "C3_decile" )
# colnames( GEND_C4 ) <- c( "subjectid", "C4_decile" )
# colnames( GEND_C5 ) <- c( "subjectid", "C5_decile" )
# colnames( GEND_C6 ) <- c( "subjectid", "C6_decile" )
# 
# colnames( GEND_N1 ) <- c( "subjectid", "N1_decile" )
# colnames( GEND_N2 ) <- c( "subjectid", "N2_decile" )
# colnames( GEND_N3 ) <- c( "subjectid", "N3_decile" )
# colnames( GEND_N4 ) <- c( "subjectid", "N4_decile" )
# colnames( GEND_N5 ) <- c( "subjectid", "N5_decile" )
# colnames( GEND_N6 ) <- c( "subjectid", "N6_decile" )
# 
# colnames( GEND_S1 ) <- c( "subjectid", "S1_decile" )
# colnames( GEND_S2 ) <- c( "subjectid", "S2_decile" )
# colnames( GEND_S3 ) <- c( "subjectid", "S3_decile" )
# colnames( GEND_S4 ) <- c( "subjectid", "S4_decile" )
# colnames( GEND_S5 ) <- c( "subjectid", "S5_decile" )
# colnames( GEND_S6 ) <- c( "subjectid", "S6_decile" )
# 
# colnames( GEND_H1 ) <- c( "subjectid", "H1_decile" )
# colnames( GEND_H2 ) <- c( "subjectid", "H2_decile" )
# colnames( GEND_H3 ) <- c( "subjectid", "H3_decile" )
# colnames( GEND_H4 ) <- c( "subjectid", "H4_decile" )
# colnames( GEND_H5 ) <- c( "subjectid", "H5_decile" )
# colnames( GEND_H6 ) <- c( "subjectid", "H6_decile" )

# Make quintiles for each PRS #####

GEND_A1 <- cbind( as.character( GEND_A1), quintiles )
GEND_A2 <- cbind( as.character( GEND_A2), quintiles )
GEND_A3 <- cbind( as.character( GEND_A3), quintiles )
GEND_A4 <- cbind( as.character( GEND_A4), quintiles )
GEND_A5 <- cbind( as.character( GEND_A5), quintiles )
GEND_A6 <- cbind( as.character( GEND_A6), quintiles )
GEND_C1 <- cbind( as.character( GEND_C1), quintiles )
GEND_C2 <- cbind( as.character( GEND_C2), quintiles )
GEND_C3 <- cbind( as.character( GEND_C3), quintiles )
GEND_C4 <- cbind( as.character( GEND_C4), quintiles )
GEND_C5 <- cbind( as.character( GEND_C5), quintiles )
GEND_C6 <- cbind( as.character( GEND_C6), quintiles )
GEND_N1 <- cbind( as.character( GEND_N1), quintiles )
GEND_N2 <- cbind( as.character( GEND_N2), quintiles )
GEND_N3 <- cbind( as.character( GEND_N3), quintiles )
GEND_N4 <- cbind( as.character( GEND_N4), quintiles )
GEND_N5 <- cbind( as.character( GEND_N5), quintiles )
GEND_N6 <- cbind( as.character( GEND_N6), quintiles )
GEND_S1 <- cbind( as.character( GEND_S1), quintiles )
GEND_S2 <- cbind( as.character( GEND_S2), quintiles )
GEND_S3 <- cbind( as.character( GEND_S3), quintiles )
GEND_S4 <- cbind( as.character( GEND_S4), quintiles )
GEND_S5 <- cbind( as.character( GEND_S5), quintiles )
GEND_S6 <- cbind( as.character( GEND_S6), quintiles )
GEND_H1 <- cbind( as.character( GEND_H1), quintiles )
GEND_H2 <- cbind( as.character( GEND_H2), quintiles )
GEND_H3 <- cbind( as.character( GEND_H3), quintiles )
GEND_H4 <- cbind( as.character( GEND_H4), quintiles )
GEND_H5 <- cbind( as.character( GEND_H5), quintiles )
GEND_H6 <- cbind( as.character( GEND_H6), quintiles )

colnames( GEND_A1 ) <- c( "subjectid", "A1_quintile" )
colnames( GEND_A2 ) <- c( "subjectid", "A2_quintile" )
colnames( GEND_A3 ) <- c( "subjectid", "A3_quintile" )
colnames( GEND_A4 ) <- c( "subjectid", "A4_quintile" )
colnames( GEND_A5 ) <- c( "subjectid", "A5_quintile" )
colnames( GEND_A6 ) <- c( "subjectid", "A6_quintile" )

colnames( GEND_C1 ) <- c( "subjectid", "C1_quintile" )
colnames( GEND_C2 ) <- c( "subjectid", "C2_quintile" )
colnames( GEND_C3 ) <- c( "subjectid", "C3_quintile" )
colnames( GEND_C4 ) <- c( "subjectid", "C4_quintile" )
colnames( GEND_C5 ) <- c( "subjectid", "C5_quintile" )
colnames( GEND_C6 ) <- c( "subjectid", "C6_quintile" )

colnames( GEND_N1 ) <- c( "subjectid", "N1_quintile" )
colnames( GEND_N2 ) <- c( "subjectid", "N2_quintile" )
colnames( GEND_N3 ) <- c( "subjectid", "N3_quintile" )
colnames( GEND_N4 ) <- c( "subjectid", "N4_quintile" )
colnames( GEND_N5 ) <- c( "subjectid", "N5_quintile" )
colnames( GEND_N6 ) <- c( "subjectid", "N6_quintile" )

colnames( GEND_S1 ) <- c( "subjectid", "S1_quintile" )
colnames( GEND_S2 ) <- c( "subjectid", "S2_quintile" )
colnames( GEND_S3 ) <- c( "subjectid", "S3_quintile" )
colnames( GEND_S4 ) <- c( "subjectid", "S4_quintile" )
colnames( GEND_S5 ) <- c( "subjectid", "S5_quintile" )
colnames( GEND_S6 ) <- c( "subjectid", "S6_quintile" )

colnames( GEND_H1 ) <- c( "subjectid", "H1_quintile" )
colnames( GEND_H2 ) <- c( "subjectid", "H2_quintile" )
colnames( GEND_H3 ) <- c( "subjectid", "H3_quintile" )
colnames( GEND_H4 ) <- c( "subjectid", "H4_quintile" )
colnames( GEND_H5 ) <- c( "subjectid", "H5_quintile" )
colnames( GEND_H6 ) <- c( "subjectid", "H6_quintile" )


# Attach decile assignments to GENDEP dataframe #####
# 
# GENDEP <- merge( GENDEP, GEND_A1, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_A2, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_A3, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_A4, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_A5, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_A6, by="subjectid" )
# 
# GENDEP <- merge( GENDEP, GEND_C1, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_C2, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_C3, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_C4, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_C5, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_C6, by="subjectid" )
# 
# GENDEP <- merge( GENDEP, GEND_N1, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_N2, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_N3, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_N4, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_N5, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_N6, by="subjectid" )
# 
# GENDEP <- merge( GENDEP, GEND_S1, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_S2, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_S3, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_S4, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_S5, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_S6, by="subjectid" )
# 
# GENDEP <- merge( GENDEP, GEND_H1, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_H2, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_H3, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_H4, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_H5, by="subjectid" )
# GENDEP <- merge( GENDEP, GEND_H6, by="subjectid" )

# Attach quintile assignments to GENDEP dataframe #####

GENDEP <- merge( GENDEP, GEND_A1, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_A2, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_A3, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_A4, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_A5, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_A6, by="subjectid" )

GENDEP <- merge( GENDEP, GEND_C1, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_C2, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_C3, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_C4, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_C5, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_C6, by="subjectid" )

GENDEP <- merge( GENDEP, GEND_N1, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_N2, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_N3, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_N4, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_N5, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_N6, by="subjectid" )

GENDEP <- merge( GENDEP, GEND_S1, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_S2, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_S3, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_S4, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_S5, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_S6, by="subjectid" )

GENDEP <- merge( GENDEP, GEND_H1, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_H2, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_H3, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_H4, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_H5, by="subjectid" )
GENDEP <- merge( GENDEP, GEND_H6, by="subjectid" )

##### Make data frames of dropouts, discontinuers, etc #####

dropouts <- subset(GENDEP,dropout==1)
nondrops <- subset(GENDEP,dropout==0)

drops_end <- subset(dropouts,week==endweek)
nondrops_end <- subset(nondrops,week==endweek)

drop_esc <- subset(dropouts,drug==0)
drop_nor <- subset(dropouts,drug==1)
ndrop_esc <- subset(nondrops,drug==0)
ndrop_nor <- subset(nondrops,drug==1)

drop_esc_end <- subset(drop_esc,week==endweek)
drop_nor_end <- subset(drop_nor,week==endweek)
ndrop_esc_end <- subset(ndrop_esc,week==endweek)
ndrop_nor_end <- subset(ndrop_nor,week==endweek)

drops_end <- drops_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
nondrops_end <- nondrops_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
drop_esc_end <- drop_esc_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
drop_nor_end <- drop_nor_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
ndrop_esc_end <- ndrop_esc_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
ndrop_nor_end <- ndrop_nor_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)

disconts <- subset(GENDEP,dropout2==1)
nondiscs <- subset(GENDEP,dropout2==0)

discs_end <- subset(disconts,week==endweek)
nondiscs_end <- subset(nondiscs,week==endweek)

disc_esc <- subset(disconts,drug==0)
disc_nor <- subset(disconts,drug==1)
ndisc_esc <- subset(nondiscs,drug==0)
ndisc_nor <- subset(nondiscs,drug==1)

disc_esc_end <- subset(disc_esc,week==endweek)
disc_nor_end <- subset(disc_nor,week==endweek)
ndisc_esc_end <- subset(ndisc_esc,week==endweek)
ndisc_nor_end <- subset(ndisc_nor,week==endweek)

discs_end <- discs_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
nondiscs_end <- nondiscs_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
disc_esc_end <- disc_esc_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
disc_nor_end <- disc_nor_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
ndisc_esc_end <- ndisc_esc_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)
ndisc_nor_end <- ndisc_nor_end %>% 
  mutate( percchange = (blmd-madrs)/blmd)

##### Plot baseline depression in different groups #####

ggplot()+
  # geom_density(aes(GEND_end$blmd))+
  # geom_density(aes(GE_esc$blmd))+
  # geom_density(aes(GE_nor$blmd))+
  # geom_density(aes(drops_end$blmd))+
  # geom_density(aes(nondrops_end$blmd))+
  geom_density(aes(drop_nor_end$blmd), col="darkgreen", size=1)+
  geom_density(aes(drop_esc_end$blmd), col="red", size=1)+
  geom_density(aes(ndrop_nor_end$blmd), col="green", size=1)+
  geom_density(aes(ndrop_esc_end$blmd), col="magenta", size=1)+
  xlab("Baseline depression score")

mean(drop_nor_end$blmd)
sd(drop_nor_end$blmd)
mean(drop_esc_end$blmd)
sd(drop_esc_end$blmd)
mean(ndrop_nor_end$blmd)
sd(ndrop_nor_end$blmd)
mean(ndrop_esc_end$blmd)
sd(ndrop_esc_end$blmd)

##### Get means for things - anorexia; both drugs #####
##### Means by deciles #####
# 
# A1_means <- GENDEP %>% 
#   group_by( A1_decile ) %>%
#   summarize( 
#     mad_z = mean( zmadrs, na.rm=TRUE ),
#     ta_z = mean( ztotasec, na.rm=TRUE ),
#     a1 = mean( asec1wk, na.rm=TRUE ),
#     a2 = mean( asec2wk, na.rm=TRUE ),
#     a3 = mean( asec3wk, na.rm=TRUE ),
#     a4 = mean( asec4wk, na.rm=TRUE ),
#     a5 = mean( asec5wk, na.rm=TRUE ),
#     a6 = mean( asec6wk, na.rm=TRUE ),
#     a7 = mean( asec7wk, na.rm=TRUE ),
#     a8 = mean( asec8wk, na.rm=TRUE ),
#     a9 = mean( asec9wk, na.rm=TRUE ),
#     a10 = mean( asec10wk, na.rm=TRUE ),
#     a11 = mean( asec11wk, na.rm=TRUE ),
#     a12 = mean( asec12wk, na.rm=TRUE ),
#     a13 = mean( asec13wk, na.rm=TRUE ),
#     a14 = mean( asec14wk, na.rm=TRUE ),
#     a15 = mean( asec15wk, na.rm=TRUE ),
#     a16 = mean( asec16wk, na.rm=TRUE ),
#     a17 = mean( asec17wk, na.rm=TRUE ),
#     a18 = mean( asec18wk, na.rm=TRUE ),
#     a19 = mean( asec19wk, na.rm=TRUE ),
#     a20 = mean( asec20wk, na.rm=TRUE ),
#     a21 = mean( asec21wk, na.rm=TRUE ),
#     f1_z = mean( zf1score, na.rm=TRUE ),
#     f2_z = mean( zf2score, na.rm=TRUE ),
#     f3_z = mean( zf3score, na.rm=TRUE ),
#     sui_z = mean( suiscore, na.rm=TRUE )
#     )
# 
# A2_means <- GENDEP %>% 
#   group_by( A2_decile ) %>%
#   summarize( 
#     mad_z = mean( zmadrs, na.rm=TRUE ),
#     ta_z = mean( ztotasec, na.rm=TRUE ),
#     a1 = mean( asec1wk, na.rm=TRUE ),
#     a2 = mean( asec2wk, na.rm=TRUE ),
#     a3 = mean( asec3wk, na.rm=TRUE ),
#     a4 = mean( asec4wk, na.rm=TRUE ),
#     a5 = mean( asec5wk, na.rm=TRUE ),
#     a6 = mean( asec6wk, na.rm=TRUE ),
#     a7 = mean( asec7wk, na.rm=TRUE ),
#     a8 = mean( asec8wk, na.rm=TRUE ),
#     a9 = mean( asec9wk, na.rm=TRUE ),
#     a10 = mean( asec10wk, na.rm=TRUE ),
#     a11 = mean( asec11wk, na.rm=TRUE ),
#     a12 = mean( asec12wk, na.rm=TRUE ),
#     a13 = mean( asec13wk, na.rm=TRUE ),
#     a14 = mean( asec14wk, na.rm=TRUE ),
#     a15 = mean( asec15wk, na.rm=TRUE ),
#     a16 = mean( asec16wk, na.rm=TRUE ),
#     a17 = mean( asec17wk, na.rm=TRUE ),
#     a18 = mean( asec18wk, na.rm=TRUE ),
#     a19 = mean( asec19wk, na.rm=TRUE ),
#     a20 = mean( asec20wk, na.rm=TRUE ),
#     a21 = mean( asec21wk, na.rm=TRUE ),
#     f1_z = mean( zf1score, na.rm=TRUE ),
#     f2_z = mean( zf2score, na.rm=TRUE ),
#     f3_z = mean( zf3score, na.rm=TRUE ),
#     sui_z = mean( suiscore, na.rm=TRUE )
#   )
# 
# A3_means <- GENDEP %>% 
#   group_by( A3_decile ) %>%
#   summarize( 
#     mad_z = mean( zmadrs, na.rm=TRUE ),
#     ta_z = mean( ztotasec, na.rm=TRUE ),
#     a1 = mean( asec1wk, na.rm=TRUE ),
#     a2 = mean( asec2wk, na.rm=TRUE ),
#     a3 = mean( asec3wk, na.rm=TRUE ),
#     a4 = mean( asec4wk, na.rm=TRUE ),
#     a5 = mean( asec5wk, na.rm=TRUE ),
#     a6 = mean( asec6wk, na.rm=TRUE ),
#     a7 = mean( asec7wk, na.rm=TRUE ),
#     a8 = mean( asec8wk, na.rm=TRUE ),
#     a9 = mean( asec9wk, na.rm=TRUE ),
#     a10 = mean( asec10wk, na.rm=TRUE ),
#     a11 = mean( asec11wk, na.rm=TRUE ),
#     a12 = mean( asec12wk, na.rm=TRUE ),
#     a13 = mean( asec13wk, na.rm=TRUE ),
#     a14 = mean( asec14wk, na.rm=TRUE ),
#     a15 = mean( asec15wk, na.rm=TRUE ),
#     a16 = mean( asec16wk, na.rm=TRUE ),
#     a17 = mean( asec17wk, na.rm=TRUE ),
#     a18 = mean( asec18wk, na.rm=TRUE ),
#     a19 = mean( asec19wk, na.rm=TRUE ),
#     a20 = mean( asec20wk, na.rm=TRUE ),
#     a21 = mean( asec21wk, na.rm=TRUE ),
#     f1_z = mean( zf1score, na.rm=TRUE ),
#     f2_z = mean( zf2score, na.rm=TRUE ),
#     f3_z = mean( zf3score, na.rm=TRUE ),
#     sui_z = mean( suiscore, na.rm=TRUE )
#   )
# 
# A4_means <- GENDEP %>% 
#   group_by( A4_decile ) %>%
#   summarize( 
#     mad_z = mean( zmadrs, na.rm=TRUE ),
#     ta_z = mean( ztotasec, na.rm=TRUE ),
#     a1 = mean( asec1wk, na.rm=TRUE ),
#     a2 = mean( asec2wk, na.rm=TRUE ),
#     a3 = mean( asec3wk, na.rm=TRUE ),
#     a4 = mean( asec4wk, na.rm=TRUE ),
#     a5 = mean( asec5wk, na.rm=TRUE ),
#     a6 = mean( asec6wk, na.rm=TRUE ),
#     a7 = mean( asec7wk, na.rm=TRUE ),
#     a8 = mean( asec8wk, na.rm=TRUE ),
#     a9 = mean( asec9wk, na.rm=TRUE ),
#     a10 = mean( asec10wk, na.rm=TRUE ),
#     a11 = mean( asec11wk, na.rm=TRUE ),
#     a12 = mean( asec12wk, na.rm=TRUE ),
#     a13 = mean( asec13wk, na.rm=TRUE ),
#     a14 = mean( asec14wk, na.rm=TRUE ),
#     a15 = mean( asec15wk, na.rm=TRUE ),
#     a16 = mean( asec16wk, na.rm=TRUE ),
#     a17 = mean( asec17wk, na.rm=TRUE ),
#     a18 = mean( asec18wk, na.rm=TRUE ),
#     a19 = mean( asec19wk, na.rm=TRUE ),
#     a20 = mean( asec20wk, na.rm=TRUE ),
#     a21 = mean( asec21wk, na.rm=TRUE ),
#     f1_z = mean( zf1score, na.rm=TRUE ),
#     f2_z = mean( zf2score, na.rm=TRUE ),
#     f3_z = mean( zf3score, na.rm=TRUE ),
#     sui_z = mean( suiscore, na.rm=TRUE )
#   )
# 
# A5_means <- GENDEP %>% 
#   group_by( A5_decile ) %>%
#   summarize( 
#     mad_z = mean( zmadrs, na.rm=TRUE ),
#     ta_z = mean( ztotasec, na.rm=TRUE ),
#     a1 = mean( asec1wk, na.rm=TRUE ),
#     a2 = mean( asec2wk, na.rm=TRUE ),
#     a3 = mean( asec3wk, na.rm=TRUE ),
#     a4 = mean( asec4wk, na.rm=TRUE ),
#     a5 = mean( asec5wk, na.rm=TRUE ),
#     a6 = mean( asec6wk, na.rm=TRUE ),
#     a7 = mean( asec7wk, na.rm=TRUE ),
#     a8 = mean( asec8wk, na.rm=TRUE ),
#     a9 = mean( asec9wk, na.rm=TRUE ),
#     a10 = mean( asec10wk, na.rm=TRUE ),
#     a11 = mean( asec11wk, na.rm=TRUE ),
#     a12 = mean( asec12wk, na.rm=TRUE ),
#     a13 = mean( asec13wk, na.rm=TRUE ),
#     a14 = mean( asec14wk, na.rm=TRUE ),
#     a15 = mean( asec15wk, na.rm=TRUE ),
#     a16 = mean( asec16wk, na.rm=TRUE ),
#     a17 = mean( asec17wk, na.rm=TRUE ),
#     a18 = mean( asec18wk, na.rm=TRUE ),
#     a19 = mean( asec19wk, na.rm=TRUE ),
#     a20 = mean( asec20wk, na.rm=TRUE ),
#     a21 = mean( asec21wk, na.rm=TRUE ),
#     f1_z = mean( zf1score, na.rm=TRUE ),
#     f2_z = mean( zf2score, na.rm=TRUE ),
#     f3_z = mean( zf3score, na.rm=TRUE ),
#     sui_z = mean( suiscore, na.rm=TRUE )
#   )
# 
# A6_means <- GENDEP %>% 
#   group_by( A6_decile ) %>%
#   summarize( 
#     mad_z = mean( zmadrs, na.rm=TRUE ),
#     ta_z = mean( ztotasec, na.rm=TRUE ),
#     a1 = mean( asec1wk, na.rm=TRUE ),
#     a2 = mean( asec2wk, na.rm=TRUE ),
#     a3 = mean( asec3wk, na.rm=TRUE ),
#     a4 = mean( asec4wk, na.rm=TRUE ),
#     a5 = mean( asec5wk, na.rm=TRUE ),
#     a6 = mean( asec6wk, na.rm=TRUE ),
#     a7 = mean( asec7wk, na.rm=TRUE ),
#     a8 = mean( asec8wk, na.rm=TRUE ),
#     a9 = mean( asec9wk, na.rm=TRUE ),
#     a10 = mean( asec10wk, na.rm=TRUE ),
#     a11 = mean( asec11wk, na.rm=TRUE ),
#     a12 = mean( asec12wk, na.rm=TRUE ),
#     a13 = mean( asec13wk, na.rm=TRUE ),
#     a14 = mean( asec14wk, na.rm=TRUE ),
#     a15 = mean( asec15wk, na.rm=TRUE ),
#     a16 = mean( asec16wk, na.rm=TRUE ),
#     a17 = mean( asec17wk, na.rm=TRUE ),
#     a18 = mean( asec18wk, na.rm=TRUE ),
#     a19 = mean( asec19wk, na.rm=TRUE ),
#     a20 = mean( asec20wk, na.rm=TRUE ),
#     a21 = mean( asec21wk, na.rm=TRUE ),
#     f1_z = mean( zf1score, na.rm=TRUE ),
#     f2_z = mean( zf2score, na.rm=TRUE ),
#     f3_z = mean( zf3score, na.rm=TRUE ),
#     sui_z = mean( suiscore, na.rm=TRUE )
#   )

##### Means by quintiles #####
# all participants #####

A1_means <- GENDEP %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A2_means <- GENDEP %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A3_means <- GENDEP %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A4_means <- GENDEP %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A5_means <- GENDEP %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A6_means <- GENDEP %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )


# dropouts #####

A1_drop_means <- dropouts %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A2_drop_means <- dropouts %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A3_drop_means <- dropouts %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A4_drop_means <- dropouts %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A5_drop_means <- dropouts %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A6_drop_means <- dropouts %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

# non-dropouts #####

A1_ndrop_means <- nondrops %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A2_ndrop_means <- nondrops %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A3_ndrop_means <- nondrops %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A4_ndrop_means <- nondrops %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A5_ndrop_means <- nondrops %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A6_ndrop_means <- nondrops %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

# non-dropouts - end week #####

A1_ndrop_means_e <- nondrops_end %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A2_ndrop_means_e <- nondrops_end %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A3_ndrop_means_e <- nondrops_end %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A4_ndrop_means_e <- nondrops_end %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A5_ndrop_means_e <- nondrops_end %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A6_ndrop_means_e <- nondrops_end %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

# dropouts - nor #####

A1_drop_nor_means <- drop_nor %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A2_drop_nor_means <- drop_nor %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A3_drop_nor_means <- drop_nor %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A4_drop_nor_means <- drop_nor %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A5_drop_nor_means <- drop_nor %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A6_drop_nor_means <- drop_nor %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

# dropouts - nor - end week #####

A1_drop_nor_means_e <- drop_nor_end %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A2_drop_nor_means_e <- drop_nor_end %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A3_drop_nor_means_e <- drop_nor_end %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A4_drop_nor_means_e <- drop_nor_end %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A5_drop_nor_means_e <- drop_nor_end %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A6_drop_nor_means_e <- drop_nor_end %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )


# dropouts - esc #####

A1_drop_esc_means <- drop_esc %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A2_drop_esc_means <- drop_esc %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A3_drop_esc_means <- drop_esc %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A4_drop_esc_means <- drop_esc %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A5_drop_esc_means <- drop_esc %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A6_drop_esc_means <- drop_esc %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

# dropouts - esc - end week #####

A1_drop_esc_means_e <- drop_esc_end %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A2_drop_esc_means_e <- drop_esc_end %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A3_drop_esc_means_e <- drop_esc_end %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A4_drop_esc_means_e <- drop_esc_end %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A5_drop_esc_means_e <- drop_esc_end %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

A6_drop_esc_means_e <- drop_esc_end %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE ),
    blmadrs_z = mean( zblmd, na.rm=TRUE )
  )

# Add PT info and column names - all #####

A1_means_2 <- cbind( A1_means, rep( "0.0001", length( A1_means$mad_z ) ) )
colnames( A1_means_2 ) <- c( "quintile",
                             colnames( A1_means_2 )[2:28],
                             "PT"
                             )
A2_means_2 <- cbind( A2_means, rep( "0.01", length( A2_means$mad_z ) ) )
colnames( A2_means_2 ) <- c( "quintile",
                             colnames( A2_means_2 )[2:28],
                             "PT"
)
A3_means_2 <- cbind( A3_means, rep( "0.05", length( A3_means$mad_z ) ) )
colnames( A3_means_2 ) <- c( "quintile",
                             colnames( A3_means_2 )[2:28],
                             "PT"
)
A4_means_2 <- cbind( A4_means, rep( "0.1", length( A4_means$mad_z ) ) )
colnames( A4_means_2 ) <- c( "quintile",
                             colnames( A4_means_2 )[2:28],
                             "PT"
)
A5_means_2 <- cbind( A5_means, rep( "0.5", length( A5_means$mad_z ) ) )
colnames( A5_means_2 ) <- c( "quintile",
                             colnames( A5_means_2 )[2:28],
                             "PT"
)
A6_means_2 <- cbind( A6_means, rep( "1.0", length( A6_means$mad_z ) ) )
colnames( A6_means_2 ) <- c( "quintile",
                             colnames( A6_means_2 )[2:28],
                             "PT"
)


A_means <- rbind(
  A1_means_2, A2_means_2, A3_means_2,
  A4_means_2, A5_means_2, A6_means_2
)

# Add PT info and column names - dropouts #####

A1_drop_means_2 <- cbind( A1_drop_means, rep( "0.0001", length( A1_drop_means$mad_z ) ) )
colnames( A1_drop_means_2 ) <- c( "quintile",
                             colnames( A1_drop_means_2 )[2:28],
                             "PT"
)
A2_drop_means_2 <- cbind( A2_drop_means, rep( "0.01", length( A2_drop_means$mad_z ) ) )
colnames( A2_drop_means_2 ) <- c( "quintile",
                             colnames( A2_drop_means_2 )[2:28],
                             "PT"
)
A3_drop_means_2 <- cbind( A3_drop_means, rep( "0.05", length( A3_drop_means$mad_z ) ) )
colnames( A3_drop_means_2 ) <- c( "quintile",
                             colnames( A3_drop_means_2 )[2:28],
                             "PT"
)
A4_drop_means_2 <- cbind( A4_drop_means, rep( "0.1", length( A4_drop_means$mad_z ) ) )
colnames( A4_drop_means_2 ) <- c( "quintile",
                             colnames( A4_drop_means_2 )[2:28],
                             "PT"
)
A5_drop_means_2 <- cbind( A5_drop_means, rep( "0.5", length( A5_drop_means$mad_z ) ) )
colnames( A5_drop_means_2 ) <- c( "quintile",
                             colnames( A5_drop_means_2 )[2:28],
                             "PT"
)
A6_drop_means_2 <- cbind( A6_drop_means, rep( "1.0", length( A6_drop_means$mad_z ) ) )
colnames( A6_drop_means_2 ) <- c( "quintile",
                             colnames( A6_drop_means_2 )[2:28],
                             "PT"
)


A_drop_means <- rbind(
  A1_drop_means_2, A2_drop_means_2, A3_drop_means_2,
  A4_drop_means_2, A5_drop_means_2, A6_drop_means_2
)

# Add PT info and column names - non-dropouts #####

A1_ndrop_means_2 <- cbind( A1_ndrop_means, rep( "0.0001", length( A1_ndrop_means$mad_z ) ) )
colnames( A1_ndrop_means_2 ) <- c( "quintile",
                             colnames( A1_ndrop_means_2 )[2:28],
                             "PT"
)
A2_ndrop_means_2 <- cbind( A2_ndrop_means, rep( "0.01", length( A2_ndrop_means$mad_z ) ) )
colnames( A2_ndrop_means_2 ) <- c( "quintile",
                             colnames( A2_ndrop_means_2 )[2:28],
                             "PT"
)
A3_ndrop_means_2 <- cbind( A3_ndrop_means, rep( "0.05", length( A3_ndrop_means$mad_z ) ) )
colnames( A3_ndrop_means_2 ) <- c( "quintile",
                             colnames( A3_ndrop_means_2 )[2:28],
                             "PT"
)
A4_ndrop_means_2 <- cbind( A4_ndrop_means, rep( "0.1", length( A4_ndrop_means$mad_z ) ) )
colnames( A4_ndrop_means_2 ) <- c( "quintile",
                             colnames( A4_ndrop_means_2 )[2:28],
                             "PT"
)
A5_ndrop_means_2 <- cbind( A5_ndrop_means, rep( "0.5", length( A5_ndrop_means$mad_z ) ) )
colnames( A5_ndrop_means_2 ) <- c( "quintile",
                             colnames( A5_ndrop_means_2 )[2:28],
                             "PT"
)
A6_ndrop_means_2 <- cbind( A6_ndrop_means, rep( "1.0", length( A6_ndrop_means$mad_z ) ) )
colnames( A6_ndrop_means_2 ) <- c( "quintile",
                             colnames( A6_ndrop_means_2 )[2:28],
                             "PT"
)


A_ndrop_means <- rbind(
  A1_ndrop_means_2, A2_ndrop_means_2, A3_ndrop_means_2,
  A4_ndrop_means_2, A5_ndrop_means_2, A6_ndrop_means_2
)

A1_ndrop_means_2_e <- cbind( A1_ndrop_means_e, rep( "0.0001", length( A1_ndrop_means_e$mad_z ) ) )
colnames( A1_ndrop_means_2_e ) <- c( "quintile",
                                   colnames( A1_ndrop_means_2_e )[2:29],
                                   "PT"
)
A2_ndrop_means_2_e <- cbind( A2_ndrop_means_e, rep( "0.01", length( A2_ndrop_means_e$mad_z ) ) )
colnames( A2_ndrop_means_2_e ) <- c( "quintile",
                                   colnames( A2_ndrop_means_2_e )[2:29],
                                   "PT"
)
A3_ndrop_means_2_e <- cbind( A3_ndrop_means_e, rep( "0.05", length( A3_ndrop_means_e$mad_z ) ) )
colnames( A3_ndrop_means_2_e ) <- c( "quintile",
                                   colnames( A3_ndrop_means_2_e )[2:29],
                                   "PT"
)
A4_ndrop_means_2_e <- cbind( A4_ndrop_means_e, rep( "0.1", length( A4_ndrop_means_e$mad_z ) ) )
colnames( A4_ndrop_means_2_e ) <- c( "quintile",
                                   colnames( A4_ndrop_means_2_e )[2:29],
                                   "PT"
)
A5_ndrop_means_2_e <- cbind( A5_ndrop_means_e, rep( "0.5", length( A5_ndrop_means_e$mad_z ) ) )
colnames( A5_ndrop_means_2_e ) <- c( "quintile",
                                   colnames( A5_ndrop_means_2_e )[2:29],
                                   "PT"
)
A6_ndrop_means_2_e <- cbind( A6_ndrop_means_e, rep( "1.0", length( A6_ndrop_means_e$mad_z ) ) )
colnames( A6_ndrop_means_2_e ) <- c( "quintile",
                                   colnames( A6_ndrop_means_2_e )[2:29],
                                   "PT"
)


A_ndrop_means_e <- rbind(
  A1_ndrop_means_2_e, A2_ndrop_means_2_e, A3_ndrop_means_2_e,
  A4_ndrop_means_2_e, A5_ndrop_means_2_e, A6_ndrop_means_2_e
)

# Add PT info and column names - nor dropouts #####

A1_drop_nor_means_2 <- cbind( A1_drop_nor_means, rep( "0.0001", length( A1_drop_nor_means$mad_z ) ) )
colnames( A1_drop_nor_means_2 ) <- c( "quintile",
                                  colnames( A1_drop_nor_means_2 )[2:28],
                                  "PT"
)
A2_drop_nor_means_2 <- cbind( A2_drop_nor_means, rep( "0.01", length( A2_drop_nor_means$mad_z ) ) )
colnames( A2_drop_nor_means_2 ) <- c( "quintile",
                                  colnames( A2_drop_nor_means_2 )[2:28],
                                  "PT"
)
A3_drop_nor_means_2 <- cbind( A3_drop_nor_means, rep( "0.05", length( A3_drop_nor_means$mad_z ) ) )
colnames( A3_drop_nor_means_2 ) <- c( "quintile",
                                  colnames( A3_drop_nor_means_2 )[2:28],
                                  "PT"
)
A4_drop_nor_means_2 <- cbind( A4_drop_nor_means, rep( "0.1", length( A4_drop_nor_means$mad_z ) ) )
colnames( A4_drop_nor_means_2 ) <- c( "quintile",
                                  colnames( A4_drop_nor_means_2 )[2:28],
                                  "PT"
)
A5_drop_nor_means_2 <- cbind( A5_drop_nor_means, rep( "0.5", length( A5_drop_nor_means$mad_z ) ) )
colnames( A5_drop_nor_means_2 ) <- c( "quintile",
                                  colnames( A5_drop_nor_means_2 )[2:28],
                                  "PT"
)
A6_drop_nor_means_2 <- cbind( A6_drop_nor_means, rep( "1.0", length( A6_drop_nor_means$mad_z ) ) )
colnames( A6_drop_nor_means_2 ) <- c( "quintile",
                                  colnames( A6_drop_nor_means_2 )[2:28],
                                  "PT"
)


A_drop_nor_means <- rbind(
  A1_drop_nor_means_2, A2_drop_nor_means_2, A3_drop_nor_means_2,
  A4_drop_nor_means_2, A5_drop_nor_means_2, A6_drop_nor_means_2
)


A1_drop_nor_means_2_e <- cbind( A1_drop_nor_means_e, rep( "0.0001", length( A1_drop_nor_means_e$mad_z ) ) )
colnames( A1_drop_nor_means_2_e ) <- c( "quintile",
                                      colnames( A1_drop_nor_means_2_e )[2:29],
                                      "PT"
)
A2_drop_nor_means_2_e <- cbind( A2_drop_nor_means_e, rep( "0.01", length( A2_drop_nor_means_e$mad_z ) ) )
colnames( A2_drop_nor_means_2_e ) <- c( "quintile",
                                      colnames( A2_drop_nor_means_2_e )[2:29],
                                      "PT"
)
A3_drop_nor_means_2_e <- cbind( A3_drop_nor_means_e, rep( "0.05", length( A3_drop_nor_means_e$mad_z ) ) )
colnames( A3_drop_nor_means_2_e ) <- c( "quintile",
                                      colnames( A3_drop_nor_means_2_e )[2:29],
                                      "PT"
)
A4_drop_nor_means_2_e <- cbind( A4_drop_nor_means_e, rep( "0.1", length( A4_drop_nor_means_e$mad_z ) ) )
colnames( A4_drop_nor_means_2_e ) <- c( "quintile",
                                      colnames( A4_drop_nor_means_2_e )[2:29],
                                      "PT"
)
A5_drop_nor_means_2_e <- cbind( A5_drop_nor_means_e, rep( "0.5", length( A5_drop_nor_means_e$mad_z ) ) )
colnames( A5_drop_nor_means_2_e ) <- c( "quintile",
                                      colnames( A5_drop_nor_means_2_e )[2:29],
                                      "PT"
)
A6_drop_nor_means_2_e <- cbind( A6_drop_nor_means_e, rep( "1.0", length( A6_drop_nor_means_e$mad_z ) ) )
colnames( A6_drop_nor_means_2_e ) <- c( "quintile",
                                      colnames( A6_drop_nor_means_2_e )[2:29],
                                      "PT"
)


A_drop_nor_means_e <- rbind(
  A1_drop_nor_means_2_e, A2_drop_nor_means_2_e, A3_drop_nor_means_2_e,
  A4_drop_nor_means_2_e, A5_drop_nor_means_2_e, A6_drop_nor_means_2_e
)

# Add PT info and column names - esc dropouts #####

A1_drop_esc_means_2 <- cbind( A1_drop_esc_means, rep( "0.0001", length( A1_drop_esc_means$mad_z ) ) )
colnames( A1_drop_esc_means_2 ) <- c( "quintile",
                                  colnames( A1_drop_esc_means_2 )[2:28],
                                  "PT"
)
A2_drop_esc_means_2 <- cbind( A2_drop_esc_means, rep( "0.01", length( A2_drop_esc_means$mad_z ) ) )
colnames( A2_drop_esc_means_2 ) <- c( "quintile",
                                  colnames( A2_drop_esc_means_2 )[2:28],
                                  "PT"
)
A3_drop_esc_means_2 <- cbind( A3_drop_esc_means, rep( "0.05", length( A3_drop_esc_means$mad_z ) ) )
colnames( A3_drop_esc_means_2 ) <- c( "quintile",
                                  colnames( A3_drop_esc_means_2 )[2:28],
                                  "PT"
)
A4_drop_esc_means_2 <- cbind( A4_drop_esc_means, rep( "0.1", length( A4_drop_esc_means$mad_z ) ) )
colnames( A4_drop_esc_means_2 ) <- c( "quintile",
                                  colnames( A4_drop_esc_means_2 )[2:28],
                                  "PT"
)
A5_drop_esc_means_2 <- cbind( A5_drop_esc_means, rep( "0.5", length( A5_drop_esc_means$mad_z ) ) )
colnames( A5_drop_esc_means_2 ) <- c( "quintile",
                                  colnames( A5_drop_esc_means_2 )[2:28],
                                  "PT"
)
A6_drop_esc_means_2 <- cbind( A6_drop_esc_means, rep( "1.0", length( A6_drop_esc_means$mad_z ) ) )
colnames( A6_drop_esc_means_2 ) <- c( "quintile",
                                  colnames( A6_drop_esc_means_2 )[2:28],
                                  "PT"
)


A_drop_esc_means <- rbind(
  A1_drop_esc_means_2, A2_drop_esc_means_2, A3_drop_esc_means_2,
  A4_drop_esc_means_2, A5_drop_esc_means_2, A6_drop_esc_means_2
)


A1_drop_esc_means_2_e <- cbind( A1_drop_esc_means_e, rep( "0.0001", length( A1_drop_esc_means_e$mad_z ) ) )
colnames( A1_drop_esc_means_2_e ) <- c( "quintile",
                                      colnames( A1_drop_esc_means_2_e )[2:29],
                                      "PT"
)
A2_drop_esc_means_2_e <- cbind( A2_drop_esc_means_e, rep( "0.01", length( A2_drop_esc_means_e$mad_z ) ) )
colnames( A2_drop_esc_means_2_e ) <- c( "quintile",
                                      colnames( A2_drop_esc_means_2_e )[2:29],
                                      "PT"
)
A3_drop_esc_means_2_e <- cbind( A3_drop_esc_means_e, rep( "0.05", length( A3_drop_esc_means_e$mad_z ) ) )
colnames( A3_drop_esc_means_2_e ) <- c( "quintile",
                                      colnames( A3_drop_esc_means_2_e )[2:29],
                                      "PT"
)
A4_drop_esc_means_2_e <- cbind( A4_drop_esc_means_e, rep( "0.1", length( A4_drop_esc_means_e$mad_z ) ) )
colnames( A4_drop_esc_means_2_e ) <- c( "quintile",
                                      colnames( A4_drop_esc_means_2_e )[2:29],
                                      "PT"
)
A5_drop_esc_means_2_e <- cbind( A5_drop_esc_means_e, rep( "0.5", length( A5_drop_esc_means_e$mad_z ) ) )
colnames( A5_drop_esc_means_2_e ) <- c( "quintile",
                                      colnames( A5_drop_esc_means_2_e )[2:29],
                                      "PT"
)
A6_drop_esc_means_2_e <- cbind( A6_drop_esc_means_e, rep( "1.0", length( A6_drop_esc_means_e$mad_z ) ) )
colnames( A6_drop_esc_means_2_e ) <- c( "quintile",
                                      colnames( A6_drop_esc_means_2_e )[2:29],
                                      "PT"
)


A_drop_esc_means_e <- rbind(
  A1_drop_esc_means_2_e, A2_drop_esc_means_2_e, A3_drop_esc_means_2_e,
  A4_drop_esc_means_2_e, A5_drop_esc_means_2_e, A6_drop_esc_means_2_e
)

# Plot mean MADRS for each obs in different groups by anorexia quintile #####

ggplot(A_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
    )+
  labs(title="Anorexia PS, MADRS")

ggplot(A_drop_nor_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, MADRS, dropouts, nortriptyline")

ggplot(A_drop_esc_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, MADRS, dropouts, escitalopram")

ggplot(A_ndrop_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, MADRS, non-dropouts")

ggplot(A_ndrop_means_e,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean end-week z-MADRS, non-dropouts, by anorexia nervosa PS quintile")+
  ylim(c(-0.7,1.2))

ggplot(A_drop_esc_means_e,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean end-week z-MADRS, escitalopram dropouts, by anorexia nervosa PS quintile")+
  ylim(c(-0.7,1.2))

ggplot(A_drop_nor_means_e,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean end-week z-MADRS, nortriptyline dropouts, by anorexia nervosa PS quintile")+
  ylim(c(-0.7,1.2))

# Plot mean total ASEC, and other outcomes, for each obs in "all" by anorexia quintile #####

ggplot(A_means,aes(x=quintile,y=ta_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, Totasec")

ggplot(A_means,aes(x=quintile,y=f1_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, F1 score (mood)")

ggplot(A_means,aes(x=quintile,y=f2_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, F2 score (cognitive)")

ggplot(A_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, F3 score (neurovegetative)")

ggplot(A_means,aes(x=quintile,y=sui_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, suicidality score")

ggplot(A_means,aes(x=quintile,y=a1,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, dry mouth")

ggplot(A_means,aes(x=quintile,y=a2,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, drowsiness")

ggplot(A_means,aes(x=quintile,y=a3,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, insomnia")

ggplot(A_means,aes(x=quintile,y=a4,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, blurred vision")

ggplot(A_means,aes(x=quintile,y=a5,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, headache")

ggplot(A_means,aes(x=quintile,y=a6,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, constipation")

ggplot(A_means,aes(x=quintile,y=a7,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, diarrhoea")

ggplot(A_means,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, increased appetite")

ggplot(A_ndrop_means_e,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean end week increased appetite, non-dropouts, by anorexia nervosa PS quintile")+
  ylim(c(0,0.5))

ggplot(A_drop_esc_means_e,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean end week increased appetite, escitalopram dropouts, by anorexia nervosa PS quintile")+
  ylim(c(0,0.5))

ggplot(A_drop_nor_means_e,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean end week increased appetite, nortriptyline dropouts, by anorexia nervosa PS quintile")+
  ylim(c(0,0.5))

# Plot baseline MADRS for non-dropouts, esc dropouts and nor dropouts

ggplot(A_ndrop_means_e,aes(x=quintile,y=blmadrs_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean baseline depression z-score, non-dropouts, by anorexia nervosa PS quintile")+
  ylim(c(-1,1))

ggplot(A_drop_esc_means_e,aes(x=quintile,y=blmadrs_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Mean baseline depression z-score, escitalopram dropouts, by anorexia nervosa PS quintile")+
  ylim(c(-1,1))

ggplot(A_drop_nor_means_e,aes(x=quintile,y=blmadrs_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7",
    "#00ACC1", "#00BCD4", "#26C6DA")
  )+
  labs(title="Mean baseline depression z-score, nortriptyline dropouts, by anorexia nervosa PS quintile")+
  ylim(c(-1,1))

ggplot(data=nondrops_end,aes(zblmd))+
  geom_density(fill="blue",alpha=0.5)+
  geom_density(data=drop_esc_end,fill="magenta",alpha=0.5)+
  geom_density(data=drop_nor_end,fill="green",alpha=0.5)

# Resume plotting side effects / ADRs

ggplot(A_means,aes(x=quintile,y=a9,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, decreased appetite")

ggplot(A_means,aes(x=quintile,y=a10,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nausea / vomiting")

ggplot(A_means,aes(x=quintile,y=a11,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, urination problems")

ggplot(A_means,aes(x=quintile,y=a12,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, sexual dysfunction")

ggplot(A_means,aes(x=quintile,y=a13,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PRS, palpitations")

ggplot(A_means,aes(x=quintile,y=a14,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, light-headedness")

ggplot(A_means,aes(x=quintile,y=a15,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, room-spinning")

ggplot(A_means,aes(x=quintile,y=a16,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, sweating")

ggplot(A_means,aes(x=quintile,y=a17,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, increased body temperature")

ggplot(A_means,aes(x=quintile,y=a18,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, tremor")

ggplot(A_means,aes(x=quintile,y=a19,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, disorientation")

ggplot(A_means,aes(x=quintile,y=a20,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, yawning")

ggplot(A_means,aes(x=quintile,y=a21,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, weight gain")

# Other anorexia PS quintile related plots #####

ggplot( GENDEP, aes(x=asec8wk,y=f3score,col=drug))+
  geom_jitter(size=1)+
  facet_wrap( ~A4_quintile )+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()+
  xlab("Increased appetite rating (0-3)")+
  ylab("Neurovegetative score, z-scored")+
  ggtitle("Neurovegetative score by increased appetite rating, drug and anorexia PS quintile (PT 0.1), jittered points")

ggplot( GENDEP, aes(x=AZ_1,y=f3score,col=drug))+
  geom_point(size=1)+
  facet_wrap( ~asec8wk )+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()+
  xlab("Anorexia PS (PT 1.0), z-scored")+
  ylab("Neurovegetative score, z-scored")+
  ggtitle("Neurovegetative score by anorexia PS (PT 1.0), drug and increased appetite rating")

##### Get means for things - anorexia; escitalopram #####

# Make GENDEP sub-tables for each drug

GENDEP_escit <- subset( GENDEP, drug==0 ) 

A1_e_means <- GENDEP_escit %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A2_e_means <- GENDEP_escit %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A3_e_means <- GENDEP_escit %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A4_e_means <- GENDEP_escit %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A5_e_means <- GENDEP_escit %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A6_e_means <- GENDEP_escit %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A1_e_means_2 <- cbind( A1_e_means, rep( "0.0001", length( A1_e_means$mad_z ) ) )
colnames( A1_e_means_2 ) <- c( "quintile",
                             colnames( A1_e_means_2 )[2:28],
                             "PT"
)
A2_e_means_2 <- cbind( A2_e_means, rep( "0.01", length( A2_e_means$mad_z ) ) )
colnames( A2_e_means_2 ) <- c( "quintile",
                             colnames( A2_e_means_2 )[2:28],
                             "PT"
)
A3_e_means_2 <- cbind( A3_e_means, rep( "0.05", length( A3_e_means$mad_z ) ) )
colnames( A3_e_means_2 ) <- c( "quintile",
                             colnames( A3_e_means_2 )[2:28],
                             "PT"
)
A4_e_means_2 <- cbind( A4_e_means, rep( "0.1", length( A4_e_means$mad_z ) ) )
colnames( A4_e_means_2 ) <- c( "quintile",
                             colnames( A4_e_means_2 )[2:28],
                             "PT"
)
A5_e_means_2 <- cbind( A5_e_means, rep( "0.5", length( A5_e_means$mad_z ) ) )
colnames( A5_e_means_2 ) <- c( "quintile",
                             colnames( A5_e_means_2 )[2:28],
                             "PT"
)
A6_e_means_2 <- cbind( A6_e_means, rep( "1.0", length( A6_e_means$mad_z ) ) )
colnames( A6_e_means_2 ) <- c( "quintile",
                             colnames( A6_e_means_2 )[2:28],
                             "PT"
)


A_e_means <- rbind(
  A1_e_means_2, A2_e_means_2, A3_e_means_2,
  A4_e_means_2, A5_e_means_2, A6_e_means_2
)


# Plot mean MADRS for each obs in escitalopram-only by anorexia quintile #####

ggplot(A_e_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, MADRS")

# Plot mean total ASEC, and other outcomes, for each obs in escitalopram-only by anorexia quintile #####

ggplot(A_e_means,aes(x=quintile,y=ta_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, Totasec")

ggplot(A_e_means,aes(x=quintile,y=f1_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, F1 score (mood)")

ggplot(A_e_means,aes(x=quintile,y=f2_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, F2 score (cognitive)")

ggplot(A_e_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, F3 score (neurovegetative)")

ggplot(A_e_means,aes(x=quintile,y=sui_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, suicidality score")

ggplot(A_e_means,aes(x=quintile,y=a1,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, dry mouth")

ggplot(A_e_means,aes(x=quintile,y=a2,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, drowsiness")

ggplot(A_e_means,aes(x=quintile,y=a3,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, insomnia")

ggplot(A_e_means,aes(x=quintile,y=a4,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, blurred vision")

ggplot(A_e_means,aes(x=quintile,y=a5,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, headache")

ggplot(A_e_means,aes(x=quintile,y=a6,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, constipation")

ggplot(A_e_means,aes(x=quintile,y=a7,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, diarrhoea")

ggplot(A_e_means,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, increased appetite")

ggplot(A_e_means,aes(x=quintile,y=a9,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, decreased appetite")

ggplot(A_e_means,aes(x=quintile,y=a10,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, nausea / vomiting")

ggplot(A_e_means,aes(x=quintile,y=a11,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, urination problems")

ggplot(A_e_means,aes(x=quintile,y=a12,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, sexual dysfunction")

ggplot(A_e_means,aes(x=quintile,y=a13,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, palpitations")

ggplot(A_e_means,aes(x=quintile,y=a14,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, light-headedness")

ggplot(A_e_means,aes(x=quintile,y=a15,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, room-spinning")

ggplot(A_e_means,aes(x=quintile,y=a16,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, sweating")

ggplot(A_e_means,aes(x=quintile,y=a17,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, increased body temperature")

ggplot(A_e_means,aes(x=quintile,y=a18,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, tremor")

ggplot(A_e_means,aes(x=quintile,y=a19,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, disorientation")

ggplot(A_e_means,aes(x=quintile,y=a20,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, yawning")

ggplot(A_e_means,aes(x=quintile,y=a21,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, escitalopram-only, weight gain")


##### Get means for things - anorexia; nortriptyline #####

GENDEP_nortrip <- subset( GENDEP, drug==1 )

A1_n_means <- GENDEP_nortrip %>% 
  group_by( A1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A2_n_means <- GENDEP_nortrip %>% 
  group_by( A2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A3_n_means <- GENDEP_nortrip %>% 
  group_by( A3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A4_n_means <- GENDEP_nortrip %>% 
  group_by( A4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A5_n_means <- GENDEP_nortrip %>% 
  group_by( A5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A6_n_means <- GENDEP_nortrip %>% 
  group_by( A6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

A1_n_means_2 <- cbind( A1_n_means, rep( "0.0001", length( A1_n_means$mad_z ) ) )
colnames( A1_n_means_2 ) <- c( "quintile",
                               colnames( A1_n_means_2 )[2:28],
                               "PT"
)
A2_n_means_2 <- cbind( A2_n_means, rep( "0.01", length( A2_n_means$mad_z ) ) )
colnames( A2_n_means_2 ) <- c( "quintile",
                               colnames( A2_n_means_2 )[2:28],
                               "PT"
)
A3_n_means_2 <- cbind( A3_n_means, rep( "0.05", length( A3_n_means$mad_z ) ) )
colnames( A3_n_means_2 ) <- c( "quintile",
                               colnames( A3_n_means_2 )[2:28],
                               "PT"
)
A4_n_means_2 <- cbind( A4_n_means, rep( "0.1", length( A4_n_means$mad_z ) ) )
colnames( A4_n_means_2 ) <- c( "quintile",
                               colnames( A4_n_means_2 )[2:28],
                               "PT"
)
A5_n_means_2 <- cbind( A5_n_means, rep( "0.5", length( A5_n_means$mad_z ) ) )
colnames( A5_n_means_2 ) <- c( "quintile",
                               colnames( A5_n_means_2 )[2:28],
                               "PT"
)
A6_n_means_2 <- cbind( A6_n_means, rep( "1.0", length( A6_n_means$mad_z ) ) )
colnames( A6_n_means_2 ) <- c( "quintile",
                               colnames( A6_n_means_2 )[2:28],
                               "PT"
)


A_n_means <- rbind(
  A1_n_means_2, A2_n_means_2, A3_n_means_2,
  A4_n_means_2, A5_n_means_2, A6_n_means_2
)

# Plot mean MADRS for each obs in nortriptyline-only by anorexia quintile #####

ggplot(A_n_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, MADRS")

# Plot mean total ASEC, and other outcomes, for each obs in nortriptyline-only by anorexia quintile #####

ggplot(A_n_means,aes(x=quintile,y=ta_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, Totasec")

ggplot(A_n_means,aes(x=quintile,y=f1_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, F1 score (mood)")

ggplot(A_n_means,aes(x=quintile,y=f2_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, F2 score (cognitive)")

ggplot(A_n_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, F3 score (neurovegetative)")

ggplot(A_n_means,aes(x=quintile,y=sui_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, suicidality score")

ggplot(A_n_means,aes(x=quintile,y=a1,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, dry mouth")

ggplot(A_n_means,aes(x=quintile,y=a2,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, drowsiness")

ggplot(A_n_means,aes(x=quintile,y=a3,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, insomnia")

ggplot(A_n_means,aes(x=quintile,y=a4,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, blurred vision")

ggplot(A_n_means,aes(x=quintile,y=a5,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, headache")

ggplot(A_n_means,aes(x=quintile,y=a6,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, constipation")

ggplot(A_n_means,aes(x=quintile,y=a7,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, diarrhoea")

ggplot(A_n_means,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, increased appetite")

ggplot(A_n_means,aes(x=quintile,y=a9,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, decreased appetite")

ggplot(A_n_means,aes(x=quintile,y=a10,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, nausea / vomiting")

ggplot(A_n_means,aes(x=quintile,y=a11,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, urination problems")

ggplot(A_n_means,aes(x=quintile,y=a12,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, sexual dysfunction")

ggplot(A_n_means,aes(x=quintile,y=a13,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, palpitations")

ggplot(A_n_means,aes(x=quintile,y=a14,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, light-headedness")

ggplot(A_n_means,aes(x=quintile,y=a15,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, room-spinning")

ggplot(A_n_means,aes(x=quintile,y=a16,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, sweating")

ggplot(A_n_means,aes(x=quintile,y=a17,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, increased body temperature")

ggplot(A_n_means,aes(x=quintile,y=a18,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, tremor")

ggplot(A_n_means,aes(x=quintile,y=a19,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, disorientation")

ggplot(A_n_means,aes(x=quintile,y=a20,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, yawning")

ggplot(A_n_means,aes(x=quintile,y=a21,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, nortriptyline-only, weight gain")

ggplot(A_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Anorexia PS, F3 score (neurovegetative)")

##### Get means for things - reaction time; both drugs #####

C1_means <- GENDEP %>% 
  group_by( C1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C2_means <- GENDEP %>% 
  group_by( C2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C3_means <- GENDEP %>% 
  group_by( C3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C4_means <- GENDEP %>% 
  group_by( C4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C5_means <- GENDEP %>% 
  group_by( C5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C6_means <- GENDEP %>% 
  group_by( C6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C1_means_2 <- cbind( C1_means, rep( "0.0001", length( C1_means$mad_z ) ) )
colnames( C1_means_2 ) <- c( "quintile",
                             colnames( C1_means_2 )[2:28],
                             "PT"
)
C2_means_2 <- cbind( C2_means, rep( "0.01", length( C2_means$mad_z ) ) )
colnames( C2_means_2 ) <- c( "quintile",
                             colnames( C2_means_2 )[2:28],
                             "PT"
)
C3_means_2 <- cbind( C3_means, rep( "0.05", length( C3_means$mad_z ) ) )
colnames( C3_means_2 ) <- c( "quintile",
                             colnames( C3_means_2 )[2:28],
                             "PT"
)
C4_means_2 <- cbind( C4_means, rep( "0.1", length( C4_means$mad_z ) ) )
colnames( C4_means_2 ) <- c( "quintile",
                             colnames( C4_means_2 )[2:28],
                             "PT"
)
C5_means_2 <- cbind( C5_means, rep( "0.5", length( C5_means$mad_z ) ) )
colnames( C5_means_2 ) <- c( "quintile",
                             colnames( C5_means_2 )[2:28],
                             "PT"
)
C6_means_2 <- cbind( C6_means, rep( "1.0", length( C6_means$mad_z ) ) )
colnames( C6_means_2 ) <- c( "quintile",
                             colnames( C6_means_2 )[2:28],
                             "PT"
)


C_means <- rbind(
  C1_means_2, C2_means_2, C3_means_2,
  C4_means_2, C5_means_2, C6_means_2
)

# Plot all outcomes for "all" by reaction time quintile #####

ggplot(C_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, MADRS")

ggplot(C_means,aes(x=quintile,y=ta_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, Totasec")

ggplot(C_means,aes(x=quintile,y=f1_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, F1 score (mood)")

ggplot(C_means,aes(x=quintile,y=f2_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, F2 score (cognitive)")

ggplot(C_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, F3 score (neurovegetative)")

ggplot(C_means,aes(x=quintile,y=sui_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, suicidality score")

ggplot(C_means,aes(x=quintile,y=a1,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, dry mouth")

ggplot(C_means,aes(x=quintile,y=a2,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, drowsiness")

ggplot(C_means,aes(x=quintile,y=a3,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, insomnia")

ggplot(C_means,aes(x=quintile,y=a4,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, blurred vision")

ggplot(C_means,aes(x=quintile,y=a5,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, headache")

ggplot(C_means,aes(x=quintile,y=a6,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, constipation")

ggplot(C_means,aes(x=quintile,y=a7,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, diarrhoea")

ggplot(C_means,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, increased appetite")

ggplot(C_means,aes(x=quintile,y=a9,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, decreased appetite")

ggplot(C_means,aes(x=quintile,y=a10,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nausea / vomiting")

ggplot(C_means,aes(x=quintile,y=a11,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, urination problems")

ggplot(C_means,aes(x=quintile,y=a12,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, sexual dysfunction")

ggplot(C_means,aes(x=quintile,y=a13,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, palpitations")

ggplot(C_means,aes(x=quintile,y=a14,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, light-headedness")

ggplot(C_means,aes(x=quintile,y=a15,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, room-spinning")

ggplot(C_means,aes(x=quintile,y=a16,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, sweating")

ggplot(C_means,aes(x=quintile,y=a17,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, increased body temperature")

ggplot(C_means,aes(x=quintile,y=a18,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, tremor")

ggplot(C_means,aes(x=quintile,y=a19,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, disorientation")

ggplot(C_means,aes(x=quintile,y=a20,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, yawning")

ggplot(C_means,aes(x=quintile,y=a21,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, weight gain")

##### Get means for things - Reaction Time; escitalopram #####

# Make GENDEP sub-tables for each drug

C1_e_means <- GENDEP_escit %>% 
  group_by( C1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C2_e_means <- GENDEP_escit %>% 
  group_by( C2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C3_e_means <- GENDEP_escit %>% 
  group_by( C3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C4_e_means <- GENDEP_escit %>% 
  group_by( C4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C5_e_means <- GENDEP_escit %>% 
  group_by( C5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C6_e_means <- GENDEP_escit %>% 
  group_by( C6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C1_e_means_2 <- cbind( C1_e_means, rep( "0.0001", length( C1_e_means$mad_z ) ) )
colnames( C1_e_means_2 ) <- c( "quintile",
                               colnames( C1_e_means_2 )[2:28],
                               "PT"
)
C2_e_means_2 <- cbind( C2_e_means, rep( "0.01", length( C2_e_means$mad_z ) ) )
colnames( C2_e_means_2 ) <- c( "quintile",
                               colnames( C2_e_means_2 )[2:28],
                               "PT"
)
C3_e_means_2 <- cbind( C3_e_means, rep( "0.05", length( C3_e_means$mad_z ) ) )
colnames( C3_e_means_2 ) <- c( "quintile",
                               colnames( C3_e_means_2 )[2:28],
                               "PT"
)
C4_e_means_2 <- cbind( C4_e_means, rep( "0.1", length( C4_e_means$mad_z ) ) )
colnames( C4_e_means_2 ) <- c( "quintile",
                               colnames( C4_e_means_2 )[2:28],
                               "PT"
)
C5_e_means_2 <- cbind( C5_e_means, rep( "0.5", length( C5_e_means$mad_z ) ) )
colnames( C5_e_means_2 ) <- c( "quintile",
                               colnames( C5_e_means_2 )[2:28],
                               "PT"
)
C6_e_means_2 <- cbind( C6_e_means, rep( "1.0", length( C6_e_means$mad_z ) ) )
colnames( C6_e_means_2 ) <- c( "quintile",
                               colnames( C6_e_means_2 )[2:28],
                               "PT"
)


C_e_means <- rbind(
  C1_e_means_2, C2_e_means_2, C3_e_means_2,
  C4_e_means_2, C5_e_means_2, C6_e_means_2
)

# Plot all outcomes for "escitalopram" by reaction time quintile #####

ggplot(C_e_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, MADRS")

ggplot(C_e_means,aes(x=quintile,y=ta_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, Totasec")

ggplot(C_e_means,aes(x=quintile,y=f1_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, F1 score (mood)")

ggplot(C_e_means,aes(x=quintile,y=f2_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, F2 score (cognitive)")

ggplot(C_e_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, F3 score (neurovegetative)")

ggplot(C_e_means,aes(x=quintile,y=sui_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, suicidality score")

ggplot(C_e_means,aes(x=quintile,y=a1,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, dry mouth")

ggplot(C_e_means,aes(x=quintile,y=a2,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, drowsiness")

ggplot(C_e_means,aes(x=quintile,y=a3,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, insomnia")

ggplot(C_e_means,aes(x=quintile,y=a4,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, blurred vision")

ggplot(C_e_means,aes(x=quintile,y=a5,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, headache")

ggplot(C_e_means,aes(x=quintile,y=a6,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, constipation")

ggplot(C_e_means,aes(x=quintile,y=a7,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, diarrhoea")

ggplot(C_e_means,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, increased appetite")

ggplot(C_e_means,aes(x=quintile,y=a9,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, decreased appetite")

ggplot(C_e_means,aes(x=quintile,y=a10,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, nausea / vomiting")

ggplot(C_e_means,aes(x=quintile,y=a11,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, urination problems")

ggplot(C_e_means,aes(x=quintile,y=a12,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, sexual dysfunction")

ggplot(C_e_means,aes(x=quintile,y=a13,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, palpitations")

ggplot(C_e_means,aes(x=quintile,y=a14,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, light-headedness")

ggplot(C_e_means,aes(x=quintile,y=a15,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, room-spinning")

ggplot(C_e_means,aes(x=quintile,y=a16,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, sweating")

ggplot(C_e_means,aes(x=quintile,y=a17,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, increased body temperature")

ggplot(C_e_means,aes(x=quintile,y=a18,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, tremor")

ggplot(C_e_means,aes(x=quintile,y=a19,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, disorientation")

ggplot(C_e_means,aes(x=quintile,y=a20,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, yawning")

ggplot(C_e_means,aes(x=quintile,y=a21,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, escitalopram-only, weight gain")


##### Get means for things - Reaction Time; nortriptyline #####

C1_n_means <- GENDEP_nortrip %>% 
  group_by( C1_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C2_n_means <- GENDEP_nortrip %>% 
  group_by( C2_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C3_n_means <- GENDEP_nortrip %>% 
  group_by( C3_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C4_n_means <- GENDEP_nortrip %>% 
  group_by( C4_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C5_n_means <- GENDEP_nortrip %>% 
  group_by( C5_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C6_n_means <- GENDEP_nortrip %>% 
  group_by( C6_quintile ) %>%
  summarize( 
    mad_z = mean( zmadrs, na.rm=TRUE ),
    ta_z = mean( ztotasec, na.rm=TRUE ),
    a1 = mean( asec1wk, na.rm=TRUE ),
    a2 = mean( asec2wk, na.rm=TRUE ),
    a3 = mean( asec3wk, na.rm=TRUE ),
    a4 = mean( asec4wk, na.rm=TRUE ),
    a5 = mean( asec5wk, na.rm=TRUE ),
    a6 = mean( asec6wk, na.rm=TRUE ),
    a7 = mean( asec7wk, na.rm=TRUE ),
    a8 = mean( asec8wk, na.rm=TRUE ),
    a9 = mean( asec9wk, na.rm=TRUE ),
    a10 = mean( asec10wk, na.rm=TRUE ),
    a11 = mean( asec11wk, na.rm=TRUE ),
    a12 = mean( asec12wk, na.rm=TRUE ),
    a13 = mean( asec13wk, na.rm=TRUE ),
    a14 = mean( asec14wk, na.rm=TRUE ),
    a15 = mean( asec15wk, na.rm=TRUE ),
    a16 = mean( asec16wk, na.rm=TRUE ),
    a17 = mean( asec17wk, na.rm=TRUE ),
    a18 = mean( asec18wk, na.rm=TRUE ),
    a19 = mean( asec19wk, na.rm=TRUE ),
    a20 = mean( asec20wk, na.rm=TRUE ),
    a21 = mean( asec21wk, na.rm=TRUE ),
    f1_z = mean( zf1score, na.rm=TRUE ),
    f2_z = mean( zf2score, na.rm=TRUE ),
    f3_z = mean( zf3score, na.rm=TRUE ),
    sui_z = mean( suiscore, na.rm=TRUE )
  )

C1_n_means_2 <- cbind( C1_n_means, rep( "0.0001", length( C1_n_means$mad_z ) ) )
colnames( C1_n_means_2 ) <- c( "quintile",
                               colnames( C1_n_means_2 )[2:28],
                               "PT"
)
C2_n_means_2 <- cbind( C2_n_means, rep( "0.01", length( C2_n_means$mad_z ) ) )
colnames( C2_n_means_2 ) <- c( "quintile",
                               colnames( C2_n_means_2 )[2:28],
                               "PT"
)
C3_n_means_2 <- cbind( C3_n_means, rep( "0.05", length( C3_n_means$mad_z ) ) )
colnames( C3_n_means_2 ) <- c( "quintile",
                               colnames( C3_n_means_2 )[2:28],
                               "PT"
)
C4_n_means_2 <- cbind( C4_n_means, rep( "0.1", length( C4_n_means$mad_z ) ) )
colnames( C4_n_means_2 ) <- c( "quintile",
                               colnames( C4_n_means_2 )[2:28],
                               "PT"
)
C5_n_means_2 <- cbind( C5_n_means, rep( "0.5", length( C5_n_means$mad_z ) ) )
colnames( C5_n_means_2 ) <- c( "quintile",
                               colnames( C5_n_means_2 )[2:28],
                               "PT"
)
C6_n_means_2 <- cbind( C6_n_means, rep( "1.0", length( C6_n_means$mad_z ) ) )
colnames( C6_n_means_2 ) <- c( "quintile",
                               colnames( C6_n_means_2 )[2:28],
                               "PT"
)


C_n_means <- rbind(
  C1_n_means_2, C2_n_means_2, C3_n_means_2,
  C4_n_means_2, C5_n_means_2, C6_n_means_2
)

# Plot all outcomes for "nortriptyline" by reaction time quintile #####

ggplot(C_n_means,aes(x=quintile,y=mad_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, MADRS")

ggplot(C_n_means,aes(x=quintile,y=ta_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, Totasec")

ggplot(C_n_means,aes(x=quintile,y=f1_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, F1 score (mood)")

ggplot(C_n_means,aes(x=quintile,y=f2_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, F2 score (cognitive)")

ggplot(C_n_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, F3 score (neurovegetative)")

ggplot(C_n_means,aes(x=quintile,y=sui_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, suicidality score")

ggplot(C_n_means,aes(x=quintile,y=a1,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, dry mouth")

ggplot(C_n_means,aes(x=quintile,y=a2,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, drowsiness")

ggplot(C_n_means,aes(x=quintile,y=a3,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, insomnia")

ggplot(C_n_means,aes(x=quintile,y=a4,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, blurred vision")

ggplot(C_n_means,aes(x=quintile,y=a5,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, headache")

ggplot(C_n_means,aes(x=quintile,y=a6,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, constipation")

ggplot(C_n_means,aes(x=quintile,y=a7,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, diarrhoea")

ggplot(C_n_means,aes(x=quintile,y=a8,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, increased appetite")

ggplot(C_n_means,aes(x=quintile,y=a9,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, decreased appetite")

ggplot(C_n_means,aes(x=quintile,y=a10,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, nausea / vomiting")

ggplot(C_n_means,aes(x=quintile,y=a11,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, urination problems")

ggplot(C_n_means,aes(x=quintile,y=a12,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, sexual dysfunction")

ggplot(C_n_means,aes(x=quintile,y=a13,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, palpitations")

ggplot(C_n_means,aes(x=quintile,y=a14,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, light-headedness")

ggplot(C_n_means,aes(x=quintile,y=a15,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, room-spinning")

ggplot(C_n_means,aes(x=quintile,y=a16,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, sweating")

ggplot(C_n_means,aes(x=quintile,y=a17,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, increased body temperature")

ggplot(C_n_means,aes(x=quintile,y=a18,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, tremor")

ggplot(C_n_means,aes(x=quintile,y=a19,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, disorientation")

ggplot(C_n_means,aes(x=quintile,y=a20,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, yawning")

ggplot(C_n_means,aes(x=quintile,y=a21,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, nortriptyline-only, weight gain")

ggplot(C_means,aes(x=quintile,y=f3_z,fill=PT))+
  geom_col(position="dodge")+
  scale_fill_manual( values=c(
    "#006064", "#00838F", "#0097A7", 
    "#00ACC1", "#00BCD4", "#26C6DA") 
  )+
  labs(title="Reaction Time PS, F3 score (neurovegetative)")

##### SURVIVAL ANALYSIS - DROPOUT OF FIRST OR SECOND DRUG #####

# SOME OF THE "NON-DROPOUTS" HERE ARE SWITCHERS, THUS
# WILL HAVE "DROPPED OUT" OF THEIR FIRST PRESCRIBED DRUG

j <- vector(mode="numeric",length=0)
count <- 0

coeffs_both_2 <- vector(mode="list", length=0)
coeffs_esc_2 <- vector(mode="list", length=0)
coeffs_nor_2 <- vector(mode="list", length=0)
HRs_both_2 <- vector(mode="list", length=0)
HRs_esc_2 <- vector(mode="list", length=0)
HRs_nor_2 <- vector(mode="list", length=0)
models_both_2 <- vector(mode = "list", length = 31)
models_esc_2 <- vector(mode = "list", length = 31)
models_nor_2 <- vector(mode = "list", length = 31)

coeffs_both_ns_2 <- vector(mode="list", length=0)
coeffs_esc_ns_2 <- vector(mode="list", length=0)
coeffs_nor_ns_2 <- vector(mode="list", length=0)
HRs_both_ns_2 <- vector(mode="list", length=0)
HRs_esc_ns_2 <- vector(mode="list", length=0)
HRs_nor_ns_2 <- vector(mode="list", length=0)
models_both_ns_2 <- vector(mode = "list", length = 31)
models_esc_ns_2 <- vector(mode = "list", length = 31)
models_nor_ns_2 <- vector(mode = "list", length = 31)

coeffs_both_rand_2 <- vector(mode="list", length=0)
coeffs_esc_rand_2 <- vector(mode="list", length=0)
coeffs_nor_rand_2 <- vector(mode="list", length=0)
HRs_both_rand_2 <- vector(mode="list", length=0)
HRs_esc_rand_2 <- vector(mode="list", length=0)
HRs_nor_rand_2 <- vector(mode="list", length=0)
models_both_rand_2 <- vector(mode = "list", length = 31)
models_esc_rand_2 <- vector(mode = "list", length = 31)
models_nor_rand_2 <- vector(mode = "list", length = 31)

coeffs_both_rand_ns_2 <- vector(mode="list", length=0)
coeffs_esc_rand_ns_2 <- vector(mode="list", length=0)
coeffs_nor_rand_ns_2 <- vector(mode="list", length=0)
HRs_both_rand_ns_2 <- vector(mode="list", length=0)
HRs_esc_rand_ns_2 <- vector(mode="list", length=0)
HRs_nor_rand_ns_2 <- vector(mode="list", length=0)
models_both_rand_ns_2 <- vector(mode = "list", length = 31)
models_esc_rand_ns_2 <- vector(mode = "list", length = 31)
models_nor_rand_ns_2 <- vector(mode = "list", length = 31)

# IMPORTANT: THIS ANALYSIS IS LOOKING "DROPOUT" AS
# "DROPPED OUT OF STUDY AFTER ONE OR TWO DRUGS,
# AND THEIR FIRST DRUG WAS X..."

for (i in PRSs){
  
  count <- count+1
  j[i] <- which( colnames(GEND_end)==i ) # store the column number for the PRS in vector j
  
  GEND_end$PRS <- GEND_end[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  GE_nonswitch$PRS <- GE_nonswitch[,j[i]]
  GE_ran$PRS <- GE_ran[,j[i]]
  GE_ran_ns$PRS <- GE_ran_ns[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GEND_end) # Fit model
  
  coeffs_both_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_2[[count]] <- coxmod
  
  GEND_end$PRS<-NULL
  
  # Analysis on non-switchers only (dropout means dropout)
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_nonswitch) # Fit model
  
  coeffs_both_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_ns_2[[count]] <- coxmod
  
  GE_nonswitch$PRS<-NULL
  
  # Analysis on randomized only 
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran) # Fit model
  
  coeffs_both_rand_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_rand_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_rand_2[[count]] <- coxmod
  
  GE_ran$PRS<-NULL
  
  # Analysis on randomized, non-switcheres only (dropout means dropout)
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_ns) # Fit model
  
  coeffs_both_rand_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_rand_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_rand_ns_2[[count]] <- coxmod
  
  GE_ran_ns$PRS<-NULL
  
  #################### ESCITALOPRAM
  
  j[i] <- which( colnames(GE_esc)==i ) # store the column number for the PRS in vector j
  
  GE_esc$PRS <- GE_esc[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_esc) # Fit model
  
  coeffs_esc_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_2[[count]] <- coxmod
  
  GE_esc$PRS<-NULL
  
  # Analysis on non-switchers only (dropout means dropout)
  
  GE_ns_esc$PRS <- GE_ns_esc[,j[i]] # Add a column for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_ns_esc) # Fit model
  
  coeffs_esc_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_ns_2[[count]] <- coxmod
  
  GE_ns_esc$PRS<-NULL
  
  # Analysis on randomized only 
  
  GE_ran_e$PRS <- GE_ran_e[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_e) # Fit model
  
  coeffs_esc_rand_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_rand_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_rand_2[[count]] <- coxmod
  
  GE_ran_e$PRS<-NULL
  
  # Analysis on randomized, non-switcheres only (dropout means dropout)
  
  GE_ran_ns_e$PRS <- GE_ran_ns_e[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_ns_e) # Fit model
  
  coeffs_esc_rand_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_rand_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_rand_ns_2[[count]] <- coxmod
  
  GE_ran_ns_e$PRS<-NULL
  
  ############## NORTRIPTYLINE
  
  j[i] <- which( colnames(GE_nor)==i ) # store the column number for the PRS in vector j
  
  GE_nor$PRS <- GE_nor[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout) ~
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_nor) # Fit model
  
  coeffs_nor_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_2[[count]] <- coxmod
  
  GE_nor$PRS<-NULL
  
  # Analysis on non-switchers only
  
  GE_ns_nor$PRS <- GE_ns_nor[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_ns_nor) # Fit model
  
  coeffs_nor_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_ns_2[[count]] <- coxmod
  
  GE_ns_nor$PRS<-NULL
  
  # Analysis on randomized only 
  
  GE_ran_n$PRS <- GE_ran_n[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_n) # Fit model
  
  coeffs_nor_rand_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_rand_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_rand_2[[count]] <- coxmod
  
  GE_ran_n$PRS<-NULL
  
  # Analysis on randomized, non-switcheres only (dropout means dropout)
  
  GE_ran_ns_n$PRS <- GE_ran_ns_n[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_ns_n) # Fit model
  
  coeffs_nor_rand_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_rand_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_rand_ns_2[[count]] <- coxmod
  
  GE_ran_ns_n$PRS<-NULL
  
}

models_both_2[[31]] <- coxph(Surv(week, dropout) ~ 
                               drug + cage + sex +
                               zblmd +
                               frailty(centreid),
                             data = GEND_end)

models_esc_2[[31]] <- coxph(Surv(week, dropout) ~ 
                              cage + sex +
                              zblmd +
                              frailty(centreid),
                            data = GE_esc)

models_nor_2[[31]] <- coxph(Surv(week, dropout) ~ 
                              cage + sex +
                              zblmd +
                              frailty(centreid),
                            data = GE_nor)

models_both_ns_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                  drug + cage + sex +
                                  zblmd +
                                  frailty(centreid),
                                data = GE_nonswitch)

models_esc_ns_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                 cage + sex +
                                 zblmd +
                                 frailty(centreid),
                               data = GE_ns_esc)

models_nor_ns_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                 cage + sex +
                                 zblmd +
                                 frailty(centreid),
                               data = GE_ns_nor)

models_both_rand_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                    drug + cage + sex +
                                    zblmd +
                                    frailty(centreid),
                                  data = GE_ran)

models_esc_rand_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                   cage + sex +
                                   zblmd +
                                   frailty(centreid),
                                 data = GE_ran_e)

models_nor_rand_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                   cage + sex +
                                   zblmd +
                                   frailty(centreid),
                                 data = GE_ran_n)

models_both_rand_ns_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                       drug + cage + sex +
                                       zblmd +
                                       frailty(centreid),
                                     data = GE_ran_ns)

models_esc_rand_ns_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                      cage + sex +
                                      zblmd +
                                      frailty(centreid),
                                    data = GE_ran_ns_e)

models_nor_rand_ns_2[[31]] <- coxph(Surv(week, dropout) ~ 
                                      cage + sex +
                                      zblmd +
                                      frailty(centreid),
                                    data = GE_ran_ns_n)

names(models_both_2) <- PRSs_31
names(models_esc_2) <- PRSs_31
names(models_nor_2) <- PRSs_31

names(models_both_ns_2) <- PRSs_31
names(models_esc_ns_2) <- PRSs_31
names(models_nor_ns_2) <- PRSs_31

names(models_both_rand_2) <- PRSs_31
names(models_esc_rand_2) <- PRSs_31
names(models_nor_rand_2) <- PRSs_31

names(models_both_rand_ns_2) <- PRSs_31
names(models_esc_rand_ns_2) <- PRSs_31
names(models_nor_rand_ns_2) <- PRSs_31

# Look at explained variance:

# N.B. cod = coefficient of determination
#      mer = measure of explained randomness
#      mev = measure of explained variation

cod_both <- c()
mer_both <- c()
mev_both <- c()
cod_esc <- c()
mer_esc <- c()
mev_esc <- c()
cod_nor <- c()
mer_nor <- c()
mev_nor <- c()

cod_both_ns <- c()
mer_both_ns <- c()
mev_both_ns <- c()
cod_esc_ns <- c()
mer_esc_ns <- c()
mev_esc_ns <- c()
cod_nor_ns <- c()
mer_nor_ns <- c()
mev_nor_ns <- c()

cod_both_ran <- c()
mer_both_ran <- c()
mev_both_ran <- c()
cod_esc_ran <- c()
mer_esc_ran <- c()
mev_esc_ran <- c()
cod_nor_ran <- c()
mer_nor_ran <- c()
mev_nor_ran <- c()

cod_both_ns_ran <- c()
mer_both_ns_ran <- c()
mev_both_ns_ran <- c()
cod_esc_ns_ran <- c()
mer_esc_ns_ran <- c()
mev_esc_ns_ran <- c()
cod_nor_ns_ran <- c()
mer_nor_ns_ran <- c()
mev_nor_ns_ran <- c()

for (i in 1:31){
  cod_both[i] <- rsq( models_both_2[[i]] )$cod
  mer_both[i] <- rsq( models_both_2[[i]] )$mer
  mev_both[i] <- rsq( models_both_2[[i]] )$mev
  cod_esc[i] <- rsq( models_esc_2[[i]] )$cod
  mer_esc[i] <- rsq( models_esc_2[[i]] )$mer
  mev_esc[i] <- rsq( models_esc_2[[i]] )$mev
  cod_nor[i] <- rsq( models_nor_2[[i]] )$cod
  mer_nor[i] <- rsq( models_nor_2[[i]] )$mer
  mev_nor[i] <- rsq( models_nor_2[[i]] )$mev
  
  cod_both_ns[i] <- rsq( models_both_ns_2[[i]] )$cod
  mer_both_ns[i] <- rsq( models_both_ns_2[[i]] )$mer
  mev_both_ns[i] <- rsq( models_both_ns_2[[i]] )$mev
  cod_esc_ns[i] <- rsq( models_esc_ns_2[[i]] )$cod
  mer_esc_ns[i] <- rsq( models_esc_ns_2[[i]] )$mer
  mev_esc_ns[i] <- rsq( models_esc_ns_2[[i]] )$mev
  cod_nor_ns[i] <- rsq( models_nor_ns_2[[i]] )$cod
  mer_nor_ns[i] <- rsq( models_nor_ns_2[[i]] )$mer
  mev_nor_ns[i] <- rsq( models_nor_ns_2[[i]] )$mev
  
  cod_both_ran[i] <- rsq( models_both_rand_2[[i]] )$cod
  mer_both_ran[i] <- rsq( models_both_rand_2[[i]] )$mer
  mev_both_ran[i] <- rsq( models_both_rand_2[[i]] )$mev
  cod_esc_ran[i] <- rsq( models_esc_rand_2[[i]] )$cod
  mer_esc_ran[i] <- rsq( models_esc_rand_2[[i]] )$mer
  mev_esc_ran[i] <- rsq( models_esc_rand_2[[i]] )$mev
  cod_nor_ran[i] <- rsq( models_nor_rand_2[[i]] )$cod
  mer_nor_ran[i] <- rsq( models_nor_rand_2[[i]] )$mer
  mev_nor_ran[i] <- rsq( models_nor_rand_2[[i]] )$mev
  
  cod_both_ns_ran[i] <- rsq( models_both_rand_ns_2[[i]] )$cod
  mer_both_ns_ran[i] <- rsq( models_both_rand_ns_2[[i]] )$mer
  mev_both_ns_ran[i] <- rsq( models_both_rand_ns_2[[i]] )$mev
  cod_esc_ns_ran[i] <- rsq( models_esc_rand_ns_2[[i]] )$cod
  mer_esc_ns_ran[i] <- rsq( models_esc_rand_ns_2[[i]] )$mer
  mev_esc_ns_ran[i] <- rsq( models_esc_rand_ns_2[[i]] )$mev
  cod_nor_ns_ran[i] <- rsq( models_nor_rand_ns_2[[i]] )$cod
  mer_nor_ns_ran[i] <- rsq( models_nor_rand_ns_2[[i]] )$mer
  mev_nor_ns_ran[i] <- rsq( models_nor_rand_ns_2[[i]] )$mev
}

rsquared_2 <- data.frame(cbind(cod_both,
                               mer_both,
                               mev_both,
                               cod_esc,
                               mer_esc,
                               mev_esc,
                               cod_nor,
                               mer_nor,
                               mev_nor),
                         row.names=PRSs_31)

rsquared_ns_2 <- data.frame(cbind(cod_both_ns,
                                  mer_both_ns,
                                  mev_both_ns,
                                  cod_esc_ns,
                                  mer_esc_ns,
                                  mev_esc_ns,
                                  cod_nor_ns,
                                  mer_nor_ns,
                                  mev_nor_ns),
                            row.names=PRSs_31)

rsquared_2_r <- data.frame(cbind(cod_both_ran,
                                 mer_both_ran,
                                 mev_both_ran,
                                 cod_esc_ran,
                                 mer_esc_ran,
                                 mev_esc_ran,
                                 cod_nor_ran,
                                 mer_nor_ran,
                                 mev_nor_ran),
                           row.names=PRSs_31)

rsquared_ns_r_2 <- data.frame(cbind(cod_both_ns_ran,
                                    mer_both_ns_ran,
                                    mev_both_ns_ran,
                                    cod_esc_ns_ran,
                                    mer_esc_ns_ran,
                                    mev_esc_ns_ran,
                                    cod_nor_ns_ran,
                                    mer_nor_ns_ran,
                                    mev_nor_ns_ran),
                              row.names=PRSs_31)

ggplot(rsquared_2,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

ggplot(rsquared_ns_2,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

ggplot(rsquared_2_r,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

ggplot(rsquared_ns_r_2,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

##### Test cox model assumptions #####

coxmod <- coxph(Surv(week, dropout) ~ 
                  CZ_0_0001 + 
                  zblmd +
                  drug + cage + sex +
                  frailty(centreid),
                data = GEND_end)

cox.zph(coxmod) # p-values should all be non-significant
ggcoxdiagnostics(coxmod, 
                 type = "dfbeta", 
                 linear.predictions = FALSE)
ggcoxdiagnostics(coxmod, 
                 type = "deviance", 
                 linear.predictions = FALSE) # Should look symmetrical around 0
ggcoxfunctional(Surv(week, dropout) ~ 
                  CZ_0_0001 + 
                  zblmd +
                  cage +
                  frailty(centreid),
                data = GEND_end)

##### Output explained variances for survival analysis #####

surv_both_betas_2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_2 <- vector(mode = "list", length = length(PRSs))

surv_both_betas_ns2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_ns2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_ns2 <- vector(mode = "list", length = length(PRSs))

surv_both_betas_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_ran2 <- vector(mode = "list", length = length(PRSs))

surv_both_betas_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_ns_ran2 <- vector(mode = "list", length = length(PRSs))

for (i in 1:length(coeffs_both_2)){
  surv_both_betas_2[[i]] <- coeffs_both_2[[i]][1,]
  surv_both_betas_ns2[[i]] <- coeffs_both_ns_2[[i]][1,]
  surv_both_betas_ran2[[i]] <- coeffs_both_rand_2[[i]][1,]
  surv_both_betas_ns_ran2[[i]] <- coeffs_both_rand_ns_2[[i]][1,]
}

for (i in 1:length(coeffs_esc_2)){
  surv_esc_betas_2[[i]] <- coeffs_esc_2[[i]][1,]
  surv_esc_betas_ns2[[i]] <- coeffs_esc_ns_2[[i]][1,]
  surv_esc_betas_ran2[[i]] <- coeffs_esc_rand_2[[i]][1,]
  surv_esc_betas_ns_ran2[[i]] <- coeffs_esc_rand_ns_2[[i]][1,]
}

for (i in 1:length(coeffs_nor_2)){
  surv_nor_betas_2[[i]] <- coeffs_nor_2[[i]][1,]
  surv_nor_betas_ns2[[i]] <- coeffs_nor_ns_2[[i]][1,]
  surv_nor_betas_ran2[[i]] <- coeffs_nor_rand_2[[i]][1,]
  surv_nor_betas_ns_ran2[[i]] <- coeffs_nor_rand_ns_2[[i]][1,]
}

names(surv_both_betas_2) <- PRSs
names(surv_esc_betas_2) <- PRSs
names(surv_nor_betas_2) <- PRSs

names(surv_both_betas_ns2) <- PRSs
names(surv_esc_betas_ns2) <- PRSs
names(surv_nor_betas_ns2) <- PRSs

names(surv_both_betas_ran2) <- PRSs
names(surv_esc_betas_ran2) <- PRSs
names(surv_nor_betas_ran2) <- PRSs

names(surv_both_betas_ns_ran2) <- PRSs
names(surv_esc_betas_ns_ran2) <- PRSs
names(surv_nor_betas_ns_ran2) <- PRSs

## Get hazard ratios (exponentiated coefficients # exp(coef)) #####

surv_both_HRs_2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_2 <- vector(mode = "list", length = length(PRSs))

surv_both_HRs_ns2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_ns2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_ns2 <- vector(mode = "list", length = length(PRSs))

surv_both_HRs_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_ran2 <- vector(mode = "list", length = length(PRSs))

surv_both_HRs_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_ns_ran2 <- vector(mode = "list", length = length(PRSs))

for (i in 1:length(HRs_both_2)){
  surv_both_HRs_2[[i]] <- HRs_both_2[[i]][1,]
  surv_both_HRs_ns2[[i]] <- HRs_both_ns_2[[i]][1,]
  surv_both_HRs_ran2[[i]] <- HRs_both_rand_2[[i]][1,]
  surv_both_HRs_ns_ran2[[i]] <- HRs_both_rand_ns_2[[i]][1,]
}

for (i in 1:length(HRs_esc_2)){
  surv_esc_HRs_2[[i]] <- HRs_esc_2[[i]][1,]
  surv_esc_HRs_ns2[[i]] <- HRs_esc_ns_2[[i]][1,]
  surv_esc_HRs_ran2[[i]] <- HRs_esc_rand_2[[i]][1,]
  surv_esc_HRs_ns_ran2[[i]] <- HRs_esc_rand_ns_2[[i]][1,]
}

for (i in 1:length(HRs_nor_2)){
  surv_nor_HRs_2[[i]] <- HRs_nor_2[[i]][1,]
  surv_nor_HRs_ns2[[i]] <- HRs_nor_ns_2[[i]][1,]
  surv_nor_HRs_ran2[[i]] <- HRs_nor_rand_2[[i]][1,]
  surv_nor_HRs_ns_ran2[[i]] <- HRs_nor_rand_ns_2[[i]][1,]
}

names(surv_both_HRs_2) <- PRSs
names(surv_esc_HRs_2) <- PRSs
names(surv_nor_HRs_2) <- PRSs

names(surv_both_HRs_ns2) <- PRSs
names(surv_esc_HRs_ns2) <- PRSs
names(surv_nor_HRs_ns2) <- PRSs

names(surv_both_HRs_ran2) <- PRSs
names(surv_esc_HRs_ran2) <- PRSs
names(surv_nor_HRs_ran2) <- PRSs

names(surv_both_HRs_ns_ran2) <- PRSs
names(surv_esc_HRs_ns_ran2) <- PRSs
names(surv_nor_HRs_ns_ran2) <- PRSs


## Make survival output file #####

surv_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_2))[,1],
    t(as.data.frame(surv_both_betas_2))[,3],
    t(as.data.frame(surv_both_betas_2))[,6],
    t(as.data.frame(surv_esc_betas_2))[,1],
    t(as.data.frame(surv_esc_betas_2))[,3],
    t(as.data.frame(surv_esc_betas_2))[,6],
    t(as.data.frame(surv_nor_betas_2))[,1],
    t(as.data.frame(surv_nor_betas_2))[,3],
    t(as.data.frame(surv_nor_betas_2))[,6]
  )
)

colnames( surv_beta_matrix_2 ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_2))[,1],
    t(as.data.frame(surv_both_HRs_2))[,3],
    t(as.data.frame(surv_both_HRs_2))[,4],
    t(as.data.frame(surv_esc_HRs_2))[,1],
    t(as.data.frame(surv_esc_HRs_2))[,3],
    t(as.data.frame(surv_esc_HRs_2))[,4],
    t(as.data.frame(surv_nor_HRs_2))[,1],
    t(as.data.frame(surv_nor_HRs_2))[,3],
    t(as.data.frame(surv_nor_HRs_2))[,4]
  )
)

colnames( surv_HR_matrix_2 ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)


sink( "survival_betas.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_2, col.names=TRUE, sep="," )

cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_2, col.names=TRUE, sep="," )

cat("\nSurvival_expvar\n,")
write.table( rsquared_2,
             col.names=TRUE, sep="," )

sink()

## Make survival output file, non-switchers only #####

surv_beta_matrix_ns_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_ns2))[,1],
    t(as.data.frame(surv_both_betas_ns2))[,3],
    t(as.data.frame(surv_both_betas_ns2))[,6],
    t(as.data.frame(surv_esc_betas_ns2))[,1],
    t(as.data.frame(surv_esc_betas_ns2))[,3],
    t(as.data.frame(surv_esc_betas_ns2))[,6],
    t(as.data.frame(surv_nor_betas_ns2))[,1],
    t(as.data.frame(surv_nor_betas_ns2))[,3],
    t(as.data.frame(surv_nor_betas_ns2))[,6]
  )
)

colnames( surv_beta_matrix_ns_2 ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_ns_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_ns2))[,1],
    t(as.data.frame(surv_both_HRs_ns2))[,3],
    t(as.data.frame(surv_both_HRs_ns2))[,4],
    t(as.data.frame(surv_esc_HRs_ns2))[,1],
    t(as.data.frame(surv_esc_HRs_ns2))[,3],
    t(as.data.frame(surv_esc_HRs_ns2))[,4],
    t(as.data.frame(surv_nor_HRs_ns2))[,1],
    t(as.data.frame(surv_nor_HRs_ns2))[,3],
    t(as.data.frame(surv_nor_HRs_ns2))[,4]
  )
)

colnames( surv_HR_matrix_ns_2 ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)


sink( "survival_nonswitchers_betas.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_ns_2, col.names=TRUE, sep="," )

cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_ns_2, col.names=TRUE, sep="," )

cat("\nSurvival_expvar\n,")
write.table( rsquared_ns_2,
             col.names=TRUE, sep="," )

sink()

## Make survival output file, random only #####

surv_beta_matrix_ran <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_ran2))[,1],
    t(as.data.frame(surv_both_betas_ran2))[,3],
    t(as.data.frame(surv_both_betas_ran2))[,6],
    t(as.data.frame(surv_esc_betas_ran2))[,1],
    t(as.data.frame(surv_esc_betas_ran2))[,3],
    t(as.data.frame(surv_esc_betas_ran2))[,6],
    t(as.data.frame(surv_nor_betas_ran2))[,1],
    t(as.data.frame(surv_nor_betas_ran2))[,3],
    t(as.data.frame(surv_nor_betas_ran2))[,6]
  )
)

colnames( surv_beta_matrix_ran ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_ran <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_ran2))[,1],
    t(as.data.frame(surv_both_HRs_ran2))[,3],
    t(as.data.frame(surv_both_HRs_ran2))[,4],
    t(as.data.frame(surv_esc_HRs_ran2))[,1],
    t(as.data.frame(surv_esc_HRs_ran2))[,3],
    t(as.data.frame(surv_esc_HRs_ran2))[,4],
    t(as.data.frame(surv_nor_HRs_ran2))[,1],
    t(as.data.frame(surv_nor_HRs_ran2))[,3],
    t(as.data.frame(surv_nor_HRs_ran2))[,4]
  )
)

colnames( surv_HR_matrix_ran ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)

sink( "survival_randomized_betas.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_ran, col.names=TRUE, sep="," )

cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_ran, col.names=TRUE, sep="," )

cat("\nSurvival_expvar\n,")
write.table( rsquared_2_r,
             col.names=TRUE, sep="," )

sink()

## Make survival output file, randomized non-switchers only #####

surv_beta_matrix_ns_ran2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_ns_ran2))[,1],
    t(as.data.frame(surv_both_betas_ns_ran2))[,3],
    t(as.data.frame(surv_both_betas_ns_ran2))[,6],
    t(as.data.frame(surv_esc_betas_ns_ran2))[,1],
    t(as.data.frame(surv_esc_betas_ns_ran2))[,3],
    t(as.data.frame(surv_esc_betas_ns_ran2))[,6],
    t(as.data.frame(surv_nor_betas_ns_ran2))[,1],
    t(as.data.frame(surv_nor_betas_ns_ran2))[,3],
    t(as.data.frame(surv_nor_betas_ns_ran2))[,6]
  )
)

colnames( surv_beta_matrix_ns_ran2 ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_ns_ran <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_ns_ran2))[,1],
    t(as.data.frame(surv_both_HRs_ns_ran2))[,3],
    t(as.data.frame(surv_both_HRs_ns_ran2))[,4],
    t(as.data.frame(surv_esc_HRs_ns_ran2))[,1],
    t(as.data.frame(surv_esc_HRs_ns_ran2))[,3],
    t(as.data.frame(surv_esc_HRs_ns_ran2))[,4],
    t(as.data.frame(surv_nor_HRs_ns_ran2))[,1],
    t(as.data.frame(surv_nor_HRs_ns_ran2))[,3],
    t(as.data.frame(surv_nor_HRs_ns_ran2))[,4]
  )
)

colnames( surv_HR_matrix_ns_ran ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)

sink( "survival_rand_nonswitchers_betas.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_ns_ran2, col.names=TRUE, sep="," )

cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_ns_ran, col.names=TRUE, sep="," )

cat("\nSurvival_expvar\n,")
write.table( rsquared_ns_r_2,
             col.names=TRUE, sep="," )

sink()

## Make survival-only heatmap #####

surv_models_betas_2 <- c("surv_both_s_B", "surv_esc_s_B", "surv_nor_s_B")
surv_models_ses_2 <- c("surv_both_s_se", "surv_esc_s_se", "surv_nor_s_se")
surv_models_ps_2 <- c("surv_both_s_p", "surv_esc_s_p", "surv_nor_s_p")

surv_models_betas_columns_2 <- which( colnames(surv_beta_matrix_2)
                                      %in% surv_models_betas_2 )
surv_models_ses_columns_2 <- which( colnames(surv_beta_matrix_2)
                                    %in% surv_models_ses_2 )
surv_models_ps_columns_2 <- which( colnames(surv_beta_matrix_2)
                                   %in% surv_models_ps_2 )

surv_betas_only_2 <- cbind(surv_beta_matrix_2[,surv_models_betas_columns_2[1]],
                           surv_beta_matrix_2[,surv_models_betas_columns_2[2]],
                           surv_beta_matrix_2[,surv_models_betas_columns_2[3]])

colnames(surv_betas_only_2) <- c("surv_both", "surv_esc","surv_nor")


surv_betas_only_2 <- cbind(surv_beta_matrix_2[,1],
                           surv_beta_matrix_2[,4],
                           surv_beta_matrix_2[,7])

colnames(surv_betas_only_2) <- c("surv_both", "surv_esc", "surv_nor")

surv_se_only_2 <- cbind(surv_beta_matrix_2[,2],
                        surv_beta_matrix_2[,5],
                        surv_beta_matrix_2[,8])

colnames(surv_se_only_2) <- c("surv_both", "surv_esc","surv_nor")

surv_p_only_2 <- cbind(surv_beta_matrix_2[,3],
                       surv_beta_matrix_2[,6],
                       surv_beta_matrix_2[,9])

colnames(surv_p_only_2) <- c("surv_both", "surv_esc","surv_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "both",
                                                         "esc",
                                                         "nor")))

lowest_ps[which(surv_p_only_2<=0.05)]<-surv_p_only_2[which(surv_p_only_2<=0.05)]

heatmap.2( surv_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps,5),notecol="black",
           main="Dropout of study"
)


##### SURVIVAL ANALYSIS - DISCONTINUATION OF FIRST DRUG #####

j <- vector(mode="numeric",length=0)
count <- 0

coeffs_both_2 <- vector(mode="list", length=0)
coeffs_esc_2 <- vector(mode="list", length=0)
coeffs_nor_2 <- vector(mode="list", length=0)
models_both_2 <- vector(mode = "list", length = 31)
models_esc_2 <- vector(mode = "list", length = 31)
models_nor_2 <- vector(mode = "list", length = 31)

coeffs_both_ns_2 <- vector(mode="list", length=0)
coeffs_esc_ns_2 <- vector(mode="list", length=0)
coeffs_nor_ns_2 <- vector(mode="list", length=0)
models_both_ns_2 <- vector(mode = "list", length = 31)
models_esc_ns_2 <- vector(mode = "list", length = 31)
models_nor_ns_2 <- vector(mode = "list", length = 31)

coeffs_both_rand_2 <- vector(mode="list", length=0)
coeffs_esc_rand_2 <- vector(mode="list", length=0)
coeffs_nor_rand_2 <- vector(mode="list", length=0)
models_both_rand_2 <- vector(mode = "list", length = 31)
models_esc_rand_2 <- vector(mode = "list", length = 31)
models_nor_rand_2 <- vector(mode = "list", length = 31)

coeffs_both_rand_ns_2 <- vector(mode="list", length=0)
coeffs_esc_rand_ns_2 <- vector(mode="list", length=0)
coeffs_nor_rand_ns_2 <- vector(mode="list", length=0)
models_both_rand_ns_2 <- vector(mode = "list", length = 31)
models_esc_rand_ns_2 <- vector(mode = "list", length = 31)
models_nor_rand_ns_2 <- vector(mode = "list", length = 31)

# IMPORTANT: THIS ANALYSIS IS LOOKING "DROPOUT" AS
# "DROPPED OUT OF FIRST DRUG"

for (i in PRSs){
  
  count <- count+1
  j[i] <- which( colnames(GEND_end)==i ) # store the column number for the PRS in vector j
  
  GEND_end$PRS <- GEND_end[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  GE_nonswitch$PRS <- GE_nonswitch[,j[i]]
  GE_ran$PRS <- GE_ran[,j[i]]
  GE_ran_ns$PRS <- GE_ran_ns[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GEND_end) # Fit model
  
  coeffs_both_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_2[[count]] <- coxmod
  
  GEND_end$PRS<-NULL
  
  # Analysis on non-switchers only (dropout2 means dropout2)
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_nonswitch) # Fit model
  
  coeffs_both_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_ns_2[[count]] <- coxmod
  
  GE_nonswitch$PRS<-NULL
  
  # Analysis on randomized only 
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran) # Fit model
  
  coeffs_both_rand_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_rand_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_rand_2[[count]] <- coxmod
  
  GE_ran$PRS<-NULL
  
  # Analysis on randomized, non-switcheres only (dropout2 means dropout2)
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_ns) # Fit model
  
  coeffs_both_rand_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_both_rand_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_both_rand_ns_2[[count]] <- coxmod
  
  GE_ran_ns$PRS<-NULL
  
  #################### ESCITALOPRAM
  
  j[i] <- which( colnames(GE_esc)==i ) # store the column number for the PRS in vector j
  
  GE_esc$PRS <- GE_esc[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_esc) # Fit model
  
  coeffs_esc_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_2[[count]] <- coxmod
  
  GE_esc$PRS<-NULL
  
  # Analysis on non-switchers only (dropout2 means dropout2)
  
  GE_ns_esc$PRS <- GE_ns_esc[,j[i]] # Add a column for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_ns_esc) # Fit model
  
  coeffs_esc_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_ns_2[[count]] <- coxmod
  
  GE_ns_esc$PRS<-NULL
  
  # Analysis on randomized only 
  
  GE_ran_e$PRS <- GE_ran_e[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_e) # Fit model
  
  coeffs_esc_rand_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_rand_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_rand_2[[count]] <- coxmod
  
  GE_ran_e$PRS<-NULL
  
  # Analysis on randomized, non-switcheres only (dropout2 means dropout2)
  
  GE_ran_ns_e$PRS <- GE_ran_ns_e[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_ns_e) # Fit model
  
  coeffs_esc_rand_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_esc_rand_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_esc_rand_ns_2[[count]] <- coxmod
  
  GE_ran_ns_e$PRS<-NULL
  
  ############## NORTRIPTYLINE
  
  j[i] <- which( colnames(GE_nor)==i ) # store the column number for the PRS in vector j
  
  GE_nor$PRS <- GE_nor[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout2) ~
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_nor) # Fit model
  
  coeffs_nor_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_2[[count]] <- coxmod
  
  GE_nor$PRS<-NULL
  
  # Analysis on non-switchers only
  
  GE_ns_nor$PRS <- GE_ns_nor[,j[i]] # Add a column to GEND_end for current PRS, to refer to below
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    cage + sex +
                    frailty(centreid),
                  data = GE_ns_nor) # Fit model
  
  coeffs_nor_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_ns_2[[count]] <- coxmod
  
  GE_ns_nor$PRS<-NULL
  
  # Analysis on randomized only 
  
  GE_ran_n$PRS <- GE_ran_n[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_n) # Fit model
  
  coeffs_nor_rand_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_rand_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_rand_2[[count]] <- coxmod
  
  GE_ran_n$PRS<-NULL
  
  # Analysis on randomized, non-switcheres only (dropout2 means dropout2)
  
  GE_ran_ns_n$PRS <- GE_ran_ns_n[,j[i]]
  
  coxmod <- coxph(Surv(week, dropout2) ~ 
                    PRS + 
                    zblmd +
                    drug + cage + sex +
                    frailty(centreid),
                  data = GE_ran_ns_n) # Fit model
  
  coeffs_nor_rand_ns_2[[i]] <- summary(coxmod)$coefficients #  Negative no. = stay on drugs longer
  HRs_nor_rand_ns_2[[i]] <- summary(coxmod)$conf.int
  
  models_nor_rand_ns_2[[count]] <- coxmod
  
  GE_ran_ns_n$PRS<-NULL
  
}

models_both_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                               drug + cage + sex +
                               zblmd +
                               frailty(centreid),
                             data = GEND_end)

models_esc_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                              cage + sex +
                              zblmd +
                              frailty(centreid),
                            data = GE_esc)

models_nor_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                              cage + sex +
                              zblmd +
                              frailty(centreid),
                            data = GE_nor)

models_both_ns_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                  drug + cage + sex +
                                  zblmd +
                                  frailty(centreid),
                                data = GE_nonswitch)

models_esc_ns_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                 cage + sex +
                                 zblmd +
                                 frailty(centreid),
                               data = GE_ns_esc)

models_nor_ns_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                 cage + sex +
                                 zblmd +
                                 frailty(centreid),
                               data = GE_ns_nor)

models_both_rand_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                    drug + cage + sex +
                                    zblmd +
                                    frailty(centreid),
                                  data = GE_ran)

models_esc_rand_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                   cage + sex +
                                   zblmd +
                                   frailty(centreid),
                                 data = GE_ran_e)

models_nor_rand_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                   cage + sex +
                                   zblmd +
                                   frailty(centreid),
                                 data = GE_ran_n)

models_both_rand_ns_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                       drug + cage + sex +
                                       zblmd +
                                       frailty(centreid),
                                     data = GE_ran_ns)

models_esc_rand_ns_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                      cage + sex +
                                      zblmd +
                                      frailty(centreid),
                                    data = GE_ran_ns_e)

models_nor_rand_ns_2[[31]] <- coxph(Surv(week, dropout2) ~ 
                                      cage + sex +
                                      zblmd +
                                      frailty(centreid),
                                    data = GE_ran_ns_n)

names(models_both_2) <- PRSs_31
names(models_esc_2) <- PRSs_31
names(models_nor_2) <- PRSs_31

names(models_both_ns_2) <- PRSs_31
names(models_esc_ns_2) <- PRSs_31
names(models_nor_ns_2) <- PRSs_31

names(models_both_rand_2) <- PRSs_31
names(models_esc_rand_2) <- PRSs_31
names(models_nor_rand_2) <- PRSs_31

names(models_both_rand_ns_2) <- PRSs_31
names(models_esc_rand_ns_2) <- PRSs_31
names(models_nor_rand_ns_2) <- PRSs_31

# Look at explained variance:

# N.B. cod = coefficient of determination
#      mer = measure of explained randomness
#      mev = measure of explained variation

cod_both <- c()
mer_both <- c()
mev_both <- c()
cod_esc <- c()
mer_esc <- c()
mev_esc <- c()
cod_nor <- c()
mer_nor <- c()
mev_nor <- c()

cod_both_ns <- c()
mer_both_ns <- c()
mev_both_ns <- c()
cod_esc_ns <- c()
mer_esc_ns <- c()
mev_esc_ns <- c()
cod_nor_ns <- c()
mer_nor_ns <- c()
mev_nor_ns <- c()

cod_both_ran <- c()
mer_both_ran <- c()
mev_both_ran <- c()
cod_esc_ran <- c()
mer_esc_ran <- c()
mev_esc_ran <- c()
cod_nor_ran <- c()
mer_nor_ran <- c()
mev_nor_ran <- c()

cod_both_ns_ran <- c()
mer_both_ns_ran <- c()
mev_both_ns_ran <- c()
cod_esc_ns_ran <- c()
mer_esc_ns_ran <- c()
mev_esc_ns_ran <- c()
cod_nor_ns_ran <- c()
mer_nor_ns_ran <- c()
mev_nor_ns_ran <- c()

for (i in 1:31){
  cod_both[i] <- rsq( models_both_2[[i]] )$cod
  mer_both[i] <- rsq( models_both_2[[i]] )$mer
  mev_both[i] <- rsq( models_both_2[[i]] )$mev
  cod_esc[i] <- rsq( models_esc_2[[i]] )$cod
  mer_esc[i] <- rsq( models_esc_2[[i]] )$mer
  mev_esc[i] <- rsq( models_esc_2[[i]] )$mev
  cod_nor[i] <- rsq( models_nor_2[[i]] )$cod
  mer_nor[i] <- rsq( models_nor_2[[i]] )$mer
  mev_nor[i] <- rsq( models_nor_2[[i]] )$mev
  
  cod_both_ns[i] <- rsq( models_both_ns_2[[i]] )$cod
  mer_both_ns[i] <- rsq( models_both_ns_2[[i]] )$mer
  mev_both_ns[i] <- rsq( models_both_ns_2[[i]] )$mev
  cod_esc_ns[i] <- rsq( models_esc_ns_2[[i]] )$cod
  mer_esc_ns[i] <- rsq( models_esc_ns_2[[i]] )$mer
  mev_esc_ns[i] <- rsq( models_esc_ns_2[[i]] )$mev
  cod_nor_ns[i] <- rsq( models_nor_ns_2[[i]] )$cod
  mer_nor_ns[i] <- rsq( models_nor_ns_2[[i]] )$mer
  mev_nor_ns[i] <- rsq( models_nor_ns_2[[i]] )$mev
  
  cod_both_ran[i] <- rsq( models_both_rand_2[[i]] )$cod
  mer_both_ran[i] <- rsq( models_both_rand_2[[i]] )$mer
  mev_both_ran[i] <- rsq( models_both_rand_2[[i]] )$mev
  cod_esc_ran[i] <- rsq( models_esc_rand_2[[i]] )$cod
  mer_esc_ran[i] <- rsq( models_esc_rand_2[[i]] )$mer
  mev_esc_ran[i] <- rsq( models_esc_rand_2[[i]] )$mev
  cod_nor_ran[i] <- rsq( models_nor_rand_2[[i]] )$cod
  mer_nor_ran[i] <- rsq( models_nor_rand_2[[i]] )$mer
  mev_nor_ran[i] <- rsq( models_nor_rand_2[[i]] )$mev
  
  cod_both_ns_ran[i] <- rsq( models_both_rand_ns_2[[i]] )$cod
  mer_both_ns_ran[i] <- rsq( models_both_rand_ns_2[[i]] )$mer
  mev_both_ns_ran[i] <- rsq( models_both_rand_ns_2[[i]] )$mev
  cod_esc_ns_ran[i] <- rsq( models_esc_rand_ns_2[[i]] )$cod
  mer_esc_ns_ran[i] <- rsq( models_esc_rand_ns_2[[i]] )$mer
  mev_esc_ns_ran[i] <- rsq( models_esc_rand_ns_2[[i]] )$mev
  cod_nor_ns_ran[i] <- rsq( models_nor_rand_ns_2[[i]] )$cod
  mer_nor_ns_ran[i] <- rsq( models_nor_rand_ns_2[[i]] )$mer
  mev_nor_ns_ran[i] <- rsq( models_nor_rand_ns_2[[i]] )$mev
}

rsquared_2 <- data.frame(cbind(cod_both,
                               mer_both,
                               mev_both,
                               cod_esc,
                               mer_esc,
                               mev_esc,
                               cod_nor,
                               mer_nor,
                               mev_nor),
                         row.names=PRSs_31)

rsquared_ns_2 <- data.frame(cbind(cod_both_ns,
                                  mer_both_ns,
                                  mev_both_ns,
                                  cod_esc_ns,
                                  mer_esc_ns,
                                  mev_esc_ns,
                                  cod_nor_ns,
                                  mer_nor_ns,
                                  mev_nor_ns),
                            row.names=PRSs_31)

rsquared_2_r <- data.frame(cbind(cod_both_ran,
                                 mer_both_ran,
                                 mev_both_ran,
                                 cod_esc_ran,
                                 mer_esc_ran,
                                 mev_esc_ran,
                                 cod_nor_ran,
                                 mer_nor_ran,
                                 mev_nor_ran),
                           row.names=PRSs_31)

rsquared_ns_r_2 <- data.frame(cbind(cod_both_ns_ran,
                                    mer_both_ns_ran,
                                    mev_both_ns_ran,
                                    cod_esc_ns_ran,
                                    mer_esc_ns_ran,
                                    mev_esc_ns_ran,
                                    cod_nor_ns_ran,
                                    mer_nor_ns_ran,
                                    mev_nor_ns_ran),
                              row.names=PRSs_31)

ggplot(rsquared_2,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

ggplot(rsquared_ns_2,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

ggplot(rsquared_2_r,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

ggplot(rsquared_ns_r_2,aes(x=PRSs_31))+
  geom_point(aes(y=cod_both),col="navy")+
  geom_point(aes(y=mer_both),col="blue")+
  geom_point(aes(y=mev_both),col="steelblue")+
  geom_hline(yintercept=cod_both[31],col="navy")+
  geom_hline(yintercept=mer_both[31],col="blue")+
  geom_hline(yintercept=mev_both[31],col="steelblue")+
  geom_point(aes(y=cod_esc),col="magenta")+
  geom_point(aes(y=mer_esc),col="pink")+
  geom_point(aes(y=mev_esc),col="deeppink")+
  geom_hline(yintercept=cod_esc[31],col="magenta")+
  geom_hline(yintercept=mer_esc[31],col="pink")+
  geom_hline(yintercept=mev_esc[31],col="deeppink")+
  geom_point(aes(y=cod_nor),col="forestgreen")+
  geom_point(aes(y=mer_nor),col="yellowgreen")+
  geom_point(aes(y=mev_nor),col="green")+
  geom_hline(yintercept=cod_nor[31],col="forestgreen")+
  geom_hline(yintercept=mer_nor[31],col="yellowgreen")+
  geom_hline(yintercept=mev_nor[31],col="green")+
  theme(axis.text.x = element_text(angle = 90))+
  coord_cartesian(ylim=c(0,0.4))

##### Test cox model assumptions #####

coxmod <- coxph(Surv(week, dropout2) ~ 
                  CZ_0_0001 + 
                  zblmd +
                  drug + cage + sex +
                  frailty(centreid),
                data = GEND_end)

cox.zph(coxmod) # p-values should all be non-significant
ggcoxdiagnostics(coxmod, 
                 type = "dfbeta", 
                 linear.predictions = FALSE)
ggcoxdiagnostics(coxmod, 
                 type = "deviance", 
                 linear.predictions = FALSE) # Should look symmetrical around 0
ggcoxfunctional(Surv(week, dropout2) ~ 
                  CZ_0_0001 + 
                  zblmd +
                  cage +
                  frailty(centreid),
                data = GEND_end)

##### Output explained variances for survival analysis #####

surv_both_betas_2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_2 <- vector(mode = "list", length = length(PRSs))

surv_both_betas_ns2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_ns2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_ns2 <- vector(mode = "list", length = length(PRSs))

surv_both_betas_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_ran2 <- vector(mode = "list", length = length(PRSs))

surv_both_betas_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_betas_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_betas_ns_ran2 <- vector(mode = "list", length = length(PRSs))

for (i in 1:length(coeffs_both_2)){
  surv_both_betas_2[[i]] <- coeffs_both_2[[i]][1,]
  surv_both_betas_ns2[[i]] <- coeffs_both_ns_2[[i]][1,]
  surv_both_betas_ran2[[i]] <- coeffs_both_rand_2[[i]][1,]
  surv_both_betas_ns_ran2[[i]] <- coeffs_both_rand_ns_2[[i]][1,]
}

for (i in 1:length(coeffs_esc_2)){
  surv_esc_betas_2[[i]] <- coeffs_esc_2[[i]][1,]
  surv_esc_betas_ns2[[i]] <- coeffs_esc_ns_2[[i]][1,]
  surv_esc_betas_ran2[[i]] <- coeffs_esc_rand_2[[i]][1,]
  surv_esc_betas_ns_ran2[[i]] <- coeffs_esc_rand_ns_2[[i]][1,]
}

for (i in 1:length(coeffs_nor_2)){
  surv_nor_betas_2[[i]] <- coeffs_nor_2[[i]][1,]
  surv_nor_betas_ns2[[i]] <- coeffs_nor_ns_2[[i]][1,]
  surv_nor_betas_ran2[[i]] <- coeffs_nor_rand_2[[i]][1,]
  surv_nor_betas_ns_ran2[[i]] <- coeffs_nor_rand_ns_2[[i]][1,]
}

names(surv_both_betas_2) <- PRSs
names(surv_esc_betas_2) <- PRSs
names(surv_nor_betas_2) <- PRSs

names(surv_both_betas_ns2) <- PRSs
names(surv_esc_betas_ns2) <- PRSs
names(surv_nor_betas_ns2) <- PRSs

names(surv_both_betas_ran2) <- PRSs
names(surv_esc_betas_ran2) <- PRSs
names(surv_nor_betas_ran2) <- PRSs

names(surv_both_betas_ns_ran2) <- PRSs
names(surv_esc_betas_ns_ran2) <- PRSs
names(surv_nor_betas_ns_ran2) <- PRSs

## Get hazard ratios (exponentiated coefficients # exp(coef)) #####

surv_both_HRs_2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_2 <- vector(mode = "list", length = length(PRSs))

surv_both_HRs_ns2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_ns2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_ns2 <- vector(mode = "list", length = length(PRSs))

surv_both_HRs_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_ran2 <- vector(mode = "list", length = length(PRSs))

surv_both_HRs_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_esc_HRs_ns_ran2 <- vector(mode = "list", length = length(PRSs))
surv_nor_HRs_ns_ran2 <- vector(mode = "list", length = length(PRSs))

for (i in 1:length(HRs_both_2)){
  surv_both_HRs_2[[i]] <- HRs_both_2[[i]][1,]
  surv_both_HRs_ns2[[i]] <- HRs_both_ns_2[[i]][1,]
  surv_both_HRs_ran2[[i]] <- HRs_both_rand_2[[i]][1,]
  surv_both_HRs_ns_ran2[[i]] <- HRs_both_rand_ns_2[[i]][1,]
}

for (i in 1:length(HRs_esc_2)){
  surv_esc_HRs_2[[i]] <- HRs_esc_2[[i]][1,]
  surv_esc_HRs_ns2[[i]] <- HRs_esc_ns_2[[i]][1,]
  surv_esc_HRs_ran2[[i]] <- HRs_esc_rand_2[[i]][1,]
  surv_esc_HRs_ns_ran2[[i]] <- HRs_esc_rand_ns_2[[i]][1,]
}

for (i in 1:length(HRs_nor_2)){
  surv_nor_HRs_2[[i]] <- HRs_nor_2[[i]][1,]
  surv_nor_HRs_ns2[[i]] <- HRs_nor_ns_2[[i]][1,]
  surv_nor_HRs_ran2[[i]] <- HRs_nor_rand_2[[i]][1,]
  surv_nor_HRs_ns_ran2[[i]] <- HRs_nor_rand_ns_2[[i]][1,]
}

names(surv_both_HRs_2) <- PRSs
names(surv_esc_HRs_2) <- PRSs
names(surv_nor_HRs_2) <- PRSs

names(surv_both_HRs_ns2) <- PRSs
names(surv_esc_HRs_ns2) <- PRSs
names(surv_nor_HRs_ns2) <- PRSs

names(surv_both_HRs_ran2) <- PRSs
names(surv_esc_HRs_ran2) <- PRSs
names(surv_nor_HRs_ran2) <- PRSs

names(surv_both_HRs_ns_ran2) <- PRSs
names(surv_esc_HRs_ns_ran2) <- PRSs
names(surv_nor_HRs_ns_ran2) <- PRSs



## Make survival output file #####

surv_beta_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_2))[,1],
    t(as.data.frame(surv_both_betas_2))[,3],
    t(as.data.frame(surv_both_betas_2))[,6],
    t(as.data.frame(surv_esc_betas_2))[,1],
    t(as.data.frame(surv_esc_betas_2))[,3],
    t(as.data.frame(surv_esc_betas_2))[,6],
    t(as.data.frame(surv_nor_betas_2))[,1],
    t(as.data.frame(surv_nor_betas_2))[,3],
    t(as.data.frame(surv_nor_betas_2))[,6]
  )
)

colnames( surv_beta_matrix_2 ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_2))[,1],
    t(as.data.frame(surv_both_HRs_2))[,3],
    t(as.data.frame(surv_both_HRs_2))[,4],
    t(as.data.frame(surv_esc_HRs_2))[,1],
    t(as.data.frame(surv_esc_HRs_2))[,3],
    t(as.data.frame(surv_esc_HRs_2))[,4],
    t(as.data.frame(surv_nor_HRs_2))[,1],
    t(as.data.frame(surv_nor_HRs_2))[,3],
    t(as.data.frame(surv_nor_HRs_2))[,4]
  )
)

colnames( surv_HR_matrix_2 ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)


sink( "survival_betasfirstdrug.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_2, col.names=TRUE, sep="," )

cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_2, col.names=TRUE, sep="," )

cat("\nSurvival_expvar\n,")
write.table( rsquared_2,
             col.names=TRUE, sep="," )

sink()

## Make survival output file, non-switchers only #####

surv_beta_matrix_ns_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_ns2))[,1],
    t(as.data.frame(surv_both_betas_ns2))[,3],
    t(as.data.frame(surv_both_betas_ns2))[,6],
    t(as.data.frame(surv_esc_betas_ns2))[,1],
    t(as.data.frame(surv_esc_betas_ns2))[,3],
    t(as.data.frame(surv_esc_betas_ns2))[,6],
    t(as.data.frame(surv_nor_betas_ns2))[,1],
    t(as.data.frame(surv_nor_betas_ns2))[,3],
    t(as.data.frame(surv_nor_betas_ns2))[,6]
  )
)

colnames( surv_beta_matrix_ns_2 ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_ns_2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_ns2))[,1],
    t(as.data.frame(surv_both_HRs_ns2))[,3],
    t(as.data.frame(surv_both_HRs_ns2))[,4],
    t(as.data.frame(surv_esc_HRs_ns2))[,1],
    t(as.data.frame(surv_esc_HRs_ns2))[,3],
    t(as.data.frame(surv_esc_HRs_ns2))[,4],
    t(as.data.frame(surv_nor_HRs_ns2))[,1],
    t(as.data.frame(surv_nor_HRs_ns2))[,3],
    t(as.data.frame(surv_nor_HRs_ns2))[,4]
  )
)

colnames( surv_HR_matrix_ns_2 ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)


sink( "survival_nonswitchers_betasfirstdrug.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_ns_2, col.names=TRUE, sep="," )
cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_ns_2, col.names=TRUE, sep="," )
cat("\nSurvival_expvar\n,")
write.table( rsquared_ns_2,
             col.names=TRUE, sep="," )

sink()

## Make survival output file, random only #####

surv_beta_matrix_ran <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_ran2))[,1],
    t(as.data.frame(surv_both_betas_ran2))[,3],
    t(as.data.frame(surv_both_betas_ran2))[,6],
    t(as.data.frame(surv_esc_betas_ran2))[,1],
    t(as.data.frame(surv_esc_betas_ran2))[,3],
    t(as.data.frame(surv_esc_betas_ran2))[,6],
    t(as.data.frame(surv_nor_betas_ran2))[,1],
    t(as.data.frame(surv_nor_betas_ran2))[,3],
    t(as.data.frame(surv_nor_betas_ran2))[,6]
  )
)

colnames( surv_beta_matrix_ran ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_ran <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_ran2))[,1],
    t(as.data.frame(surv_both_HRs_ran2))[,3],
    t(as.data.frame(surv_both_HRs_ran2))[,4],
    t(as.data.frame(surv_esc_HRs_ran2))[,1],
    t(as.data.frame(surv_esc_HRs_ran2))[,3],
    t(as.data.frame(surv_esc_HRs_ran2))[,4],
    t(as.data.frame(surv_nor_HRs_ran2))[,1],
    t(as.data.frame(surv_nor_HRs_ran2))[,3],
    t(as.data.frame(surv_nor_HRs_ran2))[,4]
  )
)

colnames( surv_HR_matrix_ran ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)

sink( "survival_randomized_betasfirstdrug.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_ran, col.names=TRUE, sep="," )
cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_ran, col.names=TRUE, sep="," )
cat("\nSurvival_expvar\n,")
write.table( rsquared_2_r,
             col.names=TRUE, sep="," )

sink()

## Make survival output file, randomized non-switchers only #####

surv_beta_matrix_ns_ran2 <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_betas_ns_ran2))[,1],
    t(as.data.frame(surv_both_betas_ns_ran2))[,3],
    t(as.data.frame(surv_both_betas_ns_ran2))[,6],
    t(as.data.frame(surv_esc_betas_ns_ran2))[,1],
    t(as.data.frame(surv_esc_betas_ns_ran2))[,3],
    t(as.data.frame(surv_esc_betas_ns_ran2))[,6],
    t(as.data.frame(surv_nor_betas_ns_ran2))[,1],
    t(as.data.frame(surv_nor_betas_ns_ran2))[,3],
    t(as.data.frame(surv_nor_betas_ns_ran2))[,6]
  )
)

colnames( surv_beta_matrix_ns_ran2 ) <- c(
  "surv_both_B", "surv_both_se", "surv_both_p",
  "surv_esc_B", "surv_esc_se", "surv_esc_p",
  "surv_nor_B", "surv_nor_se", "surv_nor_p"
)

surv_HR_matrix_ns_ran <- as.matrix(
  cbind(
    t(as.data.frame(surv_both_HRs_ns_ran2))[,1],
    t(as.data.frame(surv_both_HRs_ns_ran2))[,3],
    t(as.data.frame(surv_both_HRs_ns_ran2))[,4],
    t(as.data.frame(surv_esc_HRs_ns_ran2))[,1],
    t(as.data.frame(surv_esc_HRs_ns_ran2))[,3],
    t(as.data.frame(surv_esc_HRs_ns_ran2))[,4],
    t(as.data.frame(surv_nor_HRs_ns_ran2))[,1],
    t(as.data.frame(surv_nor_HRs_ns_ran2))[,3],
    t(as.data.frame(surv_nor_HRs_ns_ran2))[,4]
  )
)

colnames( surv_HR_matrix_ns_ran ) <- c(
  "surv_both_HR", "surv_both_LCI", "surv_both_UCI",
  "surv_esc_HR", "surv_esc_LCI", "surv_esc_UCI",
  "surv_nor_HR", "surv_nor_LCI", "surv_nor_UCI"
)

sink( "survival_rand_nonswitchers_betasfirstdrug.csv")

cat("Survival_betas\n,")
write.table( surv_beta_matrix_ns_ran2, col.names=TRUE, sep="," )

cat("\nSurvival_Hazard_Ratios\n,")
write.table( surv_HR_matrix_ns_ran, col.names=TRUE, sep="," )

cat("\nSurvival_expvar\n,")
write.table( rsquared_ns_r_2,
             col.names=TRUE, sep="," )

sink()

## Make survival-only heatmap #####

surv_models_betas_2 <- c("surv_both_s_B", "surv_esc_s_B", "surv_nor_s_B")
surv_models_ses_2 <- c("surv_both_s_se", "surv_esc_s_se", "surv_nor_s_se")
surv_models_ps_2 <- c("surv_both_s_p", "surv_esc_s_p", "surv_nor_s_p")

surv_models_betas_columns_2 <- which( colnames(surv_beta_matrix_2)
                                      %in% surv_models_betas_2 )
surv_models_ses_columns_2 <- which( colnames(surv_beta_matrix_2)
                                    %in% surv_models_ses_2 )
surv_models_ps_columns_2 <- which( colnames(surv_beta_matrix_2)
                                   %in% surv_models_ps_2 )

surv_betas_only_2 <- cbind(surv_beta_matrix_2[,surv_models_betas_columns_2[1]],
                           surv_beta_matrix_2[,surv_models_betas_columns_2[2]],
                           surv_beta_matrix_2[,surv_models_betas_columns_2[3]])

colnames(surv_betas_only_2) <- c("surv_both", "surv_esc","surv_nor")


surv_betas_only_2 <- cbind(surv_beta_matrix_2[,1],
                           surv_beta_matrix_2[,4],
                           surv_beta_matrix_2[,7])

colnames(surv_betas_only_2) <- c("surv_both", "surv_esc", "surv_nor")

surv_se_only_2 <- cbind(surv_beta_matrix_2[,2],
                        surv_beta_matrix_2[,5],
                        surv_beta_matrix_2[,8])

colnames(surv_se_only_2) <- c("surv_both", "surv_esc","surv_nor")

surv_p_only_2 <- cbind(surv_beta_matrix_2[,3],
                       surv_beta_matrix_2[,6],
                       surv_beta_matrix_2[,9])

colnames(surv_p_only_2) <- c("surv_both", "surv_esc","surv_nor")

lowest_ps <- matrix(nrow=30,ncol=3,dimnames=list(PRSs,c( "both",
                                                         "esc",
                                                         "nor")))

lowest_ps[which(surv_p_only_2<=0.05)]<-surv_p_only_2[which(surv_p_only_2<=0.05)]

heatmap.2( surv_betas_only_2, col=redblue,
           Rowv=FALSE, Colv=FALSE,
           density.info="none", trace="none",cexCol=0.7,
           cellnote=round(lowest_ps,5),notecol="black",
           main="Dropout of first drug"
)

##### IMPORT p-values TABLE #####

# I manually pasted together all of the results and saved it as a csv

pvalues <- read.csv( "p_values_for_fdr.csv" )

# Remove my "sort" column that I added for use in Excel

pvalues <- pvalues[,2:8]

hist( pvalues$pvalue )

pvalues_r_nr <- subset( pvalues, Ran_Everyone == "Everyone" )
pvalues_r <- subset( pvalues, Ran_Everyone == "Ran" )

# PLOT BH line for all #####

pdf( "False_discovery_rate.pdf" )

pvalues <- pvalues[ order( pvalues$pvalue ), ]

pvalues$Bonferroni =
  p.adjust(pvalues$pvalue,
           method = "bonferroni")

pvalues$BH =
  p.adjust(pvalues$pvalue,
           method = "BH")

pvalues$Holm =
  p.adjust(pvalues$pvalue,
           method = "holm")

pvalues$Hochberg =
  p.adjust(pvalues$pvalue,
           method = "hochberg")

pvalues$Hommel =
  p.adjust(pvalues$pvalue,
           method = "hommel")

pvalues$BY =
  p.adjust(pvalues$pvalue,
           method = "BY")

X = pvalues$pvalue
Y = cbind(pvalues$Bonferroni,
          pvalues$BH,
          pvalues$Holm,
          pvalues$Hochberg,
          pvalues$Hommel,
          pvalues$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on all results simultaneously")

# Make tables for each drug-polygenic-score combination #####

anor <- subset( pvalues, PS %in% c( "AZ_0_0001","AZ_0_01","AZ_0_05","AZ_0_1","AZ_0_5","AZ_1" ) )
react <- subset( pvalues, PS %in% c( "CZ_0_0001","CZ_0_01","CZ_0_05","CZ_0_1","CZ_0_5","CZ_1" ) )
neur <- subset( pvalues, PS %in% c( "NZ_0_0001","NZ_0_01","NZ_0_05","NZ_0_1","NZ_0_5","NZ_1" ) )
schi <- subset( pvalues, PS %in% c( "SZ_0_0001","SZ_0_01","SZ_0_05","SZ_0_1","SZ_0_5","SZ_1" ) )
hei <- subset( pvalues, PS %in% c( "HZ_0_0001","HZ_0_01","HZ_0_05","HZ_0_1","HZ_0_5","HZ_1" ) )

esc_anor <- subset( anor, Drug=="Escitalopram" )
esc_react <- subset( react, Drug=="Escitalopram" )
esc_neur <- subset( neur, Drug=="Escitalopram" )
esc_schi <- subset( schi, Drug=="Escitalopram" )
esc_hei <- subset( hei, Drug=="Escitalopram" )

nor_anor <- subset( anor, Drug=="Nortriptyline" )
nor_react <- subset( react, Drug=="Nortriptyline" )
nor_neur <- subset( neur, Drug=="Nortriptyline" )
nor_schi <- subset( schi, Drug=="Nortriptyline" )
nor_hei <- subset( hei, Drug=="Nortriptyline" )

# esc_anor #####

esc_anor <- esc_anor[ order( esc_anor$pvalue ), ]

esc_anor$Bonferroni =
  p.adjust(esc_anor$pvalue,
           method = "bonferroni")

esc_anor$BH =
  p.adjust(esc_anor$pvalue,
           method = "BH")

esc_anor$Holm =
  p.adjust(esc_anor$pvalue,
           method = "holm")

esc_anor$Hochberg =
  p.adjust(esc_anor$pvalue,
           method = "hochberg")

esc_anor$Hommel =
  p.adjust(esc_anor$pvalue,
           method = "hommel")

esc_anor$BY =
  p.adjust(esc_anor$pvalue,
           method = "BY")

X = esc_anor$pvalue
Y = cbind(esc_anor$Bonferroni,
          esc_anor$BH,
          esc_anor$Holm,
          esc_anor$Hochberg,
          esc_anor$Hommel,
          esc_anor$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-anorexia-PS results only")

# esc_react #####

esc_react <- esc_react[ order( esc_react$pvalue ), ]

esc_react$Bonferroni =
  p.adjust(esc_react$pvalue,
           method = "bonferroni")

esc_react$BH =
  p.adjust(esc_react$pvalue,
           method = "BH")

esc_react$Holm =
  p.adjust(esc_react$pvalue,
           method = "holm")

esc_react$Hochberg =
  p.adjust(esc_react$pvalue,
           method = "hochberg")

esc_react$Hommel =
  p.adjust(esc_react$pvalue,
           method = "hommel")

esc_react$BY =
  p.adjust(esc_react$pvalue,
           method = "BY")

X = esc_react$pvalue
Y = cbind(esc_react$Bonferroni,
          esc_react$BH,
          esc_react$Holm,
          esc_react$Hochberg,
          esc_react$Hommel,
          esc_react$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-reaction-time-PS results only")

# esc_neur #####

esc_neur <- esc_neur[ order( esc_neur$pvalue ), ]

esc_neur$Bonferroni =
  p.adjust(esc_neur$pvalue,
           method = "bonferroni")

esc_neur$BH =
  p.adjust(esc_neur$pvalue,
           method = "BH")

esc_neur$Holm =
  p.adjust(esc_neur$pvalue,
           method = "holm")

esc_neur$Hochberg =
  p.adjust(esc_neur$pvalue,
           method = "hochberg")

esc_neur$Hommel =
  p.adjust(esc_neur$pvalue,
           method = "hommel")

esc_neur$BY =
  p.adjust(esc_neur$pvalue,
           method = "BY")

X = esc_neur$pvalue
Y = cbind(esc_neur$Bonferroni,
          esc_neur$BH,
          esc_neur$Holm,
          esc_neur$Hochberg,
          esc_neur$Hommel,
          esc_neur$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-neuroticism-PS results only")

# esc_schi #####

esc_schi <- esc_schi[ order( esc_schi$pvalue ), ]

esc_schi$Bonferroni =
  p.adjust(esc_schi$pvalue,
           method = "bonferroni")

esc_schi$BH =
  p.adjust(esc_schi$pvalue,
           method = "BH")

esc_schi$Holm =
  p.adjust(esc_schi$pvalue,
           method = "holm")

esc_schi$Hochberg =
  p.adjust(esc_schi$pvalue,
           method = "hochberg")

esc_schi$Hommel =
  p.adjust(esc_schi$pvalue,
           method = "hommel")

esc_schi$BY =
  p.adjust(esc_schi$pvalue,
           method = "BY")

X = esc_schi$pvalue
Y = cbind(esc_schi$Bonferroni,
          esc_schi$BH,
          esc_schi$Holm,
          esc_schi$Hochberg,
          esc_schi$Hommel,
          esc_schi$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-schizophrenia-PS results only")

# esc_hei #####

esc_hei <- esc_hei[ order( esc_hei$pvalue ), ]

esc_hei$Bonferroni =
  p.adjust(esc_hei$pvalue,
           method = "bonferroni")

esc_hei$BH =
  p.adjust(esc_hei$pvalue,
           method = "BH")

esc_hei$Holm =
  p.adjust(esc_hei$pvalue,
           method = "holm")

esc_hei$Hochberg =
  p.adjust(esc_hei$pvalue,
           method = "hochberg")

esc_hei$Hommel =
  p.adjust(esc_hei$pvalue,
           method = "hommel")

esc_hei$BY =
  p.adjust(esc_hei$pvalue,
           method = "BY")

X = esc_hei$pvalue
Y = cbind(esc_hei$Bonferroni,
          esc_hei$BH,
          esc_hei$Holm,
          esc_hei$Hochberg,
          esc_hei$Hommel,
          esc_hei$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-height-PS results only")

# nor_anor #####

nor_anor <- nor_anor[ order( nor_anor$pvalue ), ]

nor_anor$Bonferroni =
  p.adjust(nor_anor$pvalue,
           method = "bonferroni")

nor_anor$BH =
  p.adjust(nor_anor$pvalue,
           method = "BH")

nor_anor$Holm =
  p.adjust(nor_anor$pvalue,
           method = "holm")

nor_anor$Hochberg =
  p.adjust(nor_anor$pvalue,
           method = "hochberg")

nor_anor$Hommel =
  p.adjust(nor_anor$pvalue,
           method = "hommel")

nor_anor$BY =
  p.adjust(nor_anor$pvalue,
           method = "BY")

X = nor_anor$pvalue
Y = cbind(nor_anor$Bonferroni,
          nor_anor$BH,
          nor_anor$Holm,
          nor_anor$Hochberg,
          nor_anor$Hommel,
          nor_anor$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-anorexia-PS results only")

# nor_react #####

nor_react <- nor_react[ order( nor_react$pvalue ), ]

nor_react$Bonferroni =
  p.adjust(nor_react$pvalue,
           method = "bonferroni")

nor_react$BH =
  p.adjust(nor_react$pvalue,
           method = "BH")

nor_react$Holm =
  p.adjust(nor_react$pvalue,
           method = "holm")

nor_react$Hochberg =
  p.adjust(nor_react$pvalue,
           method = "hochberg")

nor_react$Hommel =
  p.adjust(nor_react$pvalue,
           method = "hommel")

nor_react$BY =
  p.adjust(nor_react$pvalue,
           method = "BY")

X = nor_react$pvalue
Y = cbind(nor_react$Bonferroni,
          nor_react$BH,
          nor_react$Holm,
          nor_react$Hochberg,
          nor_react$Hommel,
          nor_react$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-reaction-time-PS results only")

# nor_neur #####

nor_neur <- nor_neur[ order( nor_neur$pvalue ), ]

nor_neur$Bonferroni =
  p.adjust(nor_neur$pvalue,
           method = "bonferroni")

nor_neur$BH =
  p.adjust(nor_neur$pvalue,
           method = "BH")

nor_neur$Holm =
  p.adjust(nor_neur$pvalue,
           method = "holm")

nor_neur$Hochberg =
  p.adjust(nor_neur$pvalue,
           method = "hochberg")

nor_neur$Hommel =
  p.adjust(nor_neur$pvalue,
           method = "hommel")

nor_neur$BY =
  p.adjust(nor_neur$pvalue,
           method = "BY")

X = nor_neur$pvalue
Y = cbind(nor_neur$Bonferroni,
          nor_neur$BH,
          nor_neur$Holm,
          nor_neur$Hochberg,
          nor_neur$Hommel,
          nor_neur$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-neuroticism-PS results only")

# nor_schi #####

nor_schi <- nor_schi[ order( nor_schi$pvalue ), ]

nor_schi$Bonferroni =
  p.adjust(nor_schi$pvalue,
           method = "bonferroni")

nor_schi$BH =
  p.adjust(nor_schi$pvalue,
           method = "BH")

nor_schi$Holm =
  p.adjust(nor_schi$pvalue,
           method = "holm")

nor_schi$Hochberg =
  p.adjust(nor_schi$pvalue,
           method = "hochberg")

nor_schi$Hommel =
  p.adjust(nor_schi$pvalue,
           method = "hommel")

nor_schi$BY =
  p.adjust(nor_schi$pvalue,
           method = "BY")

X = nor_schi$pvalue
Y = cbind(nor_schi$Bonferroni,
          nor_schi$BH,
          nor_schi$Holm,
          nor_schi$Hochberg,
          nor_schi$Hommel,
          nor_schi$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-schizophrenia-PS results only")

# nor_hei #####

nor_hei <- nor_hei[ order( nor_hei$pvalue ), ]

nor_hei$Bonferroni =
  p.adjust(nor_hei$pvalue,
           method = "bonferroni")

nor_hei$BH =
  p.adjust(nor_hei$pvalue,
           method = "BH")

nor_hei$Holm =
  p.adjust(nor_hei$pvalue,
           method = "holm")

nor_hei$Hochberg =
  p.adjust(nor_hei$pvalue,
           method = "hochberg")

nor_hei$Hommel =
  p.adjust(nor_hei$pvalue,
           method = "hommel")

nor_hei$BY =
  p.adjust(nor_hei$pvalue,
           method = "BY")

X = nor_hei$pvalue
Y = cbind(nor_hei$Bonferroni,
          nor_hei$BH,
          nor_hei$Holm,
          nor_hei$Hochberg,
          nor_hei$Hommel,
          nor_hei$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2
        ,
        xlim=c(0,1),
        ylim=c(0,1)
        )

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-height-PS results only")

dev.off()

# Getting top BH results #####

range( pvalues$BH ) # 0.866-1.000

range( esc_anor$BH ) # 0.242-0.997
range( esc_react$BH ) # 0.544-0.992
range( esc_neur$BH ) # 0.850-0.999
range( esc_schi$BH ) # 0.954-0.999
range( esc_hei$BH ) # 0.640-1.000

range( nor_anor$BH ) # 0.309-0.995
range( nor_react$BH ) # 0.735-0.999
range( nor_neur$BH ) # 0.915-1.000
range( nor_schi$BH ) # 0.194-0.987
range( nor_hei$BH ) # 0.167-0.999

lowest_pvalues <- subset( pvalues, pvalue<0.05 ) # 320

lowest_esc_anor <- subset( esc_anor, pvalue<0.05 ) # 10
# 4x light-headedness, 3x dose, 1x tremor, inc app, dec app
lowest_esc_react <- subset( esc_react, pvalue<0.05 ) # 30
# 5x cognitive, 4x tremor, 3x totasec, zmadrs, 2x dry mouth,
# drowsiness, insomnia, headache, dec app, mood, study dropout,
# 1x sweating
lowest_esc_neur <- subset( esc_neur, pvalue<0.05 ) # 20
# 6x sweating, 4x inc app, 3x disorientation, 2x urination problems,
# blurred vision, 1x sex problems, weight gain, insomnia
lowest_esc_schi <- subset( esc_schi, pvalue<0.05 ) # 12
# 10x study dropout, 1x dose, sweating
lowest_esc_hei <- subset( esc_hei, pvalue<0.05 ) # 17
# 5x cognitive, 4x dec app, 3x sweating, 2x dry mouth, naus/vom,
# 1x yawning

lowest_nor_anor <- subset( nor_anor, pvalue<0.05 ) # 53
lowest_nor_react <- subset( nor_react, pvalue<0.05 ) # 21
lowest_nor_neur <- subset( nor_neur, pvalue<0.05 ) # 11
lowest_nor_schi <- subset( nor_schi, pvalue<0.05 ) # 18
lowest_nor_hei <- subset( nor_hei, pvalue<0.05 ) # 26



# PLOT BH line for all - ran and non-ran #####

pdf( "False_discovery_rate_ran_nonran.pdf" )


pvalues_r_nr <- pvalues_r_nr[ order( pvalues_r_nr$pvalue ), ]

pvalues_r_nr$Bonferroni =
  p.adjust(pvalues_r_nr$pvalue,
           method = "bonferroni")

pvalues_r_nr$BH =
  p.adjust(pvalues_r_nr$pvalue,
           method = "BH")

pvalues_r_nr$Holm =
  p.adjust(pvalues_r_nr$pvalue,
           method = "holm")

pvalues_r_nr$Hochberg =
  p.adjust(pvalues_r_nr$pvalue,
           method = "hochberg")

pvalues_r_nr$Hommel =
  p.adjust(pvalues_r_nr$pvalue,
           method = "hommel")

pvalues_r_nr$BY =
  p.adjust(pvalues_r_nr$pvalue,
           method = "BY")

X = pvalues_r_nr$pvalue
Y = cbind(pvalues_r_nr$Bonferroni,
          pvalues_r_nr$BH,
          pvalues_r_nr$Holm,
          pvalues_r_nr$Hochberg,
          pvalues_r_nr$Hommel,
          pvalues_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on all results simultaneously, ran and non-ran")

# Make tables for each drug-polygenic-score combination #####

anor_r_nr <- subset( pvalues_r_nr, PS %in% c( "AZ_0_0001","AZ_0_01","AZ_0_05","AZ_0_1","AZ_0_5","AZ_1" ) )
react_r_nr <- subset( pvalues_r_nr, PS %in% c( "CZ_0_0001","CZ_0_01","CZ_0_05","CZ_0_1","CZ_0_5","CZ_1" ) )
neur_r_nr <- subset( pvalues_r_nr, PS %in% c( "NZ_0_0001","NZ_0_01","NZ_0_05","NZ_0_1","NZ_0_5","NZ_1" ) )
schi_r_nr <- subset( pvalues_r_nr, PS %in% c( "SZ_0_0001","SZ_0_01","SZ_0_05","SZ_0_1","SZ_0_5","SZ_1" ) )
hei_r_nr <- subset( pvalues_r_nr, PS %in% c( "HZ_0_0001","HZ_0_01","HZ_0_05","HZ_0_1","HZ_0_5","HZ_1" ) )

esc_anor_r_nr <- subset( anor_r_nr, Drug=="Escitalopram" )
esc_react_r_nr <- subset( react_r_nr, Drug=="Escitalopram" )
esc_neur_r_nr <- subset( neur_r_nr, Drug=="Escitalopram" )
esc_schi_r_nr <- subset( schi_r_nr, Drug=="Escitalopram" )
esc_hei_r_nr <- subset( hei_r_nr, Drug=="Escitalopram" )

nor_anor_r_nr <- subset( anor_r_nr, Drug=="Nortriptyline" )
nor_react_r_nr <- subset( react_r_nr, Drug=="Nortriptyline" )
nor_neur_r_nr <- subset( neur_r_nr, Drug=="Nortriptyline" )
nor_schi_r_nr <- subset( schi_r_nr, Drug=="Nortriptyline" )
nor_hei_r_nr <- subset( hei_r_nr, Drug=="Nortriptyline" )

# esc_anor_r_nr #####

esc_anor_r_nr <- esc_anor_r_nr[ order( esc_anor_r_nr$pvalue ), ]

esc_anor_r_nr$Bonferroni =
  p.adjust(esc_anor_r_nr$pvalue,
           method = "bonferroni")

esc_anor_r_nr$BH =
  p.adjust(esc_anor_r_nr$pvalue,
           method = "BH")

esc_anor_r_nr$Holm =
  p.adjust(esc_anor_r_nr$pvalue,
           method = "holm")

esc_anor_r_nr$Hochberg =
  p.adjust(esc_anor_r_nr$pvalue,
           method = "hochberg")

esc_anor_r_nr$Hommel =
  p.adjust(esc_anor_r_nr$pvalue,
           method = "hommel")

esc_anor_r_nr$BY =
  p.adjust(esc_anor_r_nr$pvalue,
           method = "BY")

X = esc_anor_r_nr$pvalue
Y = cbind(esc_anor_r_nr$Bonferroni,
          esc_anor_r_nr$BH,
          esc_anor_r_nr$Holm,
          esc_anor_r_nr$Hochberg,
          esc_anor_r_nr$Hommel,
          esc_anor_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-anorexia-PS results only, ran and non-ran")

# esc_react_r_nr #####

esc_react_r_nr <- esc_react_r_nr[ order( esc_react_r_nr$pvalue ), ]

esc_react_r_nr$Bonferroni =
  p.adjust(esc_react_r_nr$pvalue,
           method = "bonferroni")

esc_react_r_nr$BH =
  p.adjust(esc_react_r_nr$pvalue,
           method = "BH")

esc_react_r_nr$Holm =
  p.adjust(esc_react_r_nr$pvalue,
           method = "holm")

esc_react_r_nr$Hochberg =
  p.adjust(esc_react_r_nr$pvalue,
           method = "hochberg")

esc_react_r_nr$Hommel =
  p.adjust(esc_react_r_nr$pvalue,
           method = "hommel")

esc_react_r_nr$BY =
  p.adjust(esc_react_r_nr$pvalue,
           method = "BY")

X = esc_react_r_nr$pvalue
Y = cbind(esc_react_r_nr$Bonferroni,
          esc_react_r_nr$BH,
          esc_react_r_nr$Holm,
          esc_react_r_nr$Hochberg,
          esc_react_r_nr$Hommel,
          esc_react_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-reaction-time-PS results only, ran and non-ran")

# esc_neur_r_nr #####

esc_neur_r_nr <- esc_neur_r_nr[ order( esc_neur_r_nr$pvalue ), ]

esc_neur_r_nr$Bonferroni =
  p.adjust(esc_neur_r_nr$pvalue,
           method = "bonferroni")

esc_neur_r_nr$BH =
  p.adjust(esc_neur_r_nr$pvalue,
           method = "BH")

esc_neur_r_nr$Holm =
  p.adjust(esc_neur_r_nr$pvalue,
           method = "holm")

esc_neur_r_nr$Hochberg =
  p.adjust(esc_neur_r_nr$pvalue,
           method = "hochberg")

esc_neur_r_nr$Hommel =
  p.adjust(esc_neur_r_nr$pvalue,
           method = "hommel")

esc_neur_r_nr$BY =
  p.adjust(esc_neur_r_nr$pvalue,
           method = "BY")

X = esc_neur_r_nr$pvalue
Y = cbind(esc_neur_r_nr$Bonferroni,
          esc_neur_r_nr$BH,
          esc_neur_r_nr$Holm,
          esc_neur_r_nr$Hochberg,
          esc_neur_r_nr$Hommel,
          esc_neur_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-neuroticism-PS results only, ran and non-ran")

# esc_schi_r_nr #####

esc_schi_r_nr <- esc_schi_r_nr[ order( esc_schi_r_nr$pvalue ), ]

esc_schi_r_nr$Bonferroni =
  p.adjust(esc_schi_r_nr$pvalue,
           method = "bonferroni")

esc_schi_r_nr$BH =
  p.adjust(esc_schi_r_nr$pvalue,
           method = "BH")

esc_schi_r_nr$Holm =
  p.adjust(esc_schi_r_nr$pvalue,
           method = "holm")

esc_schi_r_nr$Hochberg =
  p.adjust(esc_schi_r_nr$pvalue,
           method = "hochberg")

esc_schi_r_nr$Hommel =
  p.adjust(esc_schi_r_nr$pvalue,
           method = "hommel")

esc_schi_r_nr$BY =
  p.adjust(esc_schi_r_nr$pvalue,
           method = "BY")

X = esc_schi_r_nr$pvalue
Y = cbind(esc_schi_r_nr$Bonferroni,
          esc_schi_r_nr$BH,
          esc_schi_r_nr$Holm,
          esc_schi_r_nr$Hochberg,
          esc_schi_r_nr$Hommel,
          esc_schi_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-schizophrenia-PS results only, ran and non-ran")

# esc_hei_r_nr #####

esc_hei_r_nr <- esc_hei_r_nr[ order( esc_hei_r_nr$pvalue ), ]

esc_hei_r_nr$Bonferroni =
  p.adjust(esc_hei_r_nr$pvalue,
           method = "bonferroni")

esc_hei_r_nr$BH =
  p.adjust(esc_hei_r_nr$pvalue,
           method = "BH")

esc_hei_r_nr$Holm =
  p.adjust(esc_hei_r_nr$pvalue,
           method = "holm")

esc_hei_r_nr$Hochberg =
  p.adjust(esc_hei_r_nr$pvalue,
           method = "hochberg")

esc_hei_r_nr$Hommel =
  p.adjust(esc_hei_r_nr$pvalue,
           method = "hommel")

esc_hei_r_nr$BY =
  p.adjust(esc_hei_r_nr$pvalue,
           method = "BY")

X = esc_hei_r_nr$pvalue
Y = cbind(esc_hei_r_nr$Bonferroni,
          esc_hei_r_nr$BH,
          esc_hei_r_nr$Holm,
          esc_hei_r_nr$Hochberg,
          esc_hei_r_nr$Hommel,
          esc_hei_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-height-PS results only, ran and non-ran")

# nor_anor_r_nr #####

nor_anor_r_nr <- nor_anor_r_nr[ order( nor_anor_r_nr$pvalue ), ]

nor_anor_r_nr$Bonferroni =
  p.adjust(nor_anor_r_nr$pvalue,
           method = "bonferroni")

nor_anor_r_nr$BH =
  p.adjust(nor_anor_r_nr$pvalue,
           method = "BH")

nor_anor_r_nr$Holm =
  p.adjust(nor_anor_r_nr$pvalue,
           method = "holm")

nor_anor_r_nr$Hochberg =
  p.adjust(nor_anor_r_nr$pvalue,
           method = "hochberg")

nor_anor_r_nr$Hommel =
  p.adjust(nor_anor_r_nr$pvalue,
           method = "hommel")

nor_anor_r_nr$BY =
  p.adjust(nor_anor_r_nr$pvalue,
           method = "BY")

X = nor_anor_r_nr$pvalue
Y = cbind(nor_anor_r_nr$Bonferroni,
          nor_anor_r_nr$BH,
          nor_anor_r_nr$Holm,
          nor_anor_r_nr$Hochberg,
          nor_anor_r_nr$Hommel,
          nor_anor_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-anorexia-PS results only, ran and non-ran")

# nor_react_r_nr #####

nor_react_r_nr <- nor_react_r_nr[ order( nor_react_r_nr$pvalue ), ]

nor_react_r_nr$Bonferroni =
  p.adjust(nor_react_r_nr$pvalue,
           method = "bonferroni")

nor_react_r_nr$BH =
  p.adjust(nor_react_r_nr$pvalue,
           method = "BH")

nor_react_r_nr$Holm =
  p.adjust(nor_react_r_nr$pvalue,
           method = "holm")

nor_react_r_nr$Hochberg =
  p.adjust(nor_react_r_nr$pvalue,
           method = "hochberg")

nor_react_r_nr$Hommel =
  p.adjust(nor_react_r_nr$pvalue,
           method = "hommel")

nor_react_r_nr$BY =
  p.adjust(nor_react_r_nr$pvalue,
           method = "BY")

X = nor_react_r_nr$pvalue
Y = cbind(nor_react_r_nr$Bonferroni,
          nor_react_r_nr$BH,
          nor_react_r_nr$Holm,
          nor_react_r_nr$Hochberg,
          nor_react_r_nr$Hommel,
          nor_react_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-reaction-time-PS results only, ran and non-ran")

# nor_neur_r_nr #####

nor_neur_r_nr <- nor_neur_r_nr[ order( nor_neur_r_nr$pvalue ), ]

nor_neur_r_nr$Bonferroni =
  p.adjust(nor_neur_r_nr$pvalue,
           method = "bonferroni")

nor_neur_r_nr$BH =
  p.adjust(nor_neur_r_nr$pvalue,
           method = "BH")

nor_neur_r_nr$Holm =
  p.adjust(nor_neur_r_nr$pvalue,
           method = "holm")

nor_neur_r_nr$Hochberg =
  p.adjust(nor_neur_r_nr$pvalue,
           method = "hochberg")

nor_neur_r_nr$Hommel =
  p.adjust(nor_neur_r_nr$pvalue,
           method = "hommel")

nor_neur_r_nr$BY =
  p.adjust(nor_neur_r_nr$pvalue,
           method = "BY")

X = nor_neur_r_nr$pvalue
Y = cbind(nor_neur_r_nr$Bonferroni,
          nor_neur_r_nr$BH,
          nor_neur_r_nr$Holm,
          nor_neur_r_nr$Hochberg,
          nor_neur_r_nr$Hommel,
          nor_neur_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-neuroticism-PS results only, ran and non-ran")

# nor_schi_r_nr #####

nor_schi_r_nr <- nor_schi_r_nr[ order( nor_schi_r_nr$pvalue ), ]

nor_schi_r_nr$Bonferroni =
  p.adjust(nor_schi_r_nr$pvalue,
           method = "bonferroni")

nor_schi_r_nr$BH =
  p.adjust(nor_schi_r_nr$pvalue,
           method = "BH")

nor_schi_r_nr$Holm =
  p.adjust(nor_schi_r_nr$pvalue,
           method = "holm")

nor_schi_r_nr$Hochberg =
  p.adjust(nor_schi_r_nr$pvalue,
           method = "hochberg")

nor_schi_r_nr$Hommel =
  p.adjust(nor_schi_r_nr$pvalue,
           method = "hommel")

nor_schi_r_nr$BY =
  p.adjust(nor_schi_r_nr$pvalue,
           method = "BY")

X = nor_schi_r_nr$pvalue
Y = cbind(nor_schi_r_nr$Bonferroni,
          nor_schi_r_nr$BH,
          nor_schi_r_nr$Holm,
          nor_schi_r_nr$Hochberg,
          nor_schi_r_nr$Hommel,
          nor_schi_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-schizophrenia-PS results only, ran and non-ran")

# nor_hei_r_nr #####

nor_hei_r_nr <- nor_hei_r_nr[ order( nor_hei_r_nr$pvalue ), ]

nor_hei_r_nr$Bonferroni =
  p.adjust(nor_hei_r_nr$pvalue,
           method = "bonferroni")

nor_hei_r_nr$BH =
  p.adjust(nor_hei_r_nr$pvalue,
           method = "BH")

nor_hei_r_nr$Holm =
  p.adjust(nor_hei_r_nr$pvalue,
           method = "holm")

nor_hei_r_nr$Hochberg =
  p.adjust(nor_hei_r_nr$pvalue,
           method = "hochberg")

nor_hei_r_nr$Hommel =
  p.adjust(nor_hei_r_nr$pvalue,
           method = "hommel")

nor_hei_r_nr$BY =
  p.adjust(nor_hei_r_nr$pvalue,
           method = "BY")

X = nor_hei_r_nr$pvalue
Y = cbind(nor_hei_r_nr$Bonferroni,
          nor_hei_r_nr$BH,
          nor_hei_r_nr$Holm,
          nor_hei_r_nr$Hochberg,
          nor_hei_r_nr$Hommel,
          nor_hei_r_nr$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2
        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-height-PS results only, ran and non-ran")

dev.off()

# PLOT BH line for all - ran #####

pdf( "False_discovery_rate_ran.pdf" )


pvalues_r <- pvalues_r[ order( pvalues_r$pvalue ), ]

pvalues_r$Bonferroni =
  p.adjust(pvalues_r$pvalue,
           method = "bonferroni")

pvalues_r$BH =
  p.adjust(pvalues_r$pvalue,
           method = "BH")

pvalues_r$Holm =
  p.adjust(pvalues_r$pvalue,
           method = "holm")

pvalues_r$Hochberg =
  p.adjust(pvalues_r$pvalue,
           method = "hochberg")

pvalues_r$Hommel =
  p.adjust(pvalues_r$pvalue,
           method = "hommel")

pvalues_r$BY =
  p.adjust(pvalues_r$pvalue,
           method = "BY")

X = pvalues_r$pvalue
Y = cbind(pvalues_r$Bonferroni,
          pvalues_r$BH,
          pvalues_r$Holm,
          pvalues_r$Hochberg,
          pvalues_r$Hommel,
          pvalues_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on all results simultaneously, ran")

# Make tables for each drug-polygenic-score combination #####

anor_r <- subset( pvalues_r, PS %in% c( "AZ_0_0001","AZ_0_01","AZ_0_05","AZ_0_1","AZ_0_5","AZ_1" ) )
react_r <- subset( pvalues_r, PS %in% c( "CZ_0_0001","CZ_0_01","CZ_0_05","CZ_0_1","CZ_0_5","CZ_1" ) )
neur_r <- subset( pvalues_r, PS %in% c( "NZ_0_0001","NZ_0_01","NZ_0_05","NZ_0_1","NZ_0_5","NZ_1" ) )
schi_r <- subset( pvalues_r, PS %in% c( "SZ_0_0001","SZ_0_01","SZ_0_05","SZ_0_1","SZ_0_5","SZ_1" ) )
hei_r <- subset( pvalues_r, PS %in% c( "HZ_0_0001","HZ_0_01","HZ_0_05","HZ_0_1","HZ_0_5","HZ_1" ) )

esc_anor_r <- subset( anor_r, Drug=="Escitalopram" )
esc_react_r <- subset( react_r, Drug=="Escitalopram" )
esc_neur_r <- subset( neur_r, Drug=="Escitalopram" )
esc_schi_r <- subset( schi_r, Drug=="Escitalopram" )
esc_hei_r <- subset( hei_r, Drug=="Escitalopram" )

nor_anor_r <- subset( anor_r, Drug=="Nortriptyline" )
nor_react_r <- subset( react_r, Drug=="Nortriptyline" )
nor_neur_r <- subset( neur_r, Drug=="Nortriptyline" )
nor_schi_r <- subset( schi_r, Drug=="Nortriptyline" )
nor_hei_r <- subset( hei_r, Drug=="Nortriptyline" )

# esc_anor_r #####

esc_anor_r <- esc_anor_r[ order( esc_anor_r$pvalue ), ]

esc_anor_r$Bonferroni =
  p.adjust(esc_anor_r$pvalue,
           method = "bonferroni")

esc_anor_r$BH =
  p.adjust(esc_anor_r$pvalue,
           method = "BH")

esc_anor_r$Holm =
  p.adjust(esc_anor_r$pvalue,
           method = "holm")

esc_anor_r$Hochberg =
  p.adjust(esc_anor_r$pvalue,
           method = "hochberg")

esc_anor_r$Hommel =
  p.adjust(esc_anor_r$pvalue,
           method = "hommel")

esc_anor_r$BY =
  p.adjust(esc_anor_r$pvalue,
           method = "BY")

X = esc_anor_r$pvalue
Y = cbind(esc_anor_r$Bonferroni,
          esc_anor_r$BH,
          esc_anor_r$Holm,
          esc_anor_r$Hochberg,
          esc_anor_r$Hommel,
          esc_anor_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-anorexia-PS results only, ran")

# esc_react_r #####

esc_react_r <- esc_react_r[ order( esc_react_r$pvalue ), ]

esc_react_r$Bonferroni =
  p.adjust(esc_react_r$pvalue,
           method = "bonferroni")

esc_react_r$BH =
  p.adjust(esc_react_r$pvalue,
           method = "BH")

esc_react_r$Holm =
  p.adjust(esc_react_r$pvalue,
           method = "holm")

esc_react_r$Hochberg =
  p.adjust(esc_react_r$pvalue,
           method = "hochberg")

esc_react_r$Hommel =
  p.adjust(esc_react_r$pvalue,
           method = "hommel")

esc_react_r$BY =
  p.adjust(esc_react_r$pvalue,
           method = "BY")

X = esc_react_r$pvalue
Y = cbind(esc_react_r$Bonferroni,
          esc_react_r$BH,
          esc_react_r$Holm,
          esc_react_r$Hochberg,
          esc_react_r$Hommel,
          esc_react_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-reaction-time-PS results only, ran")

# esc_neur_r #####

esc_neur_r <- esc_neur_r[ order( esc_neur_r$pvalue ), ]

esc_neur_r$Bonferroni =
  p.adjust(esc_neur_r$pvalue,
           method = "bonferroni")

esc_neur_r$BH =
  p.adjust(esc_neur_r$pvalue,
           method = "BH")

esc_neur_r$Holm =
  p.adjust(esc_neur_r$pvalue,
           method = "holm")

esc_neur_r$Hochberg =
  p.adjust(esc_neur_r$pvalue,
           method = "hochberg")

esc_neur_r$Hommel =
  p.adjust(esc_neur_r$pvalue,
           method = "hommel")

esc_neur_r$BY =
  p.adjust(esc_neur_r$pvalue,
           method = "BY")

X = esc_neur_r$pvalue
Y = cbind(esc_neur_r$Bonferroni,
          esc_neur_r$BH,
          esc_neur_r$Holm,
          esc_neur_r$Hochberg,
          esc_neur_r$Hommel,
          esc_neur_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-neuroticism-PS results only, ran")

# esc_schi_r #####

esc_schi_r <- esc_schi_r[ order( esc_schi_r$pvalue ), ]

esc_schi_r$Bonferroni =
  p.adjust(esc_schi_r$pvalue,
           method = "bonferroni")

esc_schi_r$BH =
  p.adjust(esc_schi_r$pvalue,
           method = "BH")

esc_schi_r$Holm =
  p.adjust(esc_schi_r$pvalue,
           method = "holm")

esc_schi_r$Hochberg =
  p.adjust(esc_schi_r$pvalue,
           method = "hochberg")

esc_schi_r$Hommel =
  p.adjust(esc_schi_r$pvalue,
           method = "hommel")

esc_schi_r$BY =
  p.adjust(esc_schi_r$pvalue,
           method = "BY")

X = esc_schi_r$pvalue
Y = cbind(esc_schi_r$Bonferroni,
          esc_schi_r$BH,
          esc_schi_r$Holm,
          esc_schi_r$Hochberg,
          esc_schi_r$Hommel,
          esc_schi_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-schizophrenia-PS results only, ran")

# esc_hei_r #####

esc_hei_r <- esc_hei_r[ order( esc_hei_r$pvalue ), ]

esc_hei_r$Bonferroni =
  p.adjust(esc_hei_r$pvalue,
           method = "bonferroni")

esc_hei_r$BH =
  p.adjust(esc_hei_r$pvalue,
           method = "BH")

esc_hei_r$Holm =
  p.adjust(esc_hei_r$pvalue,
           method = "holm")

esc_hei_r$Hochberg =
  p.adjust(esc_hei_r$pvalue,
           method = "hochberg")

esc_hei_r$Hommel =
  p.adjust(esc_hei_r$pvalue,
           method = "hommel")

esc_hei_r$BY =
  p.adjust(esc_hei_r$pvalue,
           method = "BY")

X = esc_hei_r$pvalue
Y = cbind(esc_hei_r$Bonferroni,
          esc_hei_r$BH,
          esc_hei_r$Holm,
          esc_hei_r$Hochberg,
          esc_hei_r$Hommel,
          esc_hei_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on escitalopram-height-PS results only, ran")

# nor_anor_r #####

nor_anor_r <- nor_anor_r[ order( nor_anor_r$pvalue ), ]

nor_anor_r$Bonferroni =
  p.adjust(nor_anor_r$pvalue,
           method = "bonferroni")

nor_anor_r$BH =
  p.adjust(nor_anor_r$pvalue,
           method = "BH")

nor_anor_r$Holm =
  p.adjust(nor_anor_r$pvalue,
           method = "holm")

nor_anor_r$Hochberg =
  p.adjust(nor_anor_r$pvalue,
           method = "hochberg")

nor_anor_r$Hommel =
  p.adjust(nor_anor_r$pvalue,
           method = "hommel")

nor_anor_r$BY =
  p.adjust(nor_anor_r$pvalue,
           method = "BY")

X = nor_anor_r$pvalue
Y = cbind(nor_anor_r$Bonferroni,
          nor_anor_r$BH,
          nor_anor_r$Holm,
          nor_anor_r$Hochberg,
          nor_anor_r$Hommel,
          nor_anor_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-anorexia-PS results only, ran")

# nor_react_r #####

nor_react_r <- nor_react_r[ order( nor_react_r$pvalue ), ]

nor_react_r$Bonferroni =
  p.adjust(nor_react_r$pvalue,
           method = "bonferroni")

nor_react_r$BH =
  p.adjust(nor_react_r$pvalue,
           method = "BH")

nor_react_r$Holm =
  p.adjust(nor_react_r$pvalue,
           method = "holm")

nor_react_r$Hochberg =
  p.adjust(nor_react_r$pvalue,
           method = "hochberg")

nor_react_r$Hommel =
  p.adjust(nor_react_r$pvalue,
           method = "hommel")

nor_react_r$BY =
  p.adjust(nor_react_r$pvalue,
           method = "BY")

X = nor_react_r$pvalue
Y = cbind(nor_react_r$Bonferroni,
          nor_react_r$BH,
          nor_react_r$Holm,
          nor_react_r$Hochberg,
          nor_react_r$Hommel,
          nor_react_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-reaction-time-PS results only, ran")

# nor_neur_r #####

nor_neur_r <- nor_neur_r[ order( nor_neur_r$pvalue ), ]

nor_neur_r$Bonferroni =
  p.adjust(nor_neur_r$pvalue,
           method = "bonferroni")

nor_neur_r$BH =
  p.adjust(nor_neur_r$pvalue,
           method = "BH")

nor_neur_r$Holm =
  p.adjust(nor_neur_r$pvalue,
           method = "holm")

nor_neur_r$Hochberg =
  p.adjust(nor_neur_r$pvalue,
           method = "hochberg")

nor_neur_r$Hommel =
  p.adjust(nor_neur_r$pvalue,
           method = "hommel")

nor_neur_r$BY =
  p.adjust(nor_neur_r$pvalue,
           method = "BY")

X = nor_neur_r$pvalue
Y = cbind(nor_neur_r$Bonferroni,
          nor_neur_r$BH,
          nor_neur_r$Holm,
          nor_neur_r$Hochberg,
          nor_neur_r$Hommel,
          nor_neur_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-neuroticism-PS results only, ran")

# nor_schi_r #####

nor_schi_r <- nor_schi_r[ order( nor_schi_r$pvalue ), ]

nor_schi_r$Bonferroni =
  p.adjust(nor_schi_r$pvalue,
           method = "bonferroni")

nor_schi_r$BH =
  p.adjust(nor_schi_r$pvalue,
           method = "BH")

nor_schi_r$Holm =
  p.adjust(nor_schi_r$pvalue,
           method = "holm")

nor_schi_r$Hochberg =
  p.adjust(nor_schi_r$pvalue,
           method = "hochberg")

nor_schi_r$Hommel =
  p.adjust(nor_schi_r$pvalue,
           method = "hommel")

nor_schi_r$BY =
  p.adjust(nor_schi_r$pvalue,
           method = "BY")

X = nor_schi_r$pvalue
Y = cbind(nor_schi_r$Bonferroni,
          nor_schi_r$BH,
          nor_schi_r$Holm,
          nor_schi_r$Hochberg,
          nor_schi_r$Hommel,
          nor_schi_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-schizophrenia-PS results only, ran")

# nor_hei_r #####

nor_hei_r <- nor_hei_r[ order( nor_hei_r$pvalue ), ]

nor_hei_r$Bonferroni =
  p.adjust(nor_hei_r$pvalue,
           method = "bonferroni")

nor_hei_r$BH =
  p.adjust(nor_hei_r$pvalue,
           method = "BH")

nor_hei_r$Holm =
  p.adjust(nor_hei_r$pvalue,
           method = "holm")

nor_hei_r$Hochberg =
  p.adjust(nor_hei_r$pvalue,
           method = "hochberg")

nor_hei_r$Hommel =
  p.adjust(nor_hei_r$pvalue,
           method = "hommel")

nor_hei_r$BY =
  p.adjust(nor_hei_r$pvalue,
           method = "BY")

X = nor_hei_r$pvalue
Y = cbind(nor_hei_r$Bonferroni,
          nor_hei_r$BH,
          nor_hei_r$Holm,
          nor_hei_r$Hochberg,
          nor_hei_r$Hommel,
          nor_hei_r$BY)

matplot(X, Y,
        xlab="Raw p-value",
        ylab="Adjusted p-value",
        type="l",
        asp=1,
        col=1:6,
        lty=1,
        lwd=2
        ,
        xlim=c(0,1),
        ylim=c(0,1)
)

legend('bottomright',
       legend = c("Bonferroni", "BH", "Holm", "Hochberg", "Hommel", "BY"),
       col = 1:6,
       cex = 1,   
       pch = 16)

abline(0, 1,
       col=1,
       lty=2,
       lwd=1)

title( main = "FDR conducted on nortriptyline-height-PS results only, ran")

dev.off()

##### Export a csv of all the results with their new, BH-adjusted p-values #####

# BH-adjusted p-values calculated all together #

sink( "pvalues_adjusted_simultaneously.csv")
cat( "," )
write.table( pvalues, col.names = TRUE, sep = ",")
sink()

# BH-adjusted p-values calculated separately - randomized only #
# N.B. no baseline analyses here #
# N.B. no "both" analyses here #
# N.B. no "randomized and non-randomized" analyses here #

pvalues_adjusted_for_ran_only_separate_drug_trait_pairs <-
  rbind(
    esc_anor_r,
    esc_react_r,
    esc_neur_r,
    esc_schi_r,
    esc_hei_r,
    nor_anor_r,
    nor_react_r,
    nor_neur_r,
    nor_schi_r,
    nor_hei_r
  )

sink( "pvalues_adjusted_separately.csv")
cat( "," )
write.table( pvalues_adjusted_for_ran_only_separate_drug_trait_pairs, col.names = TRUE, sep = ",")
sink()

# BH-adjusted p-values calculated separately - randomized and non-randomized #
# N.B. no baseline analyses here #
# N.B. no "both" analyses here #
# N.B. no "randomized only" analyses here #

pvalues_adjusted_for_ran_nonran_separate_drug_trait_pairs <-
  rbind(
    esc_anor_r_nr,
    esc_react_r_nr,
    esc_neur_r_nr,
    esc_schi_r_nr,
    esc_hei_r_nr,
    nor_anor_r_nr,
    nor_react_r_nr,
    nor_neur_r_nr,
    nor_schi_r_nr,
    nor_hei_r_nr
  )

sink( "pvalues_adjusted_separately_II.csv")
cat( "," )
write.table( pvalues_adjusted_for_ran_nonran_separate_drug_trait_pairs, col.names = TRUE, sep = ",")
sink()



##### Plot improvement (and appetite changes) in non-dropouts and dropouts, by anorexia PS #####

ggplot(nondrops_end,aes(y=zmadrs,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, non-dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(nondrops_end,aes(y=percchange,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, non-dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(subset(nondrops_end,percchange>-0.4),aes(y=percchange,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, non-dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(drop_esc_end,aes(y=zmadrs,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, escitalopram dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(drop_esc_end,aes(y=percchange,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, escitalopram dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(drop_nor_end,aes(y=zmadrs,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, nortriptyline dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)+
  scale_color_manual(values=c("green"),
                     labels=c("nortriptyline"))+
  theme_bw()

ggplot(drop_nor_end,aes(y=percchange,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, nortriptyline dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)+
  scale_color_manual(values=c("green"),
                     labels=c("nortriptyline"))+
  theme_bw()

range(nondrops_end$MADRSpercchangewk12)
hist(nondrops_end$MADRSpercchangewk12,breaks=30)
range(drop_esc_end$MADRSpercchangewk12)
hist(drop_esc_end$MADRSpercchangewk12,breaks=30)
range(drop_nor_end$MADRSpercchangewk12)
hist(drop_nor_end$MADRSpercchangewk12,breaks=30)

ggplot(nondrops,aes(y=asec8wk,x=AZ_0_1,col=drug))+
  geom_jitter()+
  geom_smooth(method="lm",size=2)+
  labs(title="Increased appetite severity, non-dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-2.5,2.5)+
  ylim(0,3)+
  scale_color_manual(values=c("magenta", "green"),
                              labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(nondrops,aes(y=asec9wk,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm",size=2)+
  labs(title="Decreased appetite severity, non-dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-3,3)+
  ylim(0,3)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(drop_esc,aes(y=asec8wk,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Increased appetite severity, escitalopram dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-2.5,2.5)+
  ylim(0,3)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(drop_esc,aes(y=asec9wk,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Decreased appetite severity, escitalopram dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-3,3)+
  ylim(0,3)+
  scale_color_manual(values=c("magenta", "green"),
                     labels=c("escitalopram", "nortriptyline"))+
  theme_bw()

ggplot(drop_nor,aes(y=asec8wk,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Increased appetite severity, nortriptyline dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-2.5,2.5)+
  ylim(0,3)+
  scale_color_manual(values=c("green"),
                     labels=c("nortriptyline"))+
  theme_bw()

ggplot(drop_nor,aes(y=asec9wk,x=AZ_0_1,col=drug))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="Decreased appetite severity, nortriptyline dropouts, by anorexia polygenic score (PT 0.1)")+
  xlim(-3,3)+
  ylim(0,3)+
  scale_color_manual(values=c("green"),
                     labels=c("nortriptyline"))+
  theme_bw()

##### Plot improvement in non-dropouts and dropouts, by reaction time PS #####

ggplot(nondrops_end,aes(y=zmadrs,x=CZ_0_05))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, non-dropouts, by reaction time polygenic score (PT 0.05)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(nondrops_end,aes(y=percchange,x=CZ_0_05))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, non-dropouts, by reaction time polygenic score (PT 0.05)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

ggplot(drop_esc_end,aes(y=zmadrs,x=CZ_0_05))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, escitalopram dropouts, by reaction time polygenic score (PT 0.05)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(drop_esc_end,aes(y=percchange,x=CZ_0_05))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, escitalopram dropouts, by reaction time polygenic score (PT 0.05)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

ggplot(drop_nor_end,aes(y=zmadrs,x=CZ_0_05))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, nortriptyline dropouts, by reaction time polygenic score (PT 0.05)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(drop_nor_end,aes(y=percchange,x=CZ_0_05))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, nortriptyline dropouts, by reaction time polygenic score (PT 0.05)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

##### Plot improvement in non-dropouts and dropouts, by schizophrenia PS #####

ggplot(nondrops_end,aes(y=zmadrs,x=SZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, non-dropouts, by schizophrenia polygenic score (PT 0.01)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(nondrops_end,aes(y=percchange,x=SZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, non-dropouts, by schizophrenia polygenic score (PT 0.01)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

ggplot(drop_esc_end,aes(y=zmadrs,x=SZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, escitalopram dropouts, by schizophrenia polygenic score (PT 0.01)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(drop_esc_end,aes(y=percchange,x=SZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, escitalopram dropouts, by schizophrenia polygenic score (PT 0.01)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

ggplot(drop_nor_end,aes(y=zmadrs,x=SZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, nortriptyline dropouts, by schizophrenia polygenic score (PT 0.01)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(drop_nor_end,aes(y=percchange,x=SZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, nortriptyline dropouts, by schizophrenia polygenic score (PT 0.01)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

##### Plot improvement in non-discontinuers and discontinuers, by height PS #####

ggplot(nondiscs_end,aes(y=zmadrs,x=HZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, non-discontinuers, by height polygenic score (PT 0.01)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(nondiscs_end,aes(y=percchange,x=HZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, non-discontinuers, by height polygenic score (PT 0.01)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

ggplot(disc_esc_end,aes(y=zmadrs,x=HZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, escitalopram discontinuers, by height polygenic score (PT 0.01)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(disc_esc_end,aes(y=percchange,x=HZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, escitalopram discontinuers, by height polygenic score (PT 0.01)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

ggplot(disc_nor_end,aes(y=zmadrs,x=HZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="End week z-MADRS, nortriptyline discontinuers, by height polygenic score (PT 0.01)")+
  xlim(-2.5,2.5)+
  ylim(-2.1,3.5)

ggplot(disc_nor_end,aes(y=percchange,x=HZ_0_01))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(title="MADRS per cent change, nortriptyline discontinuers, by height polygenic score (PT 0.01)")+
  xlim(-3,3)+
  ylim(-0.8,1.1)

##### Plot PS density plots for randomized, non-rand, dropouts, non-dropouts... #####

p1 <- ggplot()+
  geom_density(aes(GE_ran$AZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$AZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$AZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$AZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$AZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$AZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Anorexia nervosa polygenic score")+
  labs(title="Randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p2 <- ggplot()+
  geom_density(aes(GE_nran$AZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$AZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$AZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$AZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$AZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$AZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Anorexia nervosa polygenic score")+
  labs(title="Non-randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)


p3 <- ggplot()+
  geom_density(aes(GE_ran$CZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$CZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$CZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$CZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$CZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$CZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Reaction time polygenic score")+
  labs(title="Randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p4 <- ggplot()+
  geom_density(aes(GE_nran$CZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$CZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$CZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$CZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$CZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$CZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Reaction time polygenic score")+
  labs(title="Non-randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p5 <- ggplot()+
  geom_density(aes(GE_ran$NZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$NZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$NZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$NZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$NZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$NZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Neuroticism polygenic score")+
  labs(title="Randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p6 <- ggplot()+
  geom_density(aes(GE_nran$NZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$NZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$NZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$NZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$NZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$NZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Neuroticism polygenic score")+
  labs(title="Non-randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p7 <- ggplot()+
  geom_density(aes(GE_ran$SZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$SZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$SZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$SZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$SZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$SZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Schizophrenia polygenic score")+
  labs(title="Randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p8 <- ggplot()+
  geom_density(aes(GE_nran$SZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$SZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$SZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$SZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$SZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$SZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Schizophrenia polygenic score")+
  labs(title="Non-randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p9 <- ggplot()+
  geom_density(aes(GE_ran$HZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$HZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$HZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$HZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$HZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_ran$HZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Height polygenic score")+
  labs(title="Randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)

p10 <- ggplot()+
  geom_density(aes(GE_nran$HZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$HZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$HZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$HZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$HZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nran$HZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Height polygenic score")+
  labs(title="Non-randomized participants")+
  xlim(-5,5)+
  ylim(0,0.7)


p11 <- ggplot()+
  geom_density(aes(GE_drop$AZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$AZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$AZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$AZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$AZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$AZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Anorexia nervosa polygenic score")+
  labs(title="Dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p12 <- ggplot()+
  geom_density(aes(GE_nondrop$AZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$AZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$AZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$AZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$AZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$AZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Anorexia nervosa polygenic score")+
  labs(title="Non-dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)


p13 <- ggplot()+
  geom_density(aes(GE_drop$CZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$CZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$CZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$CZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$CZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$CZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Reaction time polygenic score")+
  labs(title="Dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p14 <- ggplot()+
  geom_density(aes(GE_nondrop$CZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$CZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$CZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$CZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$CZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$CZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Reaction time polygenic score")+
  labs(title="Non-dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p15 <- ggplot()+
  geom_density(aes(GE_drop$NZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$NZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$NZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$NZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$NZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$NZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Neuroticism polygenic score")+
  labs(title="Dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p16 <- ggplot()+
  geom_density(aes(GE_nondrop$NZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$NZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$NZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$NZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$NZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$NZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Neuroticism polygenic score")+
  labs(title="Non-dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p17 <- ggplot()+
  geom_density(aes(GE_drop$SZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$SZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$SZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$SZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$SZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$SZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Schizophrenia polygenic score")+
  labs(title="Dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p18 <- ggplot()+
  geom_density(aes(GE_nondrop$SZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$SZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$SZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$SZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$SZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$SZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Schizophrenia polygenic score")+
  labs(title="Non-dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p19<- ggplot()+
  geom_density(aes(GE_drop$HZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$HZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$HZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$HZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$HZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_drop$HZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Height polygenic score")+
  labs(title="Dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

p20 <- ggplot()+
  geom_density(aes(GE_nondrop$HZ_0_0001),col="red",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$HZ_0_01),col="orange",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$HZ_0_05),col="green",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$HZ_0_1),col="blue",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$HZ_0_5),col="purple",size=1,alpha=0.25)+
  geom_density(aes(GE_nondrop$HZ_1),col="magenta",size=1,alpha=0.25)+
  xlab("Height polygenic score")+
  labs(title="Non-dropouts")+
  xlim(-5,5)+
  ylim(0,0.7)

pdf( "PS_distributions.pdf")

grid.arrange( p1, p2, p11, p12)
grid.arrange( p3, p4, p13, p14 )
grid.arrange( p5, p6, p15, p16 )
grid.arrange( p7, p8, p17, p18 )
grid.arrange( p9, p10, p19, p20 )

dev.off()

