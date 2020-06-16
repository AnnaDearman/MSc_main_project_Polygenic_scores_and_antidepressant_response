library( ggplot2 )
library( avengeme )

##### PARAMETERS TO KEEP CONSTANT #####

test <- 700 # Number of individuals in test data 
weighted <- TRUE # TRUE if estimated effect sizes are used as weights
shrinkage	<- FALSE # TRUE if effect sizes are to be shrunk to BLUPs
r2gx <- 0 # Proportion of variance in ERS explained by genetic effects in training sample
corgx	<- 0 # Genetic correlation between ERS and training trait
r2xy <- 0 # Proportion of variance in training trait explained by ERS
adjustedEffects	<- FALSE # TRUE if PRS and ERS are combined as a weighted sum.
riskthresh <- 0.1 # Absolute risk threshold for calculating net reclassification
# index. Default 0.1.
lambdaS	<- NA # Sibling relative recurrence risk in training sample.
pupper <- c ( 5*10^-8, 1*10^-6, 1*10^-4, 1*10^-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0 )
# Vector of p-value thresholds for selecting markers from training sample.
nested <- TRUE # TRUE if they have the same lower bound, which is the first element of pupper.

##### TRAINING-DATA-SPECIFIC PARAMETERS #####

study <- "Major_depression"
nsnp <- 200000 # Number of SNPs in summary statistics file
cases <- 45396 # Number of cases used to generate summ stats (check ancestry in paper)
ctrls <- 97250 # Number of controls used to generate summ stats (check ancestry in paper)
binary <- c( TRUE, FALSE ) # Training trait, then target trait
prevalence <- 0.15 
logrisk	<- FALSE # TRUE if binary trait arises from log-risk model rather than 
# liability threshold.

##### TRY DIFFERENT VALUES #####

vg1	<- c(0.1, 0.2, 0.5)
cov12	<- c(0.1, vg1[1]*0.1, 0.5) # Covariance between genetic effect sizes in the two samples. If
# the effects are fully correlated then cov12<=sqrt(vg1). If the effects 
# are identical then cov12=vg1 (default). In worked example, they use vg1*vg1,
# not sqrt(vg1)   
pi0	<- c(0.01, 0.5, 0.99) # Proportion of markers with no effect on the training trait. 
# Default 0, which means all SNPs were relevant.
alpha	<- c( 0.050, 0.010, 0.001 )# Significance level for testing association of the polygenic score 
# in the target sample. Default 0.05.

##### AUTOMATICALLY CALCULATED #####

n	<- c( cases+ctrls, test ) # Number of individuals in training, test data
sampling <- c( cases/(cases+ctrls), prevalence ) 
# For a binary trait, case/control sampling fraction in the 
# training sample. By default, sampling equals the prevalence, as in a cohort 
# study. If the sampling fraction is different in the target sample, sampling 
# should be a vector with two elements for the training and target samples 
# respectively. Default = "prevalence".

# Simply produce a power calculation (edit values below):

polygenescore(75000, 2000, vg1=0.2, cov12=0.3*0.2, pi0=0.95, pupper=c(0,1),
              nested=TRUE, weighted=TRUE, binary=c(TRUE,FALSE),
              prevalence=0.01, sampling=c(1000/2000,prevalence), lambdaS=NA,
              shrinkage=FALSE, logrisk=FALSE, alpha=0.01, r2gx=0,
              corgx=0, r2xy=0, adjustedEffects=FALSE, riskthresh=0.1)$power

# Produce a PDF showing different values for different options

pdf( paste( study, "_PolyGeneScore-Plots_prev-",
            prevalence, "_nSNP-", nsnp, "_cases-",
            cases, "_ctrls-", ctrls, "_.pdf",
            sep = "" ) )

pi0_list <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[1], pupper,
                          nested, weighted, binary,
                          prevalence, sampling, lambdaS,
                          shrinkage, logrisk, alpha[1], r2gx,
                          corgx, r2xy, adjustedEffects, riskthresh)
pi0_list2 <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[2], pupper,
                           nested, weighted, binary,
                           prevalence, sampling, lambdaS,
                           shrinkage, logrisk, alpha[1], r2gx,
                           corgx, r2xy, adjustedEffects, riskthresh)
pi0_list3 <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[3], pupper,
                           nested, weighted, binary,
                           prevalence, sampling, lambdaS,
                           shrinkage, logrisk, alpha[1], r2gx,
                           corgx, r2xy, adjustedEffects, riskthresh)
pi0_df_1 <- cbind( pi0_list$R2, pi0_list$NCP, pi0_list$p, pi0_list$power,
                   pi0_list$FDR,pi0_list$MSE, as.numeric(pupper[ 2:length(pupper)]))
pi0_df_2 <- cbind( pi0_list2$R2, pi0_list2$NCP, pi0_list2$p, pi0_list2$power,
                   pi0_list2$FDR, pi0_list2$MSE, as.numeric(pupper[ 2:length(pupper)]))
pi0_df_3 <- cbind( pi0_list3$R2, pi0_list3$NCP, pi0_list3$p, pi0_list3$power,
                   pi0_list3$FDR, pi0_list3$MSE, as.numeric(pupper[ 2:length(pupper)]))
pi0_df <- data.frame( rbind( pi0_df_1, pi0_df_2, pi0_df_3 ) )
pi0_df <- cbind( pi0_df, 
                 c( rep(as.character(pi0[1]), 9), rep(as.character(pi0[2]), 9), rep(as.character(pi0[3]), 9)))
colnames( pi0_df ) <- c( "R2", "NCP", "p", "power", "FDR", "MSE", "Thresholds", "Value")

ggplot( data=pi0_df, aes( x = Thresholds, y = FDR, col = Value))+
  geom_point()+
  geom_line()+
  labs(title = paste( study, ": Expected proportion of false positives
at different estimates of proportion of informative SNPs",
                      sep = "" ) )+
  labs(subtitle = paste( "Variance explained by genetic effects in training sample: ",
                         vg1[1], "\nCovariance between genetic effect sizes in two samples: ",
                         cov12[1], "\nAlpha: ", alpha[1],
                         sep=""))

alpha_list <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[1], pupper,
                            nested, weighted, binary,
                            prevalence, sampling, lambdaS,
                            shrinkage, logrisk, alpha[1], r2gx,
                            corgx, r2xy, adjustedEffects, riskthresh)
alpha_list2 <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[1], pupper,
                             nested, weighted, binary,
                             prevalence, sampling, lambdaS,
                             shrinkage, logrisk, alpha[2], r2gx,
                             corgx, r2xy, adjustedEffects, riskthresh)
alpha_list3 <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[1], pupper,
                             nested, weighted, binary,
                             prevalence, sampling, lambdaS,
                             shrinkage, logrisk, alpha[3], r2gx,
                             corgx, r2xy, adjustedEffects, riskthresh)
alpha_df_1 <- cbind( alpha_list$R2, alpha_list$NCP, alpha_list$p, alpha_list$power,
                     alpha_list$FDR, alpha_list$MSE, as.numeric(pupper[ 2:length(pupper)]))
alpha_df_2 <- cbind( alpha_list2$R2, alpha_list2$NCP, alpha_list2$p, alpha_list2$power,
                     alpha_list2$FDR, alpha_list2$MSE, as.numeric(pupper[ 2:length(pupper)]))
alpha_df_3 <- cbind( alpha_list3$R2, alpha_list3$NCP, alpha_list3$p, alpha_list3$power,
                     alpha_list3$FDR, alpha_list3$MSE, as.numeric(pupper[ 2:length(pupper)]))
alpha_df <- data.frame( rbind( alpha_df_1, alpha_df_2, alpha_df_3 ) )
alpha_df <- cbind( alpha_df, c( rep(as.character(alpha[1]), 9), rep(as.character(alpha[2]), 9), rep(as.character(alpha[3]), 9)))
colnames( alpha_df ) <- c( "R2", "NCP", "p", "power", "FDR", "MSE", "Thresholds", "Alpha_value")


ggplot( alpha_df, aes( x = Thresholds, y = power, col = Alpha_value))+
  geom_point()+
  geom_line()+
  labs(title = paste( study, ": Expected power of chisq test of association between PRS
and target trait, at different significance values",
                      sep = "" ))+
  labs(subtitle = paste( "Variance explained by genetic effects in training sample: ",
                         vg1[1], "\nCovariance between genetic effect sizes in two samples: ",
                         cov12[1], "\nProportion of informative SNPs: ", pi0[1],
                         sep=""))

vg1_list <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[1], pupper,
                          nested, weighted, binary,
                          prevalence, sampling, lambdaS,
                          shrinkage, logrisk, alpha[1], r2gx,
                          corgx, r2xy, adjustedEffects, riskthresh)
vg1_list2 <- polygenescore(nsnp, n, vg1[2], cov12[1], pi0[1], pupper,
                           nested, weighted, binary,
                           prevalence, sampling, lambdaS,
                           shrinkage, logrisk, alpha[1], r2gx,
                           corgx, r2xy, adjustedEffects, riskthresh)
vg1_list3 <- polygenescore(nsnp, n, vg1[3], cov12[1], pi0[1], pupper,
                           nested, weighted, binary,
                           prevalence, sampling, lambdaS,
                           shrinkage, logrisk, alpha[1], r2gx,
                           corgx, r2xy, adjustedEffects, riskthresh)
vg1_df_1 <- cbind( vg1_list$R2, vg1_list$NCP, vg1_list$p, vg1_list$power,
                   vg1_list$FDR, vg1_list$MSE, as.numeric(pupper[ 2:length(pupper)]))
vg1_df_2 <- cbind( vg1_list2$R2, vg1_list2$NCP, vg1_list2$p, vg1_list2$power,
                   vg1_list2$FDR, vg1_list2$MSE, as.numeric(pupper[ 2:length(pupper)]))
vg1_df_3 <- cbind( vg1_list3$R2, vg1_list3$NCP, vg1_list3$p, vg1_list3$power,
                   vg1_list3$FDR, vg1_list3$MSE, as.numeric(pupper[ 2:length(pupper)]))
vg1_df <- data.frame( rbind( vg1_df_1, vg1_df_2, vg1_df_3 ) )
vg1_df <- cbind( vg1_df, c( rep(as.character(vg1[1]), 9), rep(as.character(vg1[2]), 9),
                            rep(as.character(vg1[3]), 9)))
colnames( vg1_df ) <- c( "R2", "NCP", "p", "power", "FDR", "MSE", "Thresholds", "Var_explained_value")


ggplot( vg1_df, aes( x = Thresholds, y = FDR, col = Var_explained_value))+
  geom_point()+
  geom_line()+
  labs(title = paste( study, ": Expected proportion of false positives
at different values of genetic heritability",
                      sep = "" ))+
  labs(subtitle = paste( "Alpha: ", alpha[1], "\nCovariance between genetic effect sizes in two samples: ",
                         cov12[1], "\nProportion of informative SNPs: ", pi0[1],
                         sep=""))

cov12_list <- polygenescore(nsnp, n, vg1[1], cov12[1], pi0[1], pupper,
                            nested, weighted, binary,
                            prevalence, sampling, lambdaS,
                            shrinkage, logrisk, alpha[1], r2gx,
                            corgx, r2xy, adjustedEffects, riskthresh)
cov12_list2 <- polygenescore(nsnp, n, vg1[1], cov12[2], pi0[1], pupper,
                             nested, weighted, binary,
                             prevalence, sampling, lambdaS,
                             shrinkage, logrisk, alpha[1], r2gx,
                             corgx, r2xy, adjustedEffects, riskthresh)
cov12_list3 <- polygenescore(nsnp, n, vg1[1], cov12[3], pi0[1], pupper,
                             nested, weighted, binary,
                             prevalence, sampling, lambdaS,
                             shrinkage, logrisk, alpha[1], r2gx,
                             corgx, r2xy, adjustedEffects, riskthresh)
cov12_df_1 <- cbind( cov12_list$R2, cov12_list$NCP, cov12_list$p, cov12_list$power,
                     cov12_list$FDR, cov12_list$MSE, as.numeric(pupper[ 2:length(pupper)]))
cov12_df_2 <- cbind( cov12_list2$R2, cov12_list2$NCP, cov12_list2$p, cov12_list2$power,
                     cov12_list2$FDR, cov12_list2$MSE, as.numeric(pupper[ 2:length(pupper)]))
cov12_df_3 <- cbind( cov12_list3$R2, cov12_list3$NCP, cov12_list3$p, cov12_list3$power,
                     cov12_list3$FDR, cov12_list3$MSE, as.numeric(pupper[ 2:length(pupper)]))
cov12_df <- data.frame( rbind( cov12_df_1, cov12_df_2, cov12_df_3 ) )
cov12_df <- cbind( cov12_df, c( rep(as.character(cov12[1]), 9), rep(as.character(cov12[2]), 9),
                                rep(as.character(cov12[3]), 9)))
colnames( cov12_df ) <- c( "R2", "NCP", "p", "power", "FDR", "MSE", "Thresholds", "Covariance_value")

ggplot( cov12_df, aes( x = Thresholds, y = R2, col = Covariance_value))+
  geom_point()+
  geom_line()+
  labs(title = paste( study, ": R2 (squared correlation between PRS and target trait)
at different covariance values for genetic effect sizes in the two samples",
                      sep = "" ))+
  labs(subtitle = paste( "Alpha: ", alpha[1], "Variance explained by genetic effects in training sample: ",
                         vg1[1], "\nProportion of informative SNPs: ", pi0[1],
                         sep=""))
ggplot( cov12_df, aes( x = Thresholds, y = NCP, col = Covariance_value))+
  geom_point()+
  geom_line()+
  labs(title = paste( study, ": NCP (Non-centrality parameter of chisq test of association
between PRS and target trait) at different covariance values for genetic effect sizes in the two samples",
                      sep = "" ))+
  labs(subtitle = paste( "Alpha: ", alpha[1], "\nVariance explained by genetic effects in training sample: ",
                         vg1[1], "\nProportion of informative SNPs: ", pi0[1],
                         sep=""))
ggplot( cov12_df, aes( x = Thresholds, y = p, col = Covariance_value))+
  geom_point()+
  geom_line()+
  labs(title = paste( study, ": Expected P-value of chisq test of association between PRS
and target trait, at different covariance values for genetic effect sizes in the two samples",
                      sep = "" ))+
  labs(subtitle = paste( "Alpha: ", alpha[1], "\nVariance explained by genetic effects in training sample: ",
                         vg1[1], "\nProportion of informative SNPs: ", pi0[1],
                         sep=""))
ggplot( cov12_df, aes( x = Thresholds, y = power, col = Covariance_value))+
  geom_point()+
  geom_line()+
  labs(title = paste( study, ": Expected power of chisq test of association between PRS
and target trait, at different covariance values for genetic effect sizes in the two samples",
                      sep = "" ))+
  labs(subtitle = paste( "Alpha: ", alpha[1], "\nVariance explained by genetic effects in training sample: ",
                         vg1[1], "\nProportion of informative SNPs: ", pi0[1],
                         sep=""))

dev.off()

