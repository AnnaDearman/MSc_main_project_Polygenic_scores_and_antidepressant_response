library( ggplot2 )
library( avengeme )

# This script produces several power calculations using different estimated values

##### BINARY TRAINING TRAIT #####

vg1 <- 0.2 
cov <- 0.3 
alpha <- 0.01
pi0 <- 0.5
pupper=c( 0, 0.5 )
cases <- c( 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000,
            30000, 40000, 50000, 75000, 100000, 200000, 300000, 400000, 500000, 750000,
            1000000 )
controls <- c( 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000,
               30000, 40000, 50000, 75000, 100000, 200000, 300000, 400000, 500000, 750000,
               1000000 )
prevalence <- c( 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 )
prev <- c( "0.01", "0.02", "0.03", "0.04", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3" )

my_df <- data.frame( Individuals=numeric(),
                     h2SNP=numeric(),
                     Covariance=numeric(),
                     Alpha=numeric(),
                     Power=numeric() )

prev_col <- vector( mode="character" )

for ( i in 1:length( prevalence ) ){
  for ( j in 1:length( cases ) ){
    power <- round( polygenescore( nsnp=200000, n=c( cases[j]+controls[j], 669 ),
                            vg1=vg1, cov12=cov*vg1, pi0=pi0, pupper=pupper,
                            nested=TRUE, weighted=TRUE, binary=c( TRUE, FALSE ),
                            prevalence=prevalence[i],
                            sampling=c( cases/( cases+controls ), prevalence ),
                            lambdaS=NA, shrinkage=FALSE, logrisk=FALSE,
                            alpha=0.01, r2gx=0, corgx=0, r2xy=0,
                            adjustedEffects=FALSE, riskthresh=0.1 )$power, 2 )
    my_df[nrow(my_df)+1,] <- c( cases[j]+controls[j], vg1, cov, alpha, power)
    prev_col <- append( prev_col, prev[i] )
  }
  }

my_df <- cbind( my_df, prev_col )
colnames( my_df ) <- c( "Individuals", "h2SNP", "Covariance",
                        "Alpha", "Power", "Prevalence" )

pdf( "Power_calculation_binary.pdf" )

ggplot( my_df, aes( x = Individuals, y = Power, col = Prevalence ) )+
  geom_point()+
  geom_line()+
  labs( title = "Predicted power as a function of n individuals in training sample\nand prevalence of training trait (binary)" )+
  labs( subtitle = paste( "h2SNP: ", vg1, "\nCovariance: ", cov, ", i.e. ", vg1, " x ", cov,
                         "\nAlpha: ", alpha, "\nProportion of uninformative SNPs: ", pi0,
                         "\nUpper bound of SNP p-value inclusion threshold: ", pupper[2], sep = "" ) )+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

##### CONTINUOUS TRAINING TRAIT #####

vg1 <- 0.2
cov <- 0.3
alpha <- 0.01
pi0 <- 0.5
pupper=c( 0, 0.5 )
cases <- c( 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 
            30000, 40000, 50000, 75000, 100000, 200000, 300000, 400000, 500000, 750000,
            1000000 )
controls <- c( 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 
               30000, 40000, 50000, 75000, 100000, 200000, 300000, 400000, 500000, 750000,
               1000000 )

my_df <- data.frame( Individuals=numeric(),
                     h2SNP=numeric(),
                     Covariance=numeric(),
                     Alpha=numeric(),
                     Power=numeric() )

prev_col <- vector( mode="character" )

for ( j in 1:length( cases ) ){
  power <- round( polygenescore( nsnp=200000, n=c( cases[j]+controls[j], 669 ), 
                                 vg1=vg1, cov12=cov*vg1, pi0=pi0, pupper=pupper,
                                 nested=TRUE, weighted=TRUE, binary=c( FALSE, FALSE ),
                                 lambdaS=NA, shrinkage=FALSE, logrisk=FALSE,
                                 alpha=0.01, r2gx=0, corgx=0, r2xy=0,
                                 adjustedEffects=FALSE, riskthresh=0.1 )$power, 2 )
  my_df[nrow(my_df)+1,] <- c( cases[j]+controls[j], vg1, cov, alpha, power)
}

pdf( "Power_calculation_continuous.pdf" )

ggplot( my_df, aes( x = Individuals, y = Power ) )+
  geom_point()+
  geom_line()+
  labs( title = "Predicted power as a function of n individuals in training sample\n(continuous training trait)" )+
  labs( subtitle = paste( "h2SNP: ", vg1, "\nCovariance: ", cov, ", i.e. ", vg1, " x ", cov,
                          "\nAlpha: ", alpha, "\nProportion of uninformative SNPs: ", pi0,
                          "\nUpper bound of SNP p-value inclusion threshold: ", pupper[2], sep = "" ) )+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

