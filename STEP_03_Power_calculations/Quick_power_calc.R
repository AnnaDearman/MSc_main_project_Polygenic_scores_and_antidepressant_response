library( avengeme )

# To quickly calculate the power, change some values below:

# Binary

cases <- 329821
controls <- 0
vg1 <- 0.15
prevalence <- 0.01

polygenescore( nsnp = 200000, n = c( cases+controls, 4000 ), vg1 = vg1, 
               cov12 = 0.5*vg1, pi0 = 0.01, 
               binary = c( FALSE, FALSE ), prevalence = prevalence,
               sampling = c( cases/( cases+controls ), prevalence ), 
               alpha = 0.01 )$power

# Continuous

n <- 10000
vg1 <- 0.10
nsnp <- 200000

polygenescore( nsnp = nsnp, n = c( n, 669 ), vg1 = vg1, 
               cov12 = 0.5*vg1, pi0 = 0.01, 
               binary = c( FALSE, FALSE ), 
               alpha = 0.01 )$power
