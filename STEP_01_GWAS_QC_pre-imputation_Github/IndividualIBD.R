# Modified from JoniColeman's original script

library( ggplot2 )

args <- commandArgs( TRUE )
root <- args[ 1 ]
sigma <- as.numeric( args[ 2 ] )
U <- read.table( paste( root, ".IBD.genome", sep = "" ), head = TRUE ) # read in table
V_ONE <- with( U, data.frame( FID1, IID1, PI_HAT ) ) # get variables of interest - reference individual
V_TWO <- with( U, data.frame( FID2, IID2, PI_HAT) ) # get variables of interest - test individual
names( V_TWO ) <- c( "FID1", "IID1", "PI_HAT" )
V <- as.data.frame( rbind( V_ONE, V_TWO ) )
names( V ) <- c( "FID1", "IID1", "PI_HAT" )
W <- aggregate( V$PI_HAT, FUN = mean, by = list( V$FID1, V$IID1 ) ) # calculate average pi hat
names( W ) <- c( "FID", "IID", "MEAN_PI_HAT" ) # rename columns
X <- mean( W$MEAN_PI_HAT ) # calculate mean of average pi hats
Y <- sd( W$MEAN_PI_HAT ) # calculate standard deviation of average pi hats
Z <- X + ( sigma * Y ) # calculate threshold
sink( paste( root, ".IBD_INDIV.txt", sep = "" ) )
W # print average pi hats

# ORIGINAL SCRIPT EDITED TO REMOVE ROW NAMES

Wout <- subset( W, W$MEAN_PI_HAT >= Z )[ , 1:2 ] # subset for outliers
write.table( Wout, file = paste( root, ".IBD_INDIV_outliers.txt", sep = "" ), quote = FALSE, row.names = FALSE ) # write outliers to file

# ADDED PLOTS TO ORIGINAL SCRIPT

# Generate plots:

# Average pi-hats

pdf( "Histogram_of_mean_pi_hats.pdf" )

ggplot( data = W, aes( W$MEAN_PI_HAT )) +
    geom_histogram( fill = "deeppink3" )+
    theme( panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
    labs( title = "Mean pi-hat values" )+
    labs( x = "Mean pi-hat", y = "Count" )

# pi-hats for outliers

# Use a sigma of three

Z2 <- X + ( 3 * Y )
Wout_sig3 <- subset( W, W$MEAN_PI_HAT >= Z2 )
names( Wout_sig3 ) <- c( "FID1", "IID1", "MEAN_PI_HAT" )
Wout_sig3$IID1 <- factor( Wout_sig3$IID1 )
indivs <- levels( Wout_sig3$IID1 )

S <- subset( V, V$IID1 %in% indivs )
names( S ) <- c( "FID1", "IID1", "PI_HAT" )

# N.B. this colour vector was appropriate for the number of outliers in my data set
# Others may need to add or subtract colours

cols = c( "red", "orangered", "gold", "limegreen", "dodgerblue1", "midnightblue", "purple2" )

# If you have used a different SD cutoff, change the filename!

pdf( "Histograms_of_pi_hats_for_sigma-3_outliers.pdf" )

for ( i in 1:length( indivs ) ){
    title <- paste( "Pi-hat values for", indivs[ i ] )
    data <- S[ S$IID1==indivs[ i ], ]
    print( ggplot( data, aes( data$PI_HAT ) )+
        geom_histogram( fill = cols[ i ])+
        theme( panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
        labs( title = title )+
        labs( x = "Pi-hat", y = "Count" )
)
}

dev.off()
