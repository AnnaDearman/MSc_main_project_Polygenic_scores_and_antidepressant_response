# Adapted from original script in JoniColeman's gwas repository

args <- commandArgs( TRUE )
root <- args[ 1 ]

ibc_file <- paste( root, ".ibc", sep = "" )

# MODIFIED FROM ORIGINAL:
# Omitted step where sdcut is assigned as 3, so user can choose
# a cutoff in the command line

# get arguments
t = commandArgs()

# If flag "--args" was used in the command line, take the arguments passed after it and save as "args"
if ( charmatch( "--args", t, nomatch = -1 ) >= 0 ) args <- t[ ( ( 1:length( t ) )[ t == "--args" ] + 1 ):length( t ) ] else args <- ""

# If "ibc_file=" was passed as an argument, overwrite ibc_file with the name passed in command line
if ( charmatch( "ibc_file=", args, nomatch = -1 ) >= 0 ) ibc_file <- strsplit( args[ charmatch( "ibc_file=", args ) ], split = "=" )[[ 1 ]][ 2 ]

# MODIFIED FROM ORIGINAL:
# If "sdcut=" was passed as an argument, use this number as the SD cutoff,
# otherwise use 3
if ( charmatch( "sdcut=", args, nomatch = -1 ) >= 0 ) sdcut <- strsplit( args[ charmatch( "sdcut=", args ) ], split = "=" )[[ 1 ]][ 2 ] else sdcut <- 3

d <- read.table( ibc_file, head = TRUE );
het_outliers <- abs( scale( d$Fhat2 ) ) > sdcut
write.table( d[ het_outliers, ], file = paste( root, ".LD_het_outliers.txt", sep = "" ), sep = "\t", quote = FALSE, row.names = FALSE );
write.table( d[ het_outliers, c( 1, 2 ) ], file = paste( root, ".LD_het_outliers_sample_exclude", sep = "" ), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE );

# EXTRA PLOTTING FUNCTION ADDED:

meen <- mean( d$Fhat2 )
sdee <- as.numeric( sd( d$Fhat2 ) )
cutoff <- as.numeric( sdcut )*as.numeric( sdee )
cutoff2 <- 0-as.numeric( cutoff )

plot <- ggplot( data = d, aes( x = d$Fhat2 ))+
  geom_histogram()+
  labs( x = "Fhat2",
        y = "Count")+
  geom_vline( xintercept = meen, colour = "blue" )+
  geom_vline( xintercept = cutoff, colour = "red" )+
  geom_vline( xintercept = cutoff2, colour = "red" )

pdf( "Heterozygosity.pdf" )
print( plot )
dev.off()
