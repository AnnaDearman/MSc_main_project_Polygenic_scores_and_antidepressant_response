# This script was written largely for fun / practise. 
# The plot is somewhat useful for seeing the distribution
# of your allele frequency differences,
# e.g. if you use the wrong reference population, your distribution
# will probably be a lot broader than if you use the correct one.

library( ggplot2 )
library( dplyr )

freqs <- read.table( "[filepath and filename of FreqPlot .txt file]" )

colnames( freqs ) <- c( "SNP", "1KG MAF", "My cohort MAF", "Diff", "col5" )

# Plotting histogram of the differences between cohorts MAFs

freq_diff <- ggplot( freqs, aes( x = freqs[ , 4 ] ) )+
  geom_histogram( fill = "#FF33FF" )+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank() )+
  labs( x = "Difference",
        y = "Count",
        title = "Differences in minor allele frequencies between cohorts" )+
  geom_vline( xintercept = 0, colour = "blue", size = 1.2 )+
  geom_vline( xintercept = -0.2, colour = "red", size = 1.2 )+
  geom_vline( xintercept = 0.2, colour = "red", size = 1.2 )

pdf( "Diff_MAF_hist.pdf" )
print( freq_diff )
dev.off()

# The frequency difference is positive if the MAF is higher in the reference cohort
# Thus we have more SNPs with smaller MAFs in our cohort than we do
# SNPs with larger MAFs

outliers <- filter( freqs, Diff > 0.2 | Diff < -0.2 )
colnames( outliers ) <- c( "SNP", "1KG MAF", "My cohort MAF", "Diff", "col5" )

outliers <- outliers[ order( outliers$Diff ),]

write.table( outliers, file = "MAF_outliers.txt", quote = FALSE, row.names = FALSE )
