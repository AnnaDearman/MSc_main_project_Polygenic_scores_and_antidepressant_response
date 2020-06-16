# Adapted from JoniColeman's gwas repository

library( ggplot2 )

args <- commandArgs( TRUE )
root <- args[ 1 ]
PCAEVEC <- read.table( paste( root, ".1kg.LD_pop_strat.pca.evec_RENAMED", sep = "" ), head = TRUE )
colnames( PCAEVEC ) <- c( "ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Population" )

# Change "My_cohort" below to your project name

pops <- c( "ASW", "CEU", "CHB", "CHS", "CLM", "FIN", "GBR", "My_cohort", "IBS", "JPT", "LWK", "MXL", "PUR", "TSI", "YRI" )

# I wanted to customise this plot so changed the script to use
# the ggplot command

pdf( "My_cohort_vs_1KG_pop_strat_PCA.pdf" )

ggplot( PCAEVEC, aes( x = PC1, y = PC2, colour = Population ) )+
    geom_point( shape = 4 )+
    scale_colour_manual( name = "", labels = c( "African-American", "Utah NW European", "Han Chinese", "S Han Chinese", "Colombian", "Finnish", "British Eng Scot", "My_cohort", "Spanish Iberian", "Japanese", "Luhya Kenyan", "Mexican-American", "Puerto Rican", "Tuscan Italian", "Yoruba Nigerian" ), values = c( "ASW" = "goldenrod1", "CEU" = "greenyellow", "CHB" = "mediumorchid2", "CHS" = "mediumpurple4", "CLM" = "royalblue3", "FIN" = "darkcyan", "GBR" = "black", "Gendep" = "chartreuse4", "IBS" = "springgreen3", "JPT" = "deeppink", "LWK" = "salmon1", "MXL" = "blue1", "PUR" = "deepskyblue2", "TSI" = "darkgreen", "YRI" = "tomato") )+
    labs( x = "PC1", y = "PC2" )+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() )+
    theme( panel.background = element_rect( fill = "white" ) )

dev.off()