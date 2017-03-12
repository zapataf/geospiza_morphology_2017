library(tidyverse)
library(magrittr)
library(clustvarsel)
library(GGally)

# Load data ----------------------------
raw_data <- read_csv("./data/CAS_Swarth_Gesopiza.csv")
raw_data
#View(original_finches_data) # View as table - interactive

# Define plot colors ----------------------
morphogroups_colors <- c("#999999", 
                          "#E69F00", 
                          "#56B4E9", 
                          "#009E73", 
                          "#F0E442", 
                          "#0072B2", 
                          "#D55E00", 
                          "#CC79A7", 
                          "#000000")

# Wrangle data ------------------------

clean_data <- 
  raw_data %>%
    filter( Sex == "Male" ) %>% # Select only Male data
    #select( -Gonys, -MiddleToe ) %>% # Remove traits not used in the analyses
    filter( !is.na( Wing ), 
            !is.na( Tail ), 
            !is.na( Blength ), 
            !is.na( Bdepth ), 
            !is.na( Bdepth ), 
            !is.na( Tarsus ) ) %>% # Remove records with NA from traits of interest
    distinct #%>% # Remove duplicated records
    #View

# Exploratory Data Analysis (EDA) -------------------

# Univariate plots
ggplot(clean_data, aes(Tail) ) + # Change trait to be plotted
  geom_histogram( binwidth = 1 ) + # May need to adjust binwidth if traits changes
  facet_wrap(~ Taxon) +
  labs( y = "Specimens" ) +
  theme_bw()

# Bivariate plots
ggplot(clean_data) +
  geom_point(aes( x = Wing, y = Tail) ) + # Change traits to of interest
  facet_wrap(~ Taxon) +
  theme_bw()

ggplot(clean_data) +
  geom_point(aes( x = Wing, y = Tail, color = Taxon, alpha = 0.7) ) +
  theme_bw()

# All pairs of traits
ggpairs(clean_data, aes( color = Taxon, alpha = 0.7 ),
        columns = c("Wing", "Tail", "Blength", "Bdepth", "Bwidth", "Tarsus"),
        upper = list(continuous = "density"),
        lower = list(continuous = "points"),
        axisLabels ="show"
      ) +
  theme_bw()

# Principal Component Analysis (PCA) on log (data) using covariance matrix -------------------

pca_results <- 
  clean_data %>%
    #select( Wing:Tarsus ) %>% # Select traits
    select( Wing:Blength, Bdepth:Tarsus ) %>%
    mutate_each( funs( log( . ) ) ) %>% # log-transform traits
    rename( LnWing = Wing, 
            LnTail = Tail, 
            LnBlength = Blength, 
            LnBdepth = Bdepth, 
            LnBwidth = Bwidth, 
            LnTarsus = Tarsus ) %>%
    prcomp(center=T, scale.=F) # PCA using covariance matrix (ie., not scaled)

clean_data_pca <- bind_cols( clean_data, tbl_df( pca_results$x ) ) # Tibble with original data plus values for all PC Axes

# Variable Selection using all PC Axes  ------------------------------
# (see Supplmentary Material for "forward" search that generated the same subset of variables)

clean_data_pca_varsel <-
  clean_data_pca %>%
    select( PC1:PC6 ) %>%
    clustvarsel( G = 1:30, search = c( "greedy" ), direction = c( "backward" ) )

#clean_data_pca_varsel$subset # Uncomment to see resulting subset (it should be: PC1, PC2, PC3, PC4)

# Gaussian Mixture Models on axes defined by PCA, and using variable selection. -------------------------

mclust.options() # Check Mclust options
OptMc <- mclust.options() # Save default
mclust.options(hcUse="VARS") # Change default as needed

clean_data_pca_varsel_gmm <-
  clean_data_pca %>%
    select( PC1:PC4 ) %>%
    Mclust( G=1:30 )

#summary( clean_data_pca_varsel_gmm ) # uncomment to see results from Mclust, number of groups, and model.
#help(mclustModelNames) # uncomment to get information about model names

