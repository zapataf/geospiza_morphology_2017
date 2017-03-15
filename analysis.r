library(tidyverse)
library(magrittr)
library(GGally)
library(clustvarsel)

# Load data -------------------------
raw_data <- read_csv( "./data/CAS_Swarth_Gesopiza.csv" )
raw_data
#View(original_finches_data) # View as table - interactive

# Wrangle data -------------------------

clean_data <- 
  raw_data %>%
    filter( Sex == "Male" ) %>% # Select only Male data
    select( -Gonys, -MiddleToe ) %>% # Remove traits not used in the analyses
    filter( !is.na( Wing ), 
            !is.na( Tail ), 
            !is.na( Blength ), 
            !is.na( Bdepth ), 
            !is.na( Bdepth ), 
            !is.na( Tarsus ) ) %>% # Remove records with NA from traits of interest
    distinct #%>% # Remove duplicated records
    #View

# Principal Component Analysis (PCA) on log (data) using covariance matrix -------------------------

pca_results <- 
  clean_data %>%
    select( Wing:Tarsus ) %>% # Select traits
    mutate_each( funs( log( . ) ) ) %>% # log-transform traits
    rename( LnWing = Wing, 
            LnTail = Tail, 
            LnBlength = Blength, 
            LnBdepth = Bdepth, 
            LnBwidth = Bwidth, 
            LnTarsus = Tarsus ) %>%
    prcomp( center = T, scale. = F )  # PCA using covariance matrix (ie., not scaled)

clean_data_pca <- bind_cols( clean_data, tbl_df( pca_results$x ) ) # Tibble with original data plus values for all PC Axes

# Variable Selection using all PC Axes -------------------------
# (Because "forward" search generated the same subset of variables (see Supplementary Material)
# here we carry out the analysis with backward search. For detail see ?clustvarsel)

clean_data_pca_varsel <-
  clean_data_pca %>%
    select( PC1:PC6 ) %>%
    clustvarsel( G = 1:30, search = c( "greedy" ), direction = c( "backward" ) )

#clean_data_pca_varsel$subset # Uncomment to see resulting subset (it should be: PC1, PC2, PC3, PC4)

# Gaussian Mixture Models on axes defined by PCA, and using variable selection -------------------------

mclust.options() # Check Mclust options
opt_mc <- mclust.options() # Save default
mclust.options( hcUse = "VARS" ) # Change default as needed

clean_data_pca_varsel_gmm <-
  clean_data_pca %>%
    select( PC1:PC4 ) %>%
    Mclust( G = 1:30 )

#summary( clean_data_pca_varsel_gmm ) # uncomment to see results from Mclust, number of groups, and model.
#help(mclustModelNames) # uncomment to get information about model names
#plot( clean_data_pca_varsel_gmm$BIC ) # uncomment to plot BIC support for all models

# Extract BIC values for the best model conditional on the number of groups
bic_best_model_per_group <- apply( clean_data_pca_varsel_gmm$BIC, 1, max, na.rm = T )

# Empirical support for the hypothesis of species limits by Lack (1947) -------------------------
h_lack <-
  clean_data_pca %>%
    select( PC1:PC4 ) %>%
    MclustDA( G = 1, class = clean_data_pca$Taxon, modelType = "MclustDA" )

#summary(h_lack) # Uncomment to see details
#h_lack$models # Uncomment to see details
#h_lack$bic # Uncomment to see BIC support

# Empirical support for the hypothesis of species limits based on current taxonomy -------------------------
h_current_taxonomy <-
  clean_data_pca %>%
    select( PC1:PC4 ) %>%
    MclustDA( G = 1, class = clean_data_pca$New_Taxonomy, modelType = "MclustDA" )

#summary(h_current_taxonomy) # Uncomment to see details
#h_current_taxonomy$models # Uncomment to see details
#h_current_taxonomy$bic # Uncomment to see BIC support

# Plot empirical support for different hypotheses of species limits -------------------------

as_tibble( bic_best_model_per_group ) %>%
              ggplot( aes( x = seq_along( value ), y = value ) ) +
                scale_x_continuous( breaks = 1:30 ) +
                xlab( "Number of morphological groups" ) +
                ylab( "Empirical support (BIC)" ) +              
                geom_point( color = "transparent" ) +
                #geom_point( subset(as_tibble( bic_best_model_per_group ), value == 3054.508), color="red" ) +
                geom_vline( xintercept = clean_data_pca_varsel_gmm$G, linetype = "dashed" ) +
                geom_vline( xintercept = length(unique(clean_data_pca$Taxon)), linetype = "dashed" ) +
                geom_vline( xintercept = length(unique(clean_data_pca$New_Taxonomy)), linetype = "dashed" ) +
                theme( axis.line = element_line( color = "black", size = 0.5), 
                       axis.title = element_text( size = 15 ), 
                       axis.text = element_text( size = 10 ),
                       panel.border = element_rect( color = "transparent", fill = NA ),
                       panel.grid.minor.y = element_line( size = 0.1, color = "grey" ),
                       panel.background = element_rect( fill = "transparent" ) )  

# Usueful summaries -------------------------

# Create tibble adding mclust classification result to full data

clean_data_pca_mclust <- 
  bind_cols( clean_data_pca, tibble( mcluster_classification = clean_data_pca_varsel_gmm$classification ) )

# How many islands per morphological group

clean_data_pca_mclust %>% 
  group_by( mcluster_classification ) %>% 
  summarise( n_distinct( Island ) ) #%>% # Uncomment here and below to see median
  #summarise(median(`n_distinct(Island)`))

# How many morphological groups per island

clean_data_pca_mclust %>% 
  group_by( Island ) %>% 
  summarise(n_distinct( mcluster_classification ) )

# Scatterplots  -------------------------

# Define colors
morphogroups_colors <- c( "#999999", 
                          "#E69F00", 
                          "#56B4E9", 
                          "#009E73", 
                          "#F0E442", 
                          "#0072B2", 
                          "#D55E00", 
                          "#CC79A7" ) 

# Plot
clean_data_pca_mclust %>% 
  ggplot( aes( x = PC1, y = PC2, color = factor( mcluster_classification ) ) ) +
    geom_point(size = 3, shape = 21, stroke = 1) +
    scale_color_manual( values = morphogroups_colors ) + 
    theme(axis.line = element_line( color = "black", size = 0.5),
          axis.title = element_text( size = 15 ), 
          axis.text = element_text( size = 14 ),
          panel.border = element_rect( color = "transparent", fill = NA ),
          panel.background = element_rect( fill = "transparent" ) )  
                                                                                                                                          
# TO DO: make this a loop to plot all by all PCs and then loadings.


# Compare assignment of specimens to groups between the best Mclust 
# model and slternative hypotheses -------------------------

clean_data_pca_mclust  %>% 
  ggplot(aes (mcluster_classification, fill = factor( mcluster_classification ) ) ) + 
    geom_histogram( binwidth = 1 ) + 
    scale_fill_manual( values = morphogroups_colors) +
    scale_x_continuous( breaks = 1:8 ) +
    facet_grid( ~ Taxon) + # Change to New_Taxonomy to generate histrogram for current taxonomy
    xlab( "Morphological groups" ) +
    ylab( "Specimens" ) +
    theme(axis.line = element_line( color = "black", size = 0.3),
          axis.text = element_text( size = 12 ),
          axis.title = element_text( size = 15 ),
          panel.border = element_rect( color = "transparent", fill = NA ),
          panel.grid.minor.y = element_line( size = 0.1, color = "grey" ),
          panel.background = element_rect( color = "grey" , fill = NA ), 
          legend.position = "none" )

# Exploratory Data Analysis (EDA) -------------------------

# Univariate plots
clean_data %>%
  ggplot(aes( Tail ) ) + # Change trait to be plotted
    geom_histogram( binwidth = 1 ) + # May need to adjust binwidth if traits changes
    facet_grid( ~ Taxon ) +
    labs( y = "Specimens" ) +
    theme(axis.line = element_line( color = "black", size = 0.3),
          axis.text = element_text( size = 14 ),
          axis.title = element_text( size = 15 ),
          panel.border = element_rect( color = "transparent", fill = NA ),
          panel.grid.minor.y = element_line( size = 0.1, color = "grey" ),
          panel.background = element_rect( color = "grey" , fill = NA)) 

# Bivariate plots
clean_data %>%
  ggplot(aes( x = Wing, y = Tail ) ) + # Change to traits of interest
    geom_point() + 
    facet_wrap( ~ Taxon ) +
    theme_bw() +
    theme( axis.text = element_text( size = 14 ),
           axis.title = element_text( size = 15 ) )

clean_data %>%
  ggplot( aes( x = Wing, y = Tail, color = Taxon ) ) + # Change to traits of interest
    geom_point() +
    theme_bw() +
    theme( axis.text = element_text( size = 14 ),
           axis.title = element_text( size = 15 ) )

# All pairs of traits
clean_data %>%  
  ggpairs(aes( color = Taxon, alpha = 0.7 ),
         columns = c( "Wing", "Tail", "Blength", "Bdepth", "Bwidth", "Tarsus" ),
         upper = list( continuous = "density" ),
         lower = list( continuous = "points" ),
         axisLabels = "show"
    ) +
    theme_bw()

# Variable Selection using all PC Axes and Forward algorithm -------------------------

clean_data_pca_varselfwd <-
  clean_data_pca %>%
  select( PC1:PC6 ) %>%
  clustvarsel( G = 1:30, search = c( "greedy" ), direction = c( "forward" ) )

#clean_data_pca_varselfwd$subset # Uncomment to see resulting subset 
