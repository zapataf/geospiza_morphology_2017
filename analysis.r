# This script includes the R code used in the analyses presented in:
#
# 	Cadena DC, Jimenez I, Zapata F. Atalantean evolution in Darwin's Finches - 
#	Issues and Perspectives in Species Delimitation using Phenotypic Data
#	http://biorxiv.org/content/early/2017/04/06/124610
#
# This code was written collaboratively by C.D. Cadena, I. Jimenez and F. Zapata
# 
#
# Code sections: preliminaries,  data organization, analyses, plots, data summaries,
# and exploratory data analyses. Each subsection is separated by -------------------------
#


###################
## PRELIMINARIES ##
###################

# libraries used in analyses and plots
library(tidyverse)
library(magrittr)
library(GGally)
library(clustvarsel)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(legendMap)

#######################
## DATA ORGANIZATION ##
#######################

# Load morphological data -------------------------

data = read_csv( "./CAS_Swarth_Geospiza.csv" )
data
#View(data) # View as table - interactive

# Load raster images for map  -------------------------

# Get the following raster files (*.tif) from the ASTER elevation model
# available at: https://reverb.echo.nasa.gov/reverb/

#N00W090 <- raster("ASTGTM2_N00W090_dem.tif")
#N00W091 <- raster("ASTGTM2_N00W091_dem.tif")
#N00W092 <- raster("ASTGTM2_N00W092_dem.tif")
#N01W092 <- raster("ASTGTM2_N01W092_dem.tif")
#N01W093 <- raster("ASTGTM2_N01W093_dem.tif")
#S01W091 <- raster("ASTGTM2_S01W091_dem.tif")
#S01W090 <- raster("ASTGTM2_S01W090_dem.tif")
#S01W092 <- raster("ASTGTM2_S01W092_dem.tif")
#S02W090 <- raster("ASTGTM2_S02W090_dem.tif")
#S02W091 <- raster("ASTGTM2_S02W091_dem.tif")
#S02W092 <- raster("ASTGTM2_S02W092_dem.tif")

# Create a single raster of islands except islands in the NW; 
# there is not enough resolution in the raster file for these islands -------------------------

allgalapagos_NO_NW = merge( N00W090,
			N00W091, 
			N00W092, 
			S01W091, 
			S01W090, 
			S01W092, 
			S02W090, 
			S02W091, 
			S02W092 )

###################
## DATA ANALYSES ##
###################

# Principal Component Analysis (PCA) on log (data) using covariance matrix -------------------------

pca_results = 
  data %>%
    dplyr::select( Wing:Tarsus ) %>% # Select traits
    mutate_each( funs( log( . ) ) ) %>% # log-transform traits
    rename( LnWing = Wing, 
            LnTail = Tail, 
            LnBlength = Blength, 
            LnBdepth = Bdepth, 
            LnBwidth = Bwidth, 
            LnTarsus = Tarsus ) %>%
    prcomp( center = T, scale. = F )  # PCA using covariance matrix (ie., not scaled)

data_pca = bind_cols( data, tbl_df( pca_results$x ) ) # Tibble with original data plus values for all PC Axes

# Variable Selection using all PC Axes -------------------------
# (Because "forward" search generated the same subset of variables (see Exploratory Analyses below),
# we carried out the analyses with backward search. For details see ?clustvarsel)

data_pca_varsel =
  data_pca %>%
    dplyr::select( PC1:PC6 ) %>%
    clustvarsel( G = 1:30, search = c( "greedy" ), direction = c( "backward" ) )

#data_pca_varsel$subset # Uncomment to see resulting subset (it should be: PC1, PC2, PC3, PC4)

# Gaussian Mixture Models on axes defined by PCA, and using variable selection -------------------------

mclust.options() # Check Mclust options. For details, see ?mclust.options 
opt_mc = mclust.options() # Save default
mclust.options( hcUse = "VARS" ) # Change default to VARS (ie, original variables), or as needed.

data_pca_varsel_gmm =
  data_pca %>%
    dplyr::select( PC1:PC4 ) %>%
    Mclust( G = 1:30 )

#summary( data_pca_varsel_gmm ) # uncomment to see results from Mclust, number of groups, and model.
#help(mclustModelNames) # uncomment to get information about model names
#plot( data_pca_varsel_gmm$BIC ) # uncomment to plot BIC support for all models

# Extract BIC values for the best model conditional on the number of groups
bic_best_model_per_group = apply( data_pca_varsel_gmm$BIC, 1, max, na.rm = T )

# Because model with 7 groups is equally supported, 
# this is the classification of specimens to 7 groups --------------------

data_pca_varsel_gmm_7groups =
  data_pca %>%
  dplyr::select( PC1:PC4 ) %>%
  Mclust( G = 7 )

#summary(data_pca_varsel_gmm_7groups) # Uncomment to see details
#data_pca_varsel_gmm_7groups$bic # 

# Empirical support for the hypothesis of species limits by Lack (1947) -------------------------
h_lack =
  data_pca %>%
    dplyr::select( PC1:PC4 ) %>%
    MclustDA( G = 1, class = data_pca$Taxon, modelType = "MclustDA" )

#summary(h_lack) # Uncomment to see details
#h_lack$models # Uncomment to see details
#h_lack$bic # Uncomment to see BIC support

# Empirical support for the hypothesis of species limits based on current taxonomy -------------------------
h_current_taxonomy =
  data_pca %>%
    dplyr::select( PC1:PC4 ) %>%
    MclustDA( G = 1, class = data_pca$New_Taxonomy, modelType = "MclustDA" )

#summary(h_current_taxonomy) # Uncomment to see details
#h_current_taxonomy$models # Uncomment to see details
#h_current_taxonomy$bic # Uncomment to see BIC support

# Consolidate results in a sigle object -------------------------

# Add best mclust classification result (with 8 and 7 groups) to data tibble 

data_pca_mclust = 
  bind_cols( data_pca, 
             tibble( mcluster_classification = data_pca_varsel_gmm$classification ), 
             tibble( mcluster_classification7 = data_pca_varsel_gmm_7groups$classification ) )


###########
## PLOTS ##
###########

# Plot empirical support for different hypotheses of species limits -------------------------

as_tibble( bic_best_model_per_group ) %>%
              ggplot( aes( x = seq_along( value ), y = value ) ) +
                scale_x_continuous( breaks = 1:30 ) +
                xlab( "Number of morphological groups" ) +
                ylab( "Empirical support (BIC)" ) +
                geom_point( size = 3, #Add all points
                            shape = 21 ) +
                geom_point( data = as_tibble( h_lack$bic ), # Add hypothesis by Lack
                            shape = 25, 
                            x = length(unique( data_pca$Taxon ) ), 
                            fill = "black", 
                            size = 4) +
                geom_segment( aes( x = length(unique( data_pca$Taxon ) ), y = min( bic_best_model_per_group ), xend = length(unique( data_pca$Taxon ) ), yend = h_lack$bic ), linetype = "dashed", size = 0.3 ) +
                geom_point( data = as_tibble( h_current_taxonomy$bic ), # Add hypothesis of current taxonomy.
                            shape = 15, 
                            x = length( unique( data_pca$New_Taxonomy ) ), 
                            fill = "black",
                            size = 4 ) +
                geom_segment( aes( x = length(unique( data_pca$New_Taxonomy ) ), y = min( bic_best_model_per_group ), xend = length(unique( data_pca$New_Taxonomy ) ), yend = h_current_taxonomy$bic ), linetype = "dashed", size = 0.3 ) +
                geom_point( data = as_tibble( min( bic_best_model_per_group ) ), # Add hypothesis McKay & Zink
                            shape = 24,
                            x = data_pca_varsel_gmm$G - 7,
                            fill = "black",
                            size = 4 ) +
                geom_segment( aes( x = data_pca_varsel_gmm$G - 7, y = min( bic_best_model_per_group ), xend = data_pca_varsel_gmm$G - 7, yend = min( bic_best_model_per_group ) ), linetype = "dashed", size = 0.3 ) +
                geom_point(data = as_tibble( max( bic_best_model_per_group ) ), # Add best mclust model
                           shape = 19,
                           x = data_pca_varsel_gmm$G,
                           fill = "black",
                           size = 4 ) +
                geom_segment( aes( x = data_pca_varsel_gmm$G, y = min( bic_best_model_per_group ), xend = data_pca_varsel_gmm$G, yend = max( bic_best_model_per_group ) ), linetype = "dashed", size = 0.3 ) +
                geom_point( data = as_tibble( data_pca_varsel_gmm_7groups$bic ), # Add equally supported mclust model
                            shape = 19,
                            x = data_pca_varsel_gmm$G - 1,
                            size = 4 ) +
                geom_segment( aes( x = data_pca_varsel_gmm$G - 1, y = min( bic_best_model_per_group ), xend = data_pca_varsel_gmm$G - 1, yend = data_pca_varsel_gmm_7groups$bic ), linetype = "dashed", size = 0.3 ) +
                theme( axis.line = element_line( color = "black", size = 0.5), 
                       axis.title = element_text( size = 15 ), 
                       axis.text = element_text( size = 10 ),
                       panel.border = element_rect( color = "transparent", fill = NA ),
                       panel.grid.minor.y = element_line( size = 0.1, color = "grey" ),
                       panel.background = element_rect( fill = "transparent" ) ) 

# Scatterplots  -------------------------

# Define colors

morphogroups_colors = c( "#999999", 
                          "#E69F00", 
                          "#56B4E9", 
                          "#009E73", 
                          "#F0E442", 
                          "#0072B2", 
                          "#D55E00", 
                          "#CC79A7" ) 

# Scatterplot PC axes -------------------------

data_pca_mclust %>% 
  ggplot( aes( x = PC1, y = PC2, color = factor( mcluster_classification ) ) ) +
    geom_point( size = 3, shape = 21, stroke = 1, alpha = 0.8 ) +
    stat_ellipse( type = "norm", level = 0.95, segments = 1000 ) + # note that this ellipse is calculated directly by ggplot modifying the code in car::ellipse (i.e., not necessarily the same ellipse as in our manuscript. For more details see car::ellipse)
    scale_color_manual( values = morphogroups_colors, name ="Morphological group", position="top" ) + 
    theme( axis.line = element_line( color = "black", size = 0.5 ),
           axis.title = element_text( size = 15 ), 
           axis.text = element_text( size = 14 ),
           panel.border = element_rect( color = "transparent", fill = NA ),
           panel.background = element_rect( fill = "transparent" ),
           legend.key = element_blank( ),
           legend.position = "right" )

# Plot loadings on PC axes --------------------

# Make unitary circle, the expected length of arrows if all traits have equal loading

angle = seq( -pi, pi, length = 50 ) 
circle_df = data.frame( x = sin( angle ), y = cos( angle ) )

# Plot
as_tibble( pca_results$rotation ) %>%
  ggplot( ) +
    geom_segment( aes( x = 0, y = 0, xend = PC1*0.9, yend = PC2*0.9 ), 
                  arrow = arrow( length = unit( 1/2, "picas" ) ), 
                  color = "grey60" ) +
    geom_text( aes( x = PC1, y = PC2,
                    label = c( "Wing length", 
                                "Tail length", 
                                "Bill length", 
                                "Bill depth", 
                                "Bill width", 
                                "Tarsus length" ) ), 
               size = 4.5, 
               check_overlap = TRUE, color = "black", fontface = "bold" ) +
  xlim(-1,1) +
  ylim(-1,1) +
  stat_ellipse( aes( x, y ), 
                data = circle_df, 
                color = "grey80", 
                type = "euclid", 
                level = 0.5 ) +
  xlab( "Loadings PC1" ) + 
  ylab( "Loadings PC2" ) +
  theme( axis.line = element_line( color = "transparent"),
          axis.title = element_text( size = 15 ), 
          axis.text = element_text( size = 13, color = "black" ),
          panel.background = element_blank( ),
          panel.border = element_blank( ),
          panel.grid.major.x = element_line( size = 0.1, color = "grey" ),
          panel.grid.major.y = element_line( size = 0.1, color = "grey" ),
          legend.position = "none" )

# Histogram comparing assignment of specimens to groups between the best Mclust 
# model (8 groups) and alternative taxonomic hypotheses -------------------------

data_pca_mclust  %>% 
  ggplot(aes (mcluster_classification, fill = factor( mcluster_classification ) ) ) + 
    geom_histogram( binwidth = 1 ) + 
    scale_fill_manual( values = morphogroups_colors ) +
    scale_x_continuous( breaks = 1:8 ) +
    facet_grid( ~ New_Taxonomy ) + # Change to Taxon to generate histrogram for Lack's taxonomy
    xlab( "Morphological group" ) +
    ylab( "Specimens" ) +
    theme( axis.line = element_line( color = "black", size = 0.3 ),
           axis.text = element_text( size = 12 ),
           axis.title = element_text( size = 15 ),
           panel.border = element_rect( color = "transparent", fill = NA ),
           panel.grid.minor.y = element_line( size = 0.1, color = "grey" ),
           panel.background = element_rect( color = "grey" , fill = NA ),
           strip.text = element_text( face = "italic" ),
           legend.position = "none" )

# Histogram comparing assignment of specimens to groups between the two equally supported 
# Mclust models (8 and 7 groups) -------------------------

data_pca_mclust  %>% 
  ggplot(aes (mcluster_classification, fill = factor( mcluster_classification ) ) ) + 
  geom_histogram( binwidth = 1 ) + 
  scale_fill_manual( values = morphogroups_colors ) +
  scale_x_continuous( breaks = 1:8 ) +
  facet_grid( ~ mcluster_classification7 ) + # 
  xlab( "Morphological groups" ) +
  ylab( "Specimens" ) +
  theme( axis.line = element_line( color = "black", size = 0.3 ),
         axis.text = element_text( size = 12 ),
         axis.title = element_text( size = 15 ),
         panel.border = element_rect( color = "transparent", fill = NA ),
         panel.grid.minor.y = element_line( size = 0.1, color = "grey" ),
         panel.background = element_rect( color = "grey" , fill = NA ),
         legend.position = "none" )


# Plot results in geographic space - MAP -------------------------

# Summarize geographic distribution in a new object to create ring plots

all_summ = data_pca_mclust %>%
  complete( mcluster_classification = full_seq( mcluster_classification, 1 ), Island ) %>%
  group_by( Island, mcluster_classification ) %>%
  mutate( Z = ifelse(is.na( Taxon ), 0, mcluster_classification ) ) %>%
  tally( Z / mcluster_classification ) %>%
  mutate( fraction = n / sum(n)) %>%
  arrange( Island ) %>%
  mutate( ymax = cumsum( fraction ) ) %>%
  mutate( ymin = c( 0, head( ymax, n=-1 ) ) )

# Create Ring plots per Island 

ring_plots <- list()

for (i in unique(all_summ$Island)) {
	ring_plots[[i]] <- ggplot( data = subset( all_summ, Island == paste(i) ),
                                aes( fill = factor(mcluster_classification), 
                                xmin = 4, 
                                xmax = 7, 
                                ymin = ymin, 
                                ymax = ymax ) ) + 
                    	geom_rect() + 
                    	coord_polar( theta = "y" ) + 
                    	scale_fill_manual( values = morphogroups_colors ) + 
                    	xlim( c( 0, 7 ) )  + 
                    	theme( panel.grid = element_blank(),
                            axis.text = element_blank(), 
                            axis.ticks = element_blank(), 
                            axis.title = element_blank(), 
                            legend.position = "none", 
                            plot.background = element_rect( fill = "transparent", color = NA ), 
                            panel.background = element_rect( fill = "transparent", color = NA ) ) + 
                    	annotate( "text", x = 0, y = 0, label = sum(subset(all_summ, Island==paste(i))$n ), size = 2.8 )
}

# Plot map --------------------------------------

gplot( allgalapagos_NO_NW, maxpixels = 1000000 ) + 
	geom_raster( aes( fill = value ) ) + 
	scale_fill_gradientn( colors = rev( gray.colors(20) ), 
        na.value = "white", 
        limits = c( 1,1800 ), 
        name = "", labels = c( "1 m", "1700 m" ), 
        breaks = c( 1,1800 ) ) + 
	scale_x_continuous( expand = c(0, 0.02) ) + 
	scale_y_continuous( expand = c(0, 0) ) + 
	xlab( "Longitude (degrees)" ) + 
	ylab( "Latitude (degrees)" ) + 
	theme( axis.line.x = element_line( color = "black", size = 0.5), 
        axis.line.y = element_line( color = "black", size = 0.5 ), 
        axis.title = element_text( size = 10 ), 
        axis.text = element_text( size = 10 ), 
        legend.position = c( 0.83, 0.95 ), 
        legend.direction = "horizontal", 
        legend.title = element_text( size = 8 ), 
        legend.text = element_text( size = 8), 
        legend.key.width = unit( 0.25, "cm" ), 
        legend.key.size = unit( 0.3, "cm" ), 
        legend.justification = 'left', 
        panel.background = element_rect( color = "transparent", fill = NA ) ) + 
    	scale_bar(lon = -91.9, 
        lat = -1.95,
        distance_lon = 25, 
        distance_lat = 2, 
        distance_legend = 10, 
        dist_unit = "km", 
        legend_size = 3, 
        legend_colour = "black", 
        rec2_fill = "black", 
        rec_colour = "black", 
        orientation = F ) + 
	annotate( "text", x = -90.38, y = -0.81, label = "Santa Cruz", size = 3 ) + 
	annotate( "text", x = -91.2, y = -1.1, label = "Isabela", size = 3 ) + 
	annotate( "text", x = -91.5, y = -0.5, label = "Fernandina", size = 3 ) + 
	annotate( "text", x = -90.1, y = -0.29, label = "Baltra", size = 3 ) + 
	annotate( "segment", x = -90.265, xend = -90.16, y = -0.46, yend = -0.41, color = "black", size = 0.2 ) + 
	annotate( "text", x = -90.36, y = -0.12, label = "Daphne", size = 3 ) + 
	annotate( "segment", x = -90.37, xend = -90.37, y = -0.42, yend = -0.34, color = "black", size = 0.2 ) + 
	annotate( "text", x = -89.7, y = -1.47, label = "Española", size = 3 ) + 
	annotate( "text", x = -90.45, y = -1.4, label = "Floreana", size = 3 ) + 
	annotate( "text", x = -89.525, y = -1.16, label = "Gardner", size = 3 ) + 
	annotate( "segment", x = -89.64, xend = -89.58, y = -1.35, yend = -1.32, color = "black", size = 0.2 ) + 
	annotate( "text", x = -89.95, y = 0.58, label = "Genovesa", size = 3 ) + 
	annotate( "text", x = -90.45, y = 0.6, label = "Marchena", size = 3 ) + 
	annotate( "text", x = -90.75, y = 0.87, label = "Pinta", size = 3 ) + 
	annotate( "text", x = -90.75, y = -0.94, label = "Pinzón", size = 3 ) + 
	annotate( "segment", x = -90.665, xend = -90.7, y = -0.6, yend = -0.73, color = "black", size = 0.2 ) + 
	annotate( "text", x = -90.81, y = -0.56, label = "Rábida", size = 3 ) + 
	annotate( "segment", x = -90.7, xend = -90.76, y = -0.418, yend = -0.43, color = "black", size = 0.2 ) + 
	annotate( "text", x = -89.3, y = -0.45, label = "San Cristóbal", size  =3 ) + 
	annotate( "text", x = -90.05, y = -0.88, label = "Santa Fé", size = 3 ) + 
	annotate( "text", x = -90.7, y = 0.05, label = "Santiago", size = 3 ) + 
	annotate( "rect", xmin = -91.975, xmax = -91.5, ymin = 0.3, ymax = 0.99, color= "black", fill="transparent", size = 0.2 ) + 
	annotate( "text", x = -91.8, y = 0.94, label = "Darwin", size = 3 ) + 
    annotate( "text", x = -91.6, y = 0.55, label = "Wolf", size = 3 )
    
# Plot rings per island. Note this assumes a new set of coordinates in the plotting canvas, 
# so if there are changes in the axis or margins in map above, the x y positions of each 
# ring must be changed here. --------------------------------------

# Baltra
baltra_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.65, 
                 y = 0.56 )
print( ring_plots$Baltra, vp = baltra_vp )                

# Daphne
daphne_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.57, 
                 y = 0.605 )
print( ring_plots$Daphne, vp = daphne_vp )  

# Espanola
espanola_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.765, 
                 y = 0.205 )
print( ring_plots$Espanola, vp = espanola_vp )

# Fernandina
fernandia_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.23, 
                	 y = 0.495 )
print( ring_plots$Fernandina, vp = fernandia_vp )

# Floreana
floreana_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.55, 
                 y = 0.227 )
print( ring_plots$Floreana, vp = floreana_vp )

# Gardner
gardner_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.826, 
                 y = 0.30 )
print( ring_plots$Gardner, vp = gardner_vp )

# Genovesa
genovesa_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.695, 
                 y = 0.815 )
print( ring_plots$Genovesa, vp = genovesa_vp )

# Isabela                
isabela_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.32, 
                 y = 0.315 )
print( ring_plots$Isabela, vp = isabela_vp )

# Marchena
marchena_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.54, 
                 y = 0.825 )
print( ring_plots$Marchena, vp = marchena_vp )

# Pinta
pinta_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.455, 
                 y = 0.9 )
print( ring_plots$Pinta, vp = pinta_vp )

# Pinzon
pinzon_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.465, 
                 y = 0.44 )
print( ring_plots$Pinzon, vp = pinzon_vp )

# Rabida
rabida_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.432, 
                 y = 0.55 )
print( ring_plots$Rabida, vp = rabida_vp )

# San Cristobal
sancristobal_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.88, 
                 y = 0.51 )
print( ring_plots$"San Cristobal", vp = sancristobal_vp )

# Santa Cruz
santacruz_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.57, 
                 y = 0.405 )
print( ring_plots$"Santa Cruz", vp = santacruz_vp )

# Santa Fe
santafe_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.665, 
                 y = 0.38 )
print( ring_plots$"Santa Fe", vp = santafe_vp )

# Santiago
santiago_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.48, 
                 y = 0.66 )
print( ring_plots$Santiago, vp = santiago_vp )

# NW ISLANDS
# Darwin
darwin_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.13, 
                 y = 0.925 )
print( ring_plots$Darwin, vp = darwin_vp )  

# Wolf
wolf_vp <- viewport( width = 0.1, 
                 height = 0.1, 
                 x = 0.2, 
                 y = 0.81 )
print( ring_plots$Wolf, vp = wolf_vp )  


####################
## DATA SUMMARIES ##
####################

# How many islands per morphological group? -------------------------

data_pca_mclust %>% 
  group_by( mcluster_classification ) %>% 
  summarise( n_distinct( Island ) ) #%>% # Uncomment here and below to estimate median
  #summarise(median(`n_distinct(Island)`))

# How many morphological groups per island? -------------------------

data_pca_mclust %>% 
  group_by( Island ) %>% 
  summarise(n_distinct( mcluster_classification ), n_distinct( `Taxon-Island` ) )

# How many specimens for each species sensu Lack are assigned to 
# different morphological groups? -------------------------

data_pca_mclust %>% 
  group_by( Taxon,  mcluster_classification ) %>% 
  tally() %>%
  mutate( fraction = n / sum(n) )

# How many specimens for each species sensu current taxonomy are 
# assigned to different morphological groups? -------------------------

data_pca_mclust %>% 
  group_by( New_Taxonomy,  mcluster_classification ) %>% 
  tally() %>%
  mutate( fraction = n / sum(n) )

####################
## DATA SUMMARIES ##
#################### 

# Univariate plots of orignal data -------------------------
data %>%
  ggplot(aes( Tail ) ) + # Change trait to be plotted
    geom_histogram( binwidth = 1 ) + # May need to adjust binwidth if traits changes
    facet_grid( ~ Taxon ) +
    labs( y = "Specimens" ) +
    theme( axis.line = element_line( color = "black", size = 0.3 ),
           axis.text = element_text( size = 14 ),
           axis.title = element_text( size = 15 ),
           panel.border = element_rect( color = "transparent", fill = NA ),
           panel.grid.minor.y = element_line( size = 0.1, color = "grey" ),
           panel.background = element_rect( color = "grey" , fill = NA ),
           strip.text = element_text( face = "italic" ) ) 

# Bivariate plots of original data -------------------------
data %>%
  ggplot(aes( x = Wing, y = Tail ) ) + # Change to traits of interest
    geom_point( ) + 
    facet_wrap( ~ Taxon ) +
    theme_bw() +
    theme( axis.text = element_text( size = 14 ),
           axis.title = element_text( size = 15 ),
           strip.text = element_text( face = "italic") )

data %>%
  ggplot( aes( x = Wing, y = Tail, color = Taxon ) ) + # Change to traits of interest
    geom_point( ) +
    theme_bw() +
    theme( axis.text = element_text( size = 14 ),
           axis.title = element_text( size = 15 ) )

# All pairs of traits of original data -------------------------
data %>%  
  ggpairs(aes( color = Taxon, alpha = 0.7 ),
         columns = c( "Wing", "Tail", "Blength", "Bdepth", "Bwidth", "Tarsus" ),
         upper = list( continuous = "density" ),
         lower = list( continuous = "points" ),
         axisLabels = "show"
    ) +
    theme_bw( )

# Variable Selection using all PC Axes and Forward algorithm -------------------------

data_pca_varselfwd =
  data_pca %>%
  dplyr::select( PC1:PC6 ) %>%
  clustvarsel( G = 1:30, search = c( "greedy" ), direction = c( "forward" ) )

#data_pca_varselfwd$subset # Uncomment to see resulting subset 
