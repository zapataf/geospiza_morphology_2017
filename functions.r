pairwise_scatterplot = function(df, mapping, ...) { 

  p <- ggplot( data = df, mapping = mapping ) +
        geom_point(size = 3, shape = 21, stroke = 1, alpha = 0.8 ) +
        stat_ellipse( type = "norm", level = 0.95, segments = 1000 ) +
        scale_color_manual( values = morphogroups_colors ) + 
        theme(axis.line = element_line( color = "black", size = 0.5),
            axis.title = element_text( size = 15 ), 
            axis.text = element_text( size = 14 ),
            panel.border = element_rect( color = "transparent", fill = NA ),
            panel.background = element_rect( fill = "transparent" ),
            legend.position = "none" )
  p
}


for (i in 1:ncol(combn(test, 2))){print(combn(test, 2)[,i])}