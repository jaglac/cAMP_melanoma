#library(PCAtools, ggplot2, viridis, ggsci, gginnards)

#make pca object with PCAtools
#add gene expression values or pathway mean z scores to pca object metadata first
#e.g. pca$metadata$MITF <- unlist(vst_df['MITF',2:ncol(vst_df)])


biplot_gene <- function(pca_object, x, y, gene, ellipse){
  #requires PCAtools, ggplot2, viridis, ggsci, gginnards
  #use previously created PCAtools pca object
  plot <- PCAtools::biplot(pca_object,
                 x = x, y = y,
                 colby = gene,
                 shape = ellipse,
                 legendPosition = 'right',
                 hline = 0, vline = 0, 
                 lab = NULL, gridlines.major = FALSE, gridlines.minor = FALSE) + 
    viridis::scale_color_viridis(option="viridis", direction = -1) +
    ggplot2::stat_ellipse(geom="polygon", aes(fill=shape), alpha = 0.25) + 
    ggsci::scale_fill_jama() + 
    ggplot2::guides(fill = guide_legend(title = "differentiation state", order=1),
           color = guide_colorbar(title = paste0(gene, " expression"), order=2), 
           shape = "none") +
    ggplot2::ggtitle(paste0(gene, " expression")) +
    #written to build subtype ellipses based on biplot function's assigned shapes, with 4 subtype differentiation states
    #would need to edit scale_shape_manual if more than 4 subtypes/ellipses categories or if you want the actual shape to change
    ggplot2::scale_shape_manual(values = c(16,16,16,16,16))
  
  #reorder plot layers so that points are on the top 
  plot <- gginnards::move_layers(plot, "GeomPolygon", position="bottom")
  plot <- gginnards::move_layers(plot, match_type = "GeomHline", position = "bottom")
  plot <- gginnards::move_layers(plot, match_type = "GeomHline", position = "bottom")
  
  #adjust plot size to fit ellipses for tsoi et all 53-cell line data
  #possibly have to change the multiplication factor for other datasets
  #if stat_ellipse polygons have missing chunks/look funny, most likely caused by axis limits being too small to fit the ellipse shape
  pb <- ggplot2::ggplot_build(plot)
  x_max <- max(pb[["data"]][[2]]$x, na.rm=T) * 1.1
  x_min <- min(pb[["data"]][[2]]$x, na.rm=T) * 1.1
  y_max <- max(pb[["data"]][[2]]$y, na.rm=T) * 1.25
  y_min <- min(pb[["data"]][[2]]$y, na.rm=T) * 1.25
  plot <- plot + ggplot2::ylim(c(y_min, y_max)) + ggplot2::xlim(c(x_min, x_max)) 
  
  plot
}


######################## pathway mean z version

        #requires PCAtools, ggplot2, viridis, ggsci, gginnards
        #use previously created PCAtools pca object with mean z dataframe appended to pca$metadata df

biplot_meanZ <- function(pca_object, x, y, geneset, ellipse){
  plot <- PCAtools::biplot(pca_object,
                           x = x, y = y,
                           colby = geneset,
                           shape = ellipse,
                           legendPosition = 'right',
                           hline = 0, vline = 0, 
                           lab = NULL, gridlines.major = FALSE, gridlines.minor = FALSE) + 
    viridis::scale_color_viridis(option="viridis", direction = -1) +
    ggplot2::stat_ellipse(geom="polygon", aes(fill=shape), alpha = 0.25) + 
    ggsci::scale_fill_jama() + 
    ggplot2::guides(fill = guide_legend(title = "differentiation state", order=1),
                    color = guide_colorbar(title = paste0(geneset, "\n", "mean z score"), order=2), 
                    shape = "none") +
    ggplot2::ggtitle(paste0(geneset, " mean z score")) +
    #written to build subtype ellipses based on biplot function's assigned shapes, with 4 subtype differentiation states
    #would need to edit scale_shape_manual if more than 4 subtypes/ellipses categories or if you want the actual shape to change
    ggplot2::scale_shape_manual(values = c(16,16,16,16,16))
  
  #reorder plot layers so that points are on the top 
  plot <- gginnards::move_layers(plot, "GeomPolygon", position="bottom")
  plot <- gginnards::move_layers(plot, match_type = "GeomHline", position = "bottom")
  plot <- gginnards::move_layers(plot, match_type = "GeomHline", position = "bottom")
  
  #adjust plot size to fit ellipses for tsoi et all 53-cell line data
  #possibly have to change the multiplication factor for other datasets
  #if stat_ellipse polygons have missing chunks/look funny, most likely caused by axis limits being too small to fit the ellipse shape
  pb <- ggplot2::ggplot_build(plot)
  x_max <- max(pb[["data"]][[2]]$x, na.rm=T) * 1.1
  x_min <- min(pb[["data"]][[2]]$x, na.rm=T) * 1.1
  y_max <- max(pb[["data"]][[2]]$y, na.rm=T) * 1.25
  y_min <- min(pb[["data"]][[2]]$y, na.rm=T) * 1.25
  plot <- plot + ggplot2::ylim(c(y_min, y_max)) + ggplot2::xlim(c(x_min, x_max)) 
  
  plot
}


############################ gene Z score version

#make pca object with PCAtools
#add gene expression values or pathway mean z scores to pca object metadata first
# e.g. pca$metadata$MITF_z <- unlist(vst_df_z['MITF',2:ncol(vst_df_z)])


biplot_geneZ <- function(pca_object, x, y, gene, ellipse){
  #requires PCAtools, ggplot2, viridis, ggsci, gginnards
  #use previously created PCAtools pca object
  plot <- PCAtools::biplot(pca_object,
                           x = x, y = y,
                           colby = gene,
                           shape = ellipse,
                           legendPosition = 'right',
                           hline = 0, vline = 0, 
                           lab = NULL, gridlines.major = FALSE, gridlines.minor = FALSE) + 
    viridis::scale_color_viridis(option="viridis", direction = -1) +
    ggplot2::stat_ellipse(geom="polygon", aes(fill=shape), alpha = 0.25) + 
    ggsci::scale_fill_jama() + 
    ggplot2::guides(fill = guide_legend(title = "differentiation state", order=1),
                    color = guide_colorbar(title = paste0(gene, " z score"), order=2), 
                    shape = "none") +
    ggplot2::ggtitle(paste0(gene, " z score")) +
    #written to build subtype ellipses based on biplot function's assigned shapes, with 4 subtype differentiation states
    #would need to edit scale_shape_manual if more than 4 subtypes/ellipses categories or if you want the actual shape to change
    ggplot2::scale_shape_manual(values = c(16,16,16,16,16))
  
  #reorder plot layers so that points are on the top 
  plot <- gginnards::move_layers(plot, "GeomPolygon", position="bottom")
  plot <- gginnards::move_layers(plot, match_type = "GeomHline", position = "bottom")
  plot <- gginnards::move_layers(plot, match_type = "GeomHline", position = "bottom")
  
  #adjust plot size to fit ellipses for tsoi et all 53-cell line data
  #possibly have to change the multiplication factor for other datasets
  #if stat_ellipse polygons have missing chunks/look funny, most likely caused by axis limits being too small to fit the ellipse shape
  pb <- ggplot2::ggplot_build(plot)
  x_max <- max(pb[["data"]][[2]]$x, na.rm=T) * 1.1
  x_min <- min(pb[["data"]][[2]]$x, na.rm=T) * 1.1
  y_max <- max(pb[["data"]][[2]]$y, na.rm=T) * 1.25
  y_min <- min(pb[["data"]][[2]]$y, na.rm=T) * 1.25
  plot <- plot + ggplot2::ylim(c(y_min, y_max)) + ggplot2::xlim(c(x_min, x_max)) 
  
  plot
}
