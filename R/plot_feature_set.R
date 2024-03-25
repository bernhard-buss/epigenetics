library(Seurat)
library(ggplot2)

plot_feature_set <- function(sobj, features, feature_set_name, ncol = 2, dpi = 72, pixelsPerRow = 300, pixelsPerCol = 300) {
  # Calculate the height dynamically based on the number of features
  num_features <- length(features)
  num_rows <- ceiling(num_features / ncol)
  plot_height_px <- num_rows * pixelsPerRow
  plot_width_px <- ncol * pixelsPerCol
  
  # Create and save the FeaturePlot
  f_plot <- FeaturePlot(sobj, features = features, ncol = ncol, order = TRUE) + NoLegend()
  feature_png_filename <- paste0('./plots/', feature_set_name, "_feature_plot.png")
  ggsave(feature_png_filename, plot = f_plot, width = plot_width_px/dpi, height = plot_height_px/dpi, dpi = dpi)
  
  # Create and save the VlnPlot
  vln_plot <- VlnPlot(sobj, features = features, ncol = ncol)
  violin_png_filename <- paste0('./plots/', feature_set_name, "_vln_plot.png")
  ggsave(violin_png_filename, plot = vln_plot, width = plot_width_px/dpi, height = plot_height_px/dpi, dpi = dpi)
}
