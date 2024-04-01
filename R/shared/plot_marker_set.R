library(Seurat)
library(ggplot2)

plot_marker_set <- function(sobj, markers, marker_set_name, ncol = 2, dpi = 72, pixelsPerRow = 300) {
  # Calculate the height dynamically based on the number of markers
  num_features <- length(markers)
  num_rows <- ceiling(num_features / 2)
  plot_height_px <- num_rows * pixelsPerRow
  
  # Create and save the FeaturePlot
  f_plot <- FeaturePlot(sobj, features = markers, ncol = ncol, order = TRUE) + NoLegend()
  feature_png_filename <- paste0('./plots/', marker_set_name, "_feature_plot.png")
  ggsave(feature_png_filename, plot = f_plot, width = 800/dpi, height = plot_height_px/dpi, dpi = dpi)
  
  # Create and save the VlnPlot
  vln_plot <- VlnPlot(sobj, features = markers, ncol = ncol)
  violin_png_filename <- paste0('./plots/', marker_set_name, "_vln_plot.png")
  ggsave(violin_png_filename, plot = vln_plot, width = 800/dpi, height = plot_height_px/dpi, dpi = dpi)
}
