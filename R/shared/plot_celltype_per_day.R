library(Seurat)
library(dplyr)
library(ggplot2)

source('R/shared/filenames.R')

plot_celltype_per_day = function (sobj) {
  # Extract metadata including Day and Cell_Type
  metadata_df = as.data.frame(sobj@meta.data)
  metadata_df$Day <- factor(metadata_df$Day, levels = rev(unique(metadata_df$Day)))
  metadata_df <- dplyr::count(metadata_df, Day, Cell_Type)

  # Plot using ggplot2
  plot <- ggplot(metadata_df, aes(x = Day, y = n, fill = Cell_Type)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = custom_palette) +
    #scale_fill_viridis_d(option = "D") +
    labs(x="Day", y="Cell Count", fill="Cell Type", title='Cell type counts for each day') +
    coord_flip() +  # Flips the x and y axes
    theme(legend.position = "bottom")
  
  # Save the plot
  ggsave(figure_filename(sobj@project.name, 'celltypes', "cell_type_composition_by_day"), plot = plot, width = 10, height = 8, dpi = 300)
}
