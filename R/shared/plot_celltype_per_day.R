library(Seurat)
library(dplyr)
library(ggplot2)

source('R/shared/filenames.R')

plot_celltype_per_day = function (sobj) {
  # Extract metadata including time points and cell types
  metadata_df <- sobj@meta.data %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell") %>%
    mutate(cell_type = as.character(Idents(sobj))) %>%
    count(time_point, cell_type)  # Replace 'time_point' with the actual column name
  
  # Plot using ggplot2
  plot <- ggplot(metadata_df, aes(x = time_point, y = n, fill = cell_type)) + 
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Time Point", y = "Cell Count", fill = "Cell Type") +
    coord_flip() +  # Flips the x and y axes
    theme(legend.position = "bottom")
  
  # View the plot
  print(plot)
  
  # Save the plot
  ggsave("cell_type_composition_by_day.png", plot = plot, width = 10, height = 8, dpi = 300)
}