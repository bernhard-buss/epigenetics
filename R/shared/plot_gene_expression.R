library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

plot_gene_expression = function(sobj, genes, topic_suffix, celltype) {
  topic = 'marker_expression'
  
  # Extract expression data (assuming this has already been done)
  sobj_subset = subset(sobj, subset = Cell_Type == celltype)
  Idents(sobj_subset) <- 'Day'

  sobj_subset <- AggregateExpression(sobj_subset, features = genes, return.seurat = TRUE)
  avg_expression = GetAssayData(sobj_subset)
  avg_expression_df = as.data.frame(avg_expression)
  
  # Add 'Gene' information
  avg_expression_df$Gene <- row.names(avg_expression_df)
  
  # Convert to long format for ggplot
  avg_expression_long <- pivot_longer(avg_expression_df, cols = -Gene, names_to = "Day", values_to = "Expression")

  # Ensure 'Day' is treated as a factor to maintain the order in the plot
  avg_expression_long$Day <- factor(avg_expression_long$Day, levels = unique(avg_expression_long$Day))
  
  # Plot all genes together
  p <- ggplot(avg_expression_long, aes(x = Day, y = Expression, color = Gene, group = Gene)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(title = paste0(celltype, ' ', topic_suffix, ": Expression of Genes Over Days"), x = "Day", y = "Average Expression")
    scale_color_viridis_d() # Use a color palette that provides distinct colors for each gene
  
  # Save the plot
  ggsave(filename = figure_filename(sobj@project.name, topic, paste0(celltype, ' ', topic_suffix, "_Gene_Expression_Over_Days")), plot = p, width = 10, height = 8, dpi = 300)
}
