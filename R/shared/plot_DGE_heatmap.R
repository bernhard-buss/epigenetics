
perform_DGE_Two_Days <- function(sobj, celltype, genes, day1, day2) {
  topic <- paste0(day1, '_', day2)
  
  print(paste0('Performing DGE for ', celltype, ', ', day1, ' vs ', day2))
  
  sobj_celltype <- subset(sobj, subset = Cell_Type == celltype)
  Idents(sobj_celltype) <- 'Day'
  
  # Validate the Day clusters
  # Check if both day1 and day2 are present in the metadata column 'Day'
  all_days_present <- all(c(day1, day2) %in% sobj_celltype@meta.data$Day)
  # If not all days are present, print a message
  if (!all_days_present) {
    missing_days <- setdiff(c(day1, day2), unique(sobj_celltype@meta.data$Day))
    cat("The following days are missing in the Seurat object:", paste(missing_days, collapse = ", "), "\n")
    return()
  }
  # Count the number of cells in each cluster
  cluster_counts <- table(sobj_celltype@meta.data$seurat_clusters)
  # Find clusters with fewer than 3 cells
  clusters_with_few_cells <- intersect(cluster_counts[cluster_counts < 3], c(day1, day2))
  # Check if there are any clusters with fewer than 3 cells
  if (length(clusters_with_few_cells) > 0) {
    cat("The following days have fewer than 3 cells:\n")
    print(clusters_with_few_cells)
    return()
  }
  
  sobj_celltype <- subset(sobj_celltype, subset = Day %in% c(day1, day2))
  
  filtered_chromatin_markers_per_day <- FindMarkers(
    object = sobj_celltype,
    features = genes,
    ident.1 = day1,
    ident.2 = day2,
    min.pct = 0.25, # Gene must be detected in at least 25% of cells within a cluster
    #min.diff.pct = 0.5,
    logfc.threshold = 0.25, # Minimum log-fold change
    test.use = 'wilcox', # Use Wilcoxon Rank Sum test
    #p.adjust.method = 'bonferroni' # Bonferroni correction for multiple testing
  )
  
  print('Generating volcano plot')
  volcano_plot <- EnhancedVolcano(
    filtered_chromatin_markers_per_day,
    x='avg_log2FC', y='p_val',
    lab='gene',
    #pCutoff = 0.05
  ) + ggtitle(paste0(celltype, ' DGE ', day1, ' vs ', day2)) +
    theme(plot.title = element_text(hjust = 0.5))  # This centers the title
  ggsave(figure_filename(sobj@project.name, topic, paste0(celltype, '_DGE_', topic, '_volcano')), plot = volcano_plot, width = 10, height = 8, dpi = 300)
  
  print('Generating significant marker list table')
  filtered_chromatin_markers_per_day_significant = filtered_chromatin_markers_per_day[filtered_chromatin_markers_per_day$p_val < 0.05,]
  write.csv(filtered_chromatin_markers_per_day_significant, file = table_filename(sobj1@project.name, topic, paste0(celltype, '_', topic, '_DGE_markers')))
  
  print('Generating heatmap')
  DGE_per_day_heatmap = DoHeatmap(sobj_celltype, features = filtered_chromatin_markers_per_day_significant$gene, group.by = 'Day', size = 3, angle = 90) +
    ggtitle(paste0(celltype, ' DGE ', day1, ' vs ', day2)) +
    theme(plot.title = element_text(hjust = 0.5))  # This centers the title
  ggsave(figure_filename(sobj@project.name, topic, paste0(celltype, '_DGE_', topic, '_heatmap')), plot = DGE_per_day_heatmap, width = 10, height = 8, dpi = 300)
}


plot_DGE_heatmap <- function(sobj, celltype, relevant_chromatin_genes) {
  #sobj1_celltype <- subset(sobj1, subset = Cell_Type_heuristic == 'astrocytes')
  sobj_celltype <- subset(sobj, subset = Cell_Type == celltype)
  Idents(sobj_celltype) <- 'Day'
  
  volcano_plot <- VariableFeaturePlot(sobj) +
    ggtitle(paste0(celltype, ' Variable features')) +
    theme(plot.title = element_text(hjust = 0.5))  # This centers the title
  ggsave(figure_filename(sobj@project.name, 'DGE', paste0(celltype, '_DGE_per_day_heatmap')), plot = DGE_per_day_heatmap, width = 10, height = 8, dpi = 300)

  
  filtered_chromatin_markers_all = FindAllMarkers(
    object = sobj_celltype,
    features = relevant_chromatin_genes,
    #only.pos = TRUE, # Consider only positive markers
    min.pct = 0.25, # Gene must be detected in at least 25% of cells within a cluster
    #min.diff.pct = 0.5,
    logfc.threshold = 0.25, # Minimum log-fold change
    test.use = 'wilcox', # Use Wilcoxon Rank Sum test
    p.adjust.method = 'bonferroni' # Bonferroni correction for multiple testing
  )
  
  # Heatmap for celltype per day
  DGE_per_day_heatmap = DoHeatmap(sobj_celltype, features = filtered_chromatin_markers_all$gene, group.by = 'Day', size = 3, angle = 90) +
    ggtitle(paste0(celltype, ' DGE per day')) +
    theme(plot.title = element_text(hjust = 0.5))  # This centers the title

  ggsave(figure_filename(sobj@project.name, 'DGE', paste0(celltype, '_DGE_per_day_heatmap')), plot = DGE_per_day_heatmap, width = 10, height = 8, dpi = 300)
}
