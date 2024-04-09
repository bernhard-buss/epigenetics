
plot_DGE_heatmap <- function(sobj, celltype, relevant_chromatin_genes) {
  #sobj1_celltype <- subset(sobj1, subset = Cell_Type_heuristic == 'astrocytes')
  sobj_celltype <- subset(sobj, subset = Cell_Type == celltype)
  Idents(sobj_celltype) <- 'Day'
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
  DGE_per_day_heatmap = DoHeatmap(sobj, features = filtered_chromatin_markers_all$gene, group.by = 'Day', size = 3, angle = 90) +
    ggtitle(paste0(celltype, ' DGE per day')) +
    theme(plot.title = element_text(hjust = 0.5))  # This centers the title

  ggsave(figure_filename(sobj@project.name, paste0('DGE/', celltype, '_DGE_per_day_heatmap')), plot = DGE_per_day_heatmap, width = 10, height = 8, dpi = 300)
}
