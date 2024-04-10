#install.packages('Seurat')
BiocManager::install('EnhancedVolcano')
library('Seurat')
library('ggplot2')
library(dplyr)
library('EnhancedVolcano')

source('R/shared/index.R')
source('R/paper1/read_paper1_data.R')
source('R/paper1/celltypes_paper1.R')


# SEURAT PIPELINE

# Create Seurat object with min genes = 500
sobj1_full <- setup_seurat(
  project='paper1',
  counts=read_paper1_data(),
  metadata=read_paper1_metadata()
)

# add custom metadata
sobj1_full <- map_metadata_column(sobj1_full, source = 'orig_ident', target = 'Day', lambda = extract_prefix)
sobj1_full <- map_metadata_column(sobj1_full, source = 'New_cellType', target = 'Cell_Type', lambda = function(x) return(x))

# Filter mito < 7.5 % && genes > 500 (& Cell_Type in celltypes)
#sobj1 <- subset(sobj1_full, subset = percent_mito < 0.075 & nFeature_RNA > 500) # includes DL (deleterious) and low quality
#sobj <- subset(sobj, subset = percent_mito < 0.075 & nFeature_RNA > 500 & Cell_Type %in% celltypes)
#sobj1 <- subset(sobj1_full, subset = percent_mito < 0.075 & nFeature_RNA > 500 & Cell_Type %in% celltypes_paper1 & Phase == 'G1')
days = c("E16", "E17", "E18", "P1", "P4")
#unique(FetchData(sobj1, 'Day')$Day)
sobj1 <- subset(sobj1_full, subset = percent_mito < 0.075 & nFeature_RNA > 500 & Day %in% days)

sobj1 <- cluster_seurat(sobj1, nfeatures=3000, resolution=1)
sobj1 <- add_celltype_metadata(sobj1, celltype_marker_lists=celltype_marker_lists)


# ANALYSIS
# Plotting & Analysis

DimPlot(sobj1, reduction = "umap", group.by = 'seurat_clusters_orig')
DimPlot(sobj1, reduction = "umap", group.by = 'seurat_clusters', label = TRUE)

origIdentPlot <- DimPlot(sobj1, reduction = "umap", group.by = 'Day')
origIdentPlot

# plot cell types over all days
UMAPPlot(sobj1, group.by = 'Cell_Type', label=TRUE, repel = TRUE)
# plot cell types for each day
cellTypePlot = UMAPPlot(sobj1, group.by = 'Cell_Type', split.by = 'Day', ncol = 3, pt.size = 0.5)
cellTypePlot
ggsave(figure_filename(sobj1@project.name, 'Cell_Type_per_day'), plot = cellTypePlot, width = 10, height = 8, dpi = 300)


# Loop over the list, create a FeaturePlot for each set of markers, and save as PNG
for (celltype in names(celltype_marker_lists)) {
  plot_feature_set(sobj1, celltype_marker_lists[[celltype]], celltype)
}

# plot all celltype markers score in one PNG file
plot_feature_set(sobj1, celltype_markers_score_norm_features(celltype_marker_lists), 'celltype_markers_scores', ncol = 6)

# plot cell types heuristic over all days
celltypeHeuristicCombinedPlot <- UMAPPlot(sobj1, group.by = 'Cell_Type_heuristic', label=TRUE, repel = TRUE)
ggsave(figure_filename(sobj1@project.name, 'Cell_Type_heuristic_combined'), plot=celltypeHeuristicCombinedPlot, width = 10, height = 8, dpi = 300)
# plot cell types heuristic for each day
celltypeHeuristicPlot <- UMAPPlot(sobj1, group.by = 'Cell_Type_heuristic', split.by = 'Day', ncol = 3, pt.size = 0.5)
ggsave(figure_filename(sobj1@project.name, 'Cell_Type_heuristic_per_day'), plot = celltypeHeuristicPlot, width = 10, height = 8, dpi = 300)


# ANALYSIS: Find all markers and plot top markers for each cluster
all_markers <- FindAllMarkers(
  object = sobj1,
  only.pos = TRUE, # Consider only positive markers
  min.pct = 0.50, #0.25, # Gene must be detected in at least 25% of cells within a cluster
  #min.diff.pct = 0.5,
  logfc.threshold = 0.7, #0.25, # Minimum log-fold change
  test.use = 'wilcox', # Use Wilcoxon Rank Sum test
  p.adjust.method = 'bonferroni' # Bonferroni correction for multiple testing
)

# Filter markers based on adjusted P-value < 0.05
significant_markers <- all_markers[all_markers$p_val_adj < 0.05, ]

# Sort the markers within each cluster by adjusted p-value and log-fold change
significant_markers_sorted <- significant_markers %>%
  arrange(cluster, desc(avg_log2FC))

# Extract the top 8 markers for each cluster
top_markers_per_cluster <- significant_markers_sorted %>%
  group_by(cluster) %>%
  slice_head(n = 8)
# Optionally, sort them back by cluster for easier analysis
top_markers_per_cluster <- top_markers_per_cluster %>%
  ungroup() %>%  # Ensure the data is ungrouped for sorting
  arrange(as.numeric(cluster))

cluster_marker_lists <- split(top_markers_per_cluster$gene, top_markers_per_cluster$cluster)

for (cluster in names(cluster_marker_lists)) {
  markers <- cluster_marker_lists[[cluster]]
  if (length(markers) > 0) {
    marker_set_name <- paste("Cluster", cluster, "TopMarkers")
    # Call the plot_marker_set function for the current set of markers
    plot_marker_set(sobj, markers, marker_set_name)
  }
}



# INTERPRETATION: gather celltype assignment lists
celltype_assignments = assign_cell_types(cluster_marker_lists, celltype_marker_lists)
celltype_assignments
print(celltype_assignments)

all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
top10$gene
DoHeatmap(sobj1, features = top10$gene) + NoLegend()
celltype_marker_lists['cpnLayer23']
DoHeatmap(sobj1, features = celltype_marker_lists['cpnLayer56']) + NoLegend()



# ANALYSIS: DGE for epigenetic modifiers/chromatin genes
# for all celltypes
relevant_chromatin_genes = intersect(as.vector(genes_chromatin$gene), Features(sobj1))
for (celltype in unique(sobj1@meta.data$Cell_Type)) {
  perform_DGE_Two_Days(sobj1, celltype, genes=relevant_chromatin_genes, day1='E18', day2='P1')
}
# Manual single celltype
# visualise the variable features being used
celltype = 'CThPN'
day1 = 'E18'
day2 = 'P1'
sobj_celltype <- subset(sobj1, subset = Cell_Type == celltype)
Idents(sobj_celltype) <- 'Day'
sobj_celltype <- subset(sobj_celltype, subset = Day %in% c(day1, day2))

relevant_chromatin_genes = intersect(as.vector(genes_chromatin$gene), Features(sobj1))

filtered_chromatin_markers_days = FindMarkers(
  object = sobj_celltype_day,
  features = relevant_chromatin_genes,
  ident.1 = day1,
  ident.2 = day2,
  min.pct = 0.25, # Gene must be detected in at least 25% of cells within a cluster
  #min.diff.pct = 0.5,
  logfc.threshold = 0.25, # Minimum log-fold change
  test.use = 'wilcox', # Use Wilcoxon Rank Sum test
  #p.adjust.method = 'bonferroni' # Bonferroni correction for multiple testing
)

volcano_plot <- EnhancedVolcano(
  filtered_chromatin_markers_days,
  x='avg_log2FC', y='p_val',
  lab=filtered_chromatin_markers_days$gene,
  #pCutoff = 0.05
) + ggtitle(paste0(celltype, ' DGE ', day1, ' vs ', day2)) +
  theme(plot.title = element_text(hjust = 0.5))  # This centers the title
#ggsave(figure_filename(sobj@project.name, 'DGE', paste0(celltype, '_DGE_per_day_volcano')), plot = volcano_plot, width = 10, height = 8, dpi = 300)

filtered_chromatin_markers_days_significant = filtered_chromatin_markers_days[filtered_chromatin_markers_days$p_val < 0.05,]
#write.csv(filtered_chromatin_markers_E18_P1_significant, file = table_filename(sobj1@project.name, 'E18_P1', paste0(celltype, '_E18_P1_DGE_markers')))

# Heatmaps for epigenetic modifiers
relevant_chromatin_genes = intersect(as.vector(genes_chromatin$gene), Features(sobj1))
for (celltype in unique(sobj1@meta.data$Cell_Type)) {
  plot_DGE_heatmap(sobj1, celltype, relevant_chromatin_genes)
}

dim(filtered_chromatin_markers_all)
filtered_chromatin_markers_all_ordered <- filtered_chromatin_markers_all[order(filtered_chromatin_markers_all$cluster),]
chromatin_markers_astrocytes = FindMarkers(
  object = sobj1,
  features = relevant_chromatin_genes,
  ident.1='Astrocytes',
  group.by='New_cellType',
  logfc.threshold = 1
)
rownames(chromatin_markers_astrocytes)

DoHeatmap(sobj1, features = epigenetic_modifiers) + NoLegend()

DoHeatmap(sobj1, features = epigenetic_modifiers, group.by = 'New_cellType', size = 3, angle = 90)

DoHeatmap(sobj1, features = epigenetic_modifiers, group.by = 'Cell_Type_heuristic', size = 3, angle = 90)


DoHeatmap(sobj1, features = filtered_chromatin_markers_all$gene, group.by = 'New_cellType', size = 3, angle = 90)
DoHeatmap(sobj1, features = rownames(chromatin_markers_astrocytes), group.by = 'New_cellType', size = 3, angle = 90)


DoHeatmap(sobj1, features = as.vector(genes_chromatin$gene), group.by = 'Cell_Type_heuristic', size = 3, angle = 90)



# Analyzing scale data 
scale_data = GetAssayData(sobj1, layer = 'scale.data')
intersect(rownames(scale_data), epigenetic_modifiers)

counts_data = GetAssayData(sobj1, layer = 'counts')
counts_data_epigenetic_modifiers = intersect(rownames(counts_data), epigenetic_modifiers)
dim(counts_data[counts_data_epigenetic_modifiers,])
sanitized_counts = as.matrix(counts_data[counts_data_epigenetic_modifiers,])
#sanitized_counts = apply(, function(x) return(ifelse(is.numeric(x), x, 0)))
rowSums(sanitized_counts)
count(sanitized_counts)
non_zero_counts_per_row <- apply(sanitized_counts, 1, function(row) sum(row != 0))
non_zero_counts_per_row

DoHeatmap(sobj1, features = epigenetic_modifiers, group.by = 'day') + NoLegend()

sobj1_migrating <- subset(sobj1, subset = New_cellType == 'Migrating neurons')

DoHeatmap(sobj1_migrating, features = epigenetic_modifiers, group.by = 'Day', size = 3, angle = 90)

sobj1_e11 <- subset(sobj1, subset = orig_ident == 'E11')

DoHeatmap(sobj1_e11, features = epigenetic_modifiers, group.by = 'New_cellType')

