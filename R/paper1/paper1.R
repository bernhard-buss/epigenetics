#install.packages('Seurat')

library('Seurat')
library('ggplot2')
library(dplyr)

source('R/shared/index.R')
source('R/paper1/read_paper1_data.R')
source('R/paper1/celltypes_paper1.R')

# Create Seurat object with min genes = 500
sobj1 <- setup_seurat(
  project='paper1',
  counts=read_paper1_data(),
  metadata=read_paper1_metadata()
)

# Filter mito < 7.5 % && genes > 500 (& New_cellType in celltypes)
#sobj <- subset(sobj, subset = percent_mito < 0.075 & nFeature_RNA > 500) # includes DL (deleterious) and low quality
#sobj <- subset(sobj, subset = percent_mito < 0.075 & nFeature_RNA > 500 & New_cellType %in% celltypes)
sobj1 <- subset(sobj1, subset = percent_mito < 0.075 & nFeature_RNA > 500 & New_cellType %in% celltypes & Phase == 'G1')

sobj1 <- cluster_seurat(sobj1, nfeatures=3000, resolution=1)

DimPlot(sobj1, reduction = "umap", group.by = 'seurat_clusters_orig')
DimPlot(sobj1, reduction = "umap", group.by = 'seurat_clusters', label = TRUE)

origIdentPlot <- DimPlot(sobj1, reduction = "umap", group.by = 'orig_ident')
origIdentPlot

cellTypePlot <- DimPlot(sobj1, reduction = "umap", group.by = 'New_cellType', label = TRUE, repel = TRUE) + NoLegend()
cellTypePlot
ggsave(figure_filename(sobj1@project.name, 'New_cellType_plot'), plot = cellTypePlot, width = 10, height = 8, dpi = 300)


# Look at cluster IDs of the first 5 cells
head(Idents(sobj1), 5)

# Loop over the list, create a FeaturePlot for each set of markers, and save as PNG
for (celltype in names(celltype_marker_lists)) {
  plot_feature_set(sobj1, celltype_marker_lists[[celltype]], celltype)
}

# add markers_score for all celltypes
for (celltype in names(celltype_marker_lists)) {
  sobj1 <- add_celltype_markers_score(sobj1, celltype_marker_lists[[celltype]], celltype)
}
# plot all celltype markers score in one PNG file
celltype_markers_score_features = sapply(names(celltype_marker_lists), celltype_markers_score_col_name)
plot_feature_set(sobj1, celltype_markers_score_features, 'celltype_markers_scores', ncol = 6)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(sobj1, ident.1 = 2)
head(cluster1.markers, n = 5)

# Find all markers
all_markers <- FindAllMarkers(
  object = sobj1,
  only.pos = TRUE, # Consider only positive markers
  min.pct = 0.50, #0.25, # Gene must be detected in at least 25% of cells within a cluster
  #min.diff.pct = 0.5,
  logfc.threshold = 0.7, #0.25, # Minimum log-fold change
  test.use = 'wilcox', # Use Wilcoxon Rank Sum test
  p.adjust.method = 'bonferroni' # Bonferroni correction for multiple testing
)

head(all_markers)

# Filter markers based on adjusted P-value < 0.05
significant_markers <- all_markers[all_markers$p_val_adj < 0.05, ]
head(significant_markers)

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
top_markers_per_cluster

cluster_marker_lists <- split(top_markers_per_cluster$gene, top_markers_per_cluster$cluster)

for (cluster in names(cluster_marker_lists)) {
  markers <- cluster_marker_lists[[cluster]]
  if (length(markers) > 0) {
    marker_set_name <- paste("Cluster", cluster, "TopMarkers")
    # Call the plot_marker_set function for the current set of markers
    plot_marker_set(sobj, markers, marker_set_name)
  }
}

# gather celltype assignment lists
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

# Heatmap for epigenetic modifiers
DoHeatmap(sobj1, features = epigenetic_modifiers) + NoLegend()
