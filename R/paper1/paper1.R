#install.packages('Seurat')

library('Seurat')
library('ggplot2')
library(dplyr)

source('R/shared/index.R')
source('R/paper1/read_paper1_data.R')

project = 'paper1'

# read experiment data and metadata
gene_data = read_paper1_data()
metadata = read_paper1_metadata()

# Create Seurat object with min genes = 500
sobj <- CreateSeuratObject(counts = gene_data, project = project)
sobj <- AddMetaData(sobj, metadata=metadata)

# Filter mito < 7.5 % && genes > 500 (& New_cellType in celltypes)
#sobj <- subset(sobj, subset = percent_mito < 0.075 & nFeature_RNA > 500) # includes DL (deleterious) and low quality
#sobj <- subset(sobj, subset = percent_mito < 0.075 & nFeature_RNA > 500 & New_cellType %in% celltypes)
sobj <- subset(sobj, subset = percent_mito < 0.075 & nFeature_RNA > 500 & New_cellType %in% celltypes & Phase == 'G1')

# Normalize
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

# Find Variable Features
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sobj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Scale data
#sobj <- ScaleData(sobj, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent_mito', 'CC_Difference'), do.centre=TRUE, do.scale=TRUE)
#saveRDS(sobj, file="sobj_scaled_full_paper.rds")
#sobj <- ScaleData(sobj, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent_mito'), do.scale=TRUE)
sobj <- ScaleData(sobj)


# PCA
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
DimPlot(sobj, reduction = "pca") + NoLegend()
ElbowPlot(sobj, ndims=50)

sobj <- FindNeighbors(sobj, dims = 1:50)
sobj <- FindClusters(sobj, algorithm = 1, method = 'matrix', resolution=1)

# UMAP
sobj <- RunUMAP(sobj, dims = 1:50)
DimPlot(sobj, reduction = "umap", group.by = 'seurat_clusters_orig')
DimPlot(sobj, reduction = "umap", group.by = 'seurat_clusters', label = TRUE)

origIdentPlot <- DimPlot(sobj, reduction = "umap", group.by = 'orig_ident')
origIdentPlot

cellTypePlot <- DimPlot(sobj, reduction = "umap", group.by = 'New_cellType', label = TRUE, repel = TRUE) + NoLegend()
cellTypePlot
ggsave(figure_filename(sobj@project.name, 'New_cellType_plot'), plot = cellTypePlot, width = 10, height = 8, dpi = 300)


# Look at cluster IDs of the first 5 cells
head(Idents(sobj), 5)

# Loop over the list, create a FeaturePlot for each set of markers, and save as PNG
for (celltype in names(celltype_marker_lists)) {
  plot_feature_set(sobj, celltype_marker_lists[[celltype]], celltype)
}

# add markers_score for all celltypes
for (celltype in names(celltype_marker_lists)) {
  sobj <- add_celltype_markers_score(sobj, celltype_marker_lists[[celltype]], celltype)
}
# plot all celltype markers score in one PNG file
celltype_markers_score_features = sapply(names(celltype_marker_lists), celltype_markers_score_col_name)
plot_feature_set(sobj, celltype_markers_score_features, 'celltype_markers_scores', ncol = 6)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(sobj, ident.1 = 2)
head(cluster1.markers, n = 5)

# Find all markers
all_markers <- FindAllMarkers(
  object = sobj,
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
DoHeatmap(sobj, features = top10$gene) + NoLegend()
celltype_marker_lists['cpnLayer23']
DoHeatmap(sobj, features = celltype_marker_lists['cpnLayer56']) + NoLegend()

# Heatmap for epigenetic modifiers
DoHeatmap(sobj, features = epigenetic_modifiers) + NoLegend()
