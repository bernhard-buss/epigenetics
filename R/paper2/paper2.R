library(readr)
library(dplyr) 
library(ggplot2)
library(Seurat)
library(patchwork)
library(Matrix)

source('R/paper2/read_paper2_data.R')
source('R/shared/index.R')

sobj2 <- setup_seurat(
  project = 'paper2',
  counts = read_paper2_data(),
  metadata = read_paper2_metadata()
)


# Filter mito %, genes # & nUMI
sobj2 <- subset(sobj2, subset = percent.mito < 0.1 & nGene < 6000 & nGene > 1000 & nUMI < 15000)

# Clustering (Normalizing, Scaling, PCA, Neighbours)
sobj2 <- cluster_seurat(sobj2, dims=1:35, resolution=0.3, nn.eps=0.5)

# Cluster Identification
# DimPlot(sobj2, reduction="umap", group.by = "CellType")
# DimPlot(sobj2, reduction="umap", group.by = "timePoint")
# DimPlot(sobj2, reduction="umap", group.by = "combId")
# DimPlot(sobj2, reduction="umap", group.by = "seurat_clusters", label=T)
CT_TP_plot <- DimPlot(sobj2, reduction="umap", group.by = "CellType", split.by="timePoint", label=T, repel=TRUE) + NoLegend() 
CT_TP_plot
ggsave(figure_filename(sobj2@project.name, 'CellTypeTimePoint_plot'), plot = CT_TP_plot, width = 10, height = 8, dpi = 300)

# Look at cluster IDs of the first 5 cells
# head(Idents(sobj2), 5)

# Loop over the list, create a FeaturePlot for each set of markers & save as PNG
for (celltype in names(celltype_marker_lists)) {
  plot_feature_set(sobj2, celltype_marker_lists[[celltype]], celltype)
}

# add markers_score for all celltypes
for (celltype in names(celltype_marker_lists)) {
  sobj2 <- add_celltype_markers_score(sobj2, celltype_marker_lists[[celltype]], celltype)
}

# plot all celltype markers score in one PNG file
celltype_markers_score_features = sapply(names(celltype_marker_lists), celltype_markers_score_col_name)
plot_feature_set(sobj2, celltype_markers_score_features, 'celltype_markers_scores', ncol = 6)

# find all markers of cluster 1
# cluster1.markers <- FindMarkers(sobj2, ident.1 = 2)
# head(cluster1.markers, n = 5)

# Find all markers of each cluster
# Old list of markers, incorporated into total list
# CPN_markers <- c("Cux2", "Lcorl", "Rora")
# CThPN_markers <- c("Bcl11b", "Sox8", "Tle4", "Fezf2")
# PN_markers <- c("Neurod2")
# IN_markers <- c("Gad2", 'Tubb3', 'Dlx2', 'Gad1')
# Oligo_markers <- c("Pdgfra", 'Olig1', 'Cst3', 'Olig2')
# Astro_markers <- c("Aqp4", 'Apoe', 'Slc1a3', 'Cst3', 'Aldh1l1', 'Olig1')
all_markers <- FindAllMarkers(
  object = sobj2,
  only.pos = TRUE, # Consider only positive markers
  min.pct = 0.50, #0.25, # Gene must be detected in at least 25% of cells within a cluster
  #min.diff.pct = 0.5,
  logfc.threshold = 0.7, #0.25, # Minimum log-fold change
  test.use = 'wilcox', # Use Wilcoxon Rank Sum test
  p.adjust.method = 'bonferroni' # Bonferroni correction for multiple testing
)
# head(all_markers)
# print(rownames(all_markers))

# Filter markers based on adjusted P-value < 0.05
significant_markers <- all_markers[all_markers$p_val_adj < 0.05, ]
# head(significant_markers)

# Sort markers within each cluster by adj. p-value & log-fold change
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
    plot_marker_set(sobj2, markers, marker_set_name)
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








FeaturePlot(sobj2, features = CPN_markers , ncol = 2, order = TRUE)
VlnPlot(sobj2, features = CPN_markers, ncol = 2)

#Function (marker_sets), does featureplots & violinplots & saves as image
plot_feature_set(sobj2, CPN_markers, "CPN")
plot_feature_set(sobj2, CThPN_markers, "CThPN")
plot_feature_set(sobj2, PN_markers, "PN")
plot_feature_set(sobj2, IN_markers, "IN")
plot_feature_set(sobj2, Oligo_markers, "Oligo")
plot_feature_set(sobj2, Astros_markers, "Astro")

# umaps of gene expression
sobj2 <- add_celltype_markers_score(sobj2, CPN_markers, "CPN")
FeaturePlot(sobj2, features=c(celltype_markers_score_col_name("CPN")))

sobj2 <- add_celltype_markers_score(sobj2, CThPN_markers, "CThPN")
FeaturePlot(sobj2, features=c(celltype_markers_score_col_name("CThPN")))

sobj2 <- add_celltype_markers_score(sobj2, PN_markers, "PN")
FeaturePlot(sobj2, features=c(celltype_markers_score_col_name("PN")))

sobj2 <- add_celltype_markers_score(sobj2, IN_markers, "IN")
FeaturePlot(sobj2, features=c(celltype_markers_score_col_name("IN")))

sobj2 <- add_celltype_markers_score(sobj2, PN_markers, "Oligo")
FeaturePlot(sobj2, features=c(celltype_markers_score_col_name("Oligo")))

sobj2 <- add_celltype_markers_score(sobj2, PN_markers, "Astro")
FeaturePlot(sobj2, features=c(celltype_markers_score_col_name("Astro")))

# Heatmaps
DoHeatmap(object=sobj2, features = CPN_markers) + NoLegend()
DoHeatmap(object=sobj2, features = CThPN_markers) + NoLegend()
DoHeatmap(object=sobj2, features = PN_markers) + NoLegend()
DoHeatmap(object=sobj2, features = IN_markers) + NoLegend()
DoHeatmap(object=sobj2, features = Oligo_markers) + NoLegend()
DoHeatmap(object=sobj2, features = Astro_markers) + NoLegend()

# Marker gene expression by cell
