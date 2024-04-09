celltype_marker_lists2 <- list(
  ProjectionNeurons = c('Neurod2', 'Neurod6', 'Tubb3', 'Tbr1', 'Cux2', 'Cux1', 'Satb2', 'Bhlhe22', 'Hspb3', 'Lmo4', 'Rorb',
                        'Foxp1', 'Fezf2', 'Pcp4', 'Ldb2', 'Etv1', 'Bcl11b', 'Tle4', 'Syt6', 'Foxp2', 'Zfpm2', 'Crym'), 
  # unterteile?
  InhibitoryNeurons = c('Gad1', 'Gad2', 'Slc32a1', 'Erbb4', 'Sst', 'Calb2', 'Vip'),
  OligodendrocytePrecursorCells = c('Fabp7', 'Ccnd1', 'Cspg4', 'Mdk', 'Ednrb', 'Pdgfra'),
  CommittedOligodendrocytePrecursors = c('Brca1', 'Pak4', 'Mycl', 'Pdcd4', 'Fyn', 'Bmp4', 'Epcam'),
  NewlyFormedOligodendrocytes = c('Sema4d', 'Mob3b', 'Ddc', 'Cnksr3', 'Prom1', 'Fam107b', 'Tmem2', 'Rras2', 'Plekha1',
                                  'Kndc1', 'Slc9a3r2', 'Tmem141', 'Man1a', 'Prr5l', 'Aspa', 'Anln', 'Ndrg1'),
  Astrocytes = c('Aldh1l1', 'Slc1a3', 'Apoe', 'Gfap', 'Aqp4'),
  Microglia = c('Tmem119', 'Aif1'),
  Pericytes = c('Pdgfrb')
)








##################################################################################################
# Cell-Type Identification:

library(readr)
library(dplyr) 
library(ggplot2)
library(Seurat)
library(patchwork)
library(Matrix)
library(tidyr)

source('R/paper2/read_paper2_data.R')
source('R/shared/index.R')


sobj2_full <- setup_seurat(
  project = 'paper2',
  counts = read_paper2_data(),
  metadata = read_paper2_metadata()
)

sobj2 <- add_celltype_metadata(sobj2_full, day_column="timePoint", celltype_column="CellType", celltype_marker_lists=celltype_marker_lists2)
# Filter mito %, genes # & nUMI
sobj2 <- subset(sobj2, subset = percent.mito < 0.1 & nGene < 6000 & nGene > 1000 & nUMI < 15000)

# Clustering (Normalizing, Scaling, PCA, Neighbours)
sobj2 <- cluster_seurat(sobj2, dims=1:35, resolution=0.3, nn.eps=0.5)


# genes list by Rodrigo
genes_chromatin <- read.csv("input/genes_chromatin_mm39.csv")

# Delete 2nd column & change header
print(names(genes_chromatin))
genes_chromatin <- genes_chromatin[, !names(genes_chromatin) %in% c('X')]
new_column_names <- c('gene', 'ID', 'whole_gene_name', 'function')
names(genes_chromatin) <- new_column_names
View(genes_chromatin)


#################################################################################################
# Barplot of Cell-Type Composition per Time Point:
head(sobj2)



timePoint <- sobj2@meta.data$timePoint
CellType <- sobj2@meta.data$CellType
topics <- sobj2@meta.data$topics

# Create a data frame with relevant information
df <- data.frame(Timepoint = timePoint, Cell_Type = CellType, Topic = topics)

# Optionally, you might want to summarize your data to get counts of cells per combination of timepoint, cell type, and topic
# For example, using dplyr:
library(dplyr)
df_summary <- df %>%
  group_by(Timepoint, Cell_Type, Topic) %>%
  summarise(Count = n())

# Create the bar plot
ggplot(df_summary, aes(x = Timepoint, y = Count, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Topic) +
  labs(title = "Bar Plot of Cell Types Across Timepoints",
       x = "Timepoint",
       y = "Count",
       fill = "Cell Type") +
  theme_minimal()






CT_TP_plot <- DimPlot(sobj2, reduction="umap", group.by = "CellType", split.by="timePoint", label=T, repel=TRUE) + NoLegend() 

umap_coords <- CT_TP_plot$sobj2

ggplot(umap_coords, aes(x='CellType', y='timePoint', fill='CellType')) +
  geom_bar(stat='identity') +
  labs(title='Cell-Types at Different Time Points',
       x='Cell Types',
       y= 'time Points') +
  theme_minimal()



ggplot(umap_coords, aes(x = "timePoint", y = "CellType")) +
  geom_bin2d() +
  labs(title = "Cell-Types at Different Time Points",
       x = "timePoint",
       y = "CellType") +
  theme_minimal()

##################################################################################################
# Feature Map:
# Cluster Identification
# DimPlot(sobj2, reduction="umap", group.by = "CellType")
# DimPlot(sobj2, reduction="umap", group.by = "timePoint")
# DimPlot(sobj2, reduction="umap", group.by = "combId")
# DimPlot(sobj2, reduction="umap", group.by = "seurat_clusters", label=T)
CT_TP_plot <- DimPlot(sobj2, reduction="umap", group.by = "CellType", split.by="timePoint", label=T, repel=TRUE) + NoLegend() 
CT_TP_plot
ggsave(figure_filename(sobj2@project.name, 'CellTypeTimePoint_plot'), plot = CT_TP_plot, width = 10, height = 8, dpi = 300)

#MG_TP_plot <- DimPlot(sobj2, reduction="umap", group.by = "matches", split.by="timePoint", label=T, repel=TRUE) + NoLegend() 
#MG_TP_plot
#ggsave(figure_filename(sobj2@project.name, 'MargerGenesTimePoint_plot'), plot = CT_TP_plot, width = 10, height = 8, dpi = 300)
# Look at cluster IDs of the first 5 cells
# head(Idents(sobj2), 5)

# Loop over the list, create a FeaturePlot for each set of markers & save as PNG
for (celltype_markers in names(matches)) {
  plot_feature_set(sobj2, matches[[celltype_markers]], celltype_markers)
}

# add markers_score for all celltypes
for (celltype in names(matches)) {
  sobj2 <- add_celltype_markers_score(sobj2, matches[[celltype]], celltype)
}

# plot all celltype markers score in one PNG file
celltype_markers_score_features = sapply(names(celltype_marker_lists2), celltype_markers_score_col_name)
plot_feature_set(sobj2, celltype_markers_score_features, 'celltype_markers_scores', ncol = 6)

# find all markers of cluster 1
# cluster1.markers <- FindMarkers(sobj2, ident.1 = 2)
# head(cluster1.markers, n = 5)

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


##################################################################################################
# Heat Maps:

for (cluster in names(genes_chromatin)) {
  markers <- cluster_marker_lists[[cluster]]
  if (length(markers) > 0) {
    marker_set_name <- paste("Cluster", cluster, "TopMarkers")
    # Call the plot_marker_set function for the current set of markers
    plot_marker_set(sobj2, markers, marker_set_name)
  }
}

# gather celltype assignment lists
celltype_assignments = assign_cell_types(cluster_marker_lists2, celltype_marker_lists2)
celltype_assignments
print(celltype_assignments)

all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
top10$gene
DoHeatmap(sobj2, features = top10$gene) + NoLegend()
celltype_marker_lists['cpnLayer23']
DoHeatmap(sobj2, features = celltype_marker_lists['cpnLayer56']) + NoLegend()

# Heatmap for epigenetic modifiers
DoHeatmap(sobj2, features = genes) + NoLegend()
genes <- as.vector(genes_chromatin$gene)

DoHeatmap(object=sobj2, features = as.vector(genes_chromatin$gene)) + NoLegend()

          

##################################################################################################
# Differential Expression: Of Astrocytes against most abundant 

# same TPs differ. CTs

# Plot heatmap, sufficient # of genes, visually clear

# 1 c. type over differ. days








#FeaturePlot(sobj2, features = CPN_markers , ncol = 2, order = TRUE)
#VlnPlot(sobj2, features = CPN_markers, ncol = 2)

#Function (marker_sets), does featureplots & violinplots & saves as image
#plot_feature_set(sobj2, CPN_markers, "CPN")

# umaps of gene expression
#FeaturePlot(sobj2, features=c(celltype_markers_score_col_name("CPN")))
#sobj2 <- add_celltype_markers_score(sobj2, CPN_markers, "CPN")

# Heatmaps
DoHeatmap(object=sobj2, features = as.vector(genes_chromatin$genes)) + NoLegend()

# Marker gene expression by cell
