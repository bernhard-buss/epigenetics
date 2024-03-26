library(readr)
library(dplyr) 
library(ggplot2)
library(Seurat)
library(patchwork)
library(Matrix)

source('R/paper2/read_paper2_data.R')
source('R/shared/index.R')

sobj2 <- setup_seurat(
  project='paper2',
  counts=read_paper2_data(),
  metadata=read_paper2_metadata()
)

sobj2 <- subset(sobj2, subset=percent.mito<0.1 & nGene<6000 & nGene>1000 & nUMI<15000)

sobj2 <- cluster_seurat(sobj2, dims=1:35, resolution=0.3, nn.eps=0.5)

# Cluster Identification
DimPlot(sobj2, reduction="umap", group.by = "CellType", split.by="timePoint", label=T, repel=TRUE) + NoLegend() 
DimPlot(sobj2, reduction="umap", group.by = "timePoint")
DimPlot(sobj2, reduction="umap", group.by = "combId")
DimPlot(sobj2, reduction="umap", group.by = "seurat_clusters", label=T)
# Look at cluster IDs of the first 5 cells ?
# head(Idents(sobj2), 5) ?

# Average expression for each cluster
# -> Compute joint cell type scores, for each cluster & potential identity 

# Find all markers of each cluster
#markers <- FindAllMarkers(sobj2, only.pos=T, min.pct=0.5, logfc.threshold=0.7, test.use='wilcox', p.adjust.method='bonferroni', min.diff.pct = 0.5)
#saveRDS(markers,"Documents/Neuroepigenetics/markers.rds")
markers <- readRDS("tmp/markers.rds")
#genenames <- rownames(markers)
#genenames

CPN_markers <- c("Cux2", "Lcorl", "Rora")
CThPN_markers <- c("Bcl11b", "Sox8", "Tle4", "Fezf2")
PN_markers <- c("Neurod2")
IN_markers <- c("Gad2", 'Tubb3', 'Dlx2', 'Gad1')
Oligo_markers <- c("Pdgfra", 'Olig1', 'Cst3', 'Olig2')
Astro_markers <- c("Aqp4", 'Apoe', 'Slc1a3', 'Cst3', 'Aldh1l1', 'Olig1')

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
