library(readr)
library(dplyr) 
library(ggplot2)
library(Seurat)
library(patchwork)
library(Matrix)

#Importing dataset:
metad <- read.table("GSE204759_mouse_scRNA_metadata.txt", sep="\t", header=T, quote='')
head(metad)
rownames(metad)<-metad$cellId

matrix <- readRDS("GSE204759_mouse_scRNA_raw_counts.rds")
head(matrix)


#Setup the Seurat Object & filtering: (min genes = 500, mito<5%)
seurato <- CreateSeuratObject(counts=matrix, project="seurat")

seurato <- AddMetaData(seurato, metadata=metad)

seurato <- subset(seurato, subset=percent.mito<0.1 & nGene<6000 & nGene>1000 & nUMI<15000)
head(seurato@meta.data)

# Normalize
seurato <- NormalizeData(seurato, normalization.method = "LogNormalize", scale.factor = 10000)

# Find Variable Features
seurato <- FindVariableFeatures(seurato, selection.method = "vst", nfeatures = 1500)
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(seurato), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(seurato)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

# Scale and PCA
seurato <- ScaleData(seurato)
seurato <- RunPCA(seurato, features = VariableFeatures(object = seurato))

# Visualization
#DimPlot(seurato, reduction = "pca") + NoLegend()
#ElbowPlot(seurato, ndims=50)

# UMAP: To assign identities, examination of expression of all marker genes individually
seurato <- RunUMAP(seurato, dims = 1:35, min.dist = 0.2)
head(seurato)

# Cluster Identification
seurato <- FindNeighbors(seurato, ndims=1:10)
seurato <- FindClusters(seurato, algorithm=1, method='matrix', resolution=0.3, nn.eps=0.5)
DimPlot(seurato, reduction="umap", group.by = "CellType") 
DimPlot(seurato, reduction="umap", group.by = "timePoint") 
DimPlot(seurato, reduction="umap", group.by = "combId")
DimPlot(seurato, reduction="umap", group.by = "seurat_clusters", label=T)
# Look at cluster IDs of the first 5 cells ?
# head(Idents(seurato), 5) ?

# Average expression for each cluster
# -> Compute joint cell type scores, for each cluster & potential identity 

# Find all markers of each cluster
#markers <- FindAllMarkers(seurato, only.pos=T, min.pct=0.5, logfc.threshold=0.7, test.use='wilcox', p.adjust.method='bonferroni', min.diff.pct = 0.5)
#saveRDS(markers,"Documents/Neuroepigenetics/markers.rds")
markers <- readRDS("markers.rds")
#genenames <- rownames(markers)
#genenames

CPN_markers <- c("Cux2", "Lcorl", "Rora")
CThPN_markers <- c("Bcl11b", "Sox8", "Tle4", "Fezf2")
PN_markers <- c("Neurod2")
IN_markers <- c("Gad2", 'Tubb3', 'Dlx2', 'Gad1')
Oligo_markers <- c("Pdgfra", 'Olig1', 'Cst3', 'Olig2')
Astro_markers <- c("Aqp4", 'Apoe', 'Slc1a3', 'Cst3', 'Aldh1l1', 'Olig1')

FeaturePlot(seurato, features = CPN_markers , ncol = 2, order = TRUE)
VlnPlot(seurato, features = CPN_markers, ncol = 2)

#Function (marker_sets), does featureplots & violinplots & saves as image
source("plot_marker_set.R")
plot_marker_set(seurato, CPN_markers, "CPN")
plot_marker_set(seurato, CThPN_markers, "CThPN")
plot_marker_set(seurato, PN_markers, "PN")
plot_marker_set(seurato, IN_markers, "IN")
plot_marker_set(seurato, Oligo_markers, "Oligo")
plot_marker_set(seurato, Astros_markers, "Astro")

# umaps of gene expression
source("add_celltype_markers_score.R")
seurato <- add_celltype_markers_score(seurato, CPN_markers, "CPN")
FeaturePlot(seurato, features=c(celltype_markers_score_col_name("CPN")))

seurato <- add_celltype_markers_score(seurato, CThPN_markers, "CThPN")
FeaturePlot(seurato, features=c(celltype_markers_score_col_name("CThPN")))

seurato <- add_celltype_markers_score(seurato, PN_markers, "PN")
FeaturePlot(seurato, features=c(celltype_markers_score_col_name("PN")))

seurato <- add_celltype_markers_score(seurato, IN_markers, "IN")
FeaturePlot(seurato, features=c(celltype_markers_score_col_name("IN")))

seurato <- add_celltype_markers_score(seurato, PN_markers, "Oligo")
FeaturePlot(seurato, features=c(celltype_markers_score_col_name("Oligo")))

seurato <- add_celltype_markers_score(seurato, PN_markers, "Astro")
FeaturePlot(seurato, features=c(celltype_markers_score_col_name("Astro")))

# Heatmaps
DoHeatmap(object=seurato, features = CPN_markers) + NoLegend()
DoHeatmap(object=seurato, features = CThPN_markers) + NoLegend()
DoHeatmap(object=seurato, features = PN_markers) + NoLegend()
DoHeatmap(object=seurato, features = IN_markers) + NoLegend()
DoHeatmap(object=seurato, features = Oligo_markers) + NoLegend()
DoHeatmap(object=seurato, features = Astro_markers) + NoLegend()




# find all markers of cluster 1
#cluster1.markers <- FindMarkers(seurato, ident.1 = 2)
#head(cluster1.markers, n = 5)

# subset giving cluster IDs -> differential cluster analysis











# supplementary matrix, umaps per time point
# genetic modifiers: Genes list by Rodruigo
# negativemarkers (f.e. from progenitors)
# plot that simply shows as much info at once as possible (time, cell typs)
# define specific goals

