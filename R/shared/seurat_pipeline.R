library(Seurat)

# Create SeuratObject and setup with metadata.
setup_seurat <- function(project, counts, metadata) {
  #Setup the Seurat Object & filtering: (min genes = 500, mito<5%)
  sobj <- CreateSeuratObject(counts=counts, project=project)
  sobj <- AddMetaData(sobj, metadata=metadata)
  
  return(sobj)
}

# Normalize, FindVariableFeatures. Scale and PCA.
cluster_seurat <- function(sobj, nfeatures = 1500, dims=1:50, min.dist=0.3, resolution=0.5, nn.eps=0.5) {
  # Normalize
  sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find Variable Features
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = nfeatures)
  
  # Scale and PCA
  sobj <- ScaleData(sobj)
  #sobj <- ScaleData(sobj, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent_mito', 'CC_Difference'), do.centre=TRUE, do.scale=TRUE)
  #sobj <- ScaleData(sobj, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent_mito'), do.scale=TRUE)

  sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
  
  # Identify the 10 most highly variable genes
  #top10 <- head(VariableFeatures(seurato), 10)
  
  # plot variable features with and without labels
  #plot1 <- VariableFeaturePlot(seurato)
  #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  #plot1 + plot2
  
  # Visualization
  #DimPlot(seurato, reduction = "pca") + NoLegend()
  #ElbowPlot(seurato, ndims=50)
  
  # UMAP: To assign identities, examination of expression of all marker genes individually
  sobj <- RunUMAP(sobj, dims = dims, min.dist = min.dist)
  
  sobj <- FindNeighbors(sobj, ndims=dims)
  sobj <- FindClusters(sobj, algorithm=1, method='matrix', resolution=resolution, nn.eps=nn.eps)
  
  return(sobj)
}
