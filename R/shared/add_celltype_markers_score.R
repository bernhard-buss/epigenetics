library(Seurat)

# return the markers score column name for the given celltype
celltype_markers_score_col_name <- function(celltype) {
  return(paste0(celltype, '_markers_score'))
}

# return the normalized markers score column name for the given celltype
celltype_markers_score_norm_col_name <- function(celltype) {
  return(paste0(celltype, '_markers_score_norm'))
}

# calculate and add celltype markers score columns to the seurat object
add_celltype_markers_score <- function(sobj, markers, celltype, fc_threshold = 1) {
  markers_passing_threshold = FetchData(sobj, vars = markers) > fc_threshold
  # Count the number of genes that exceed the threshold for each cell
  celltype_markers_score <- rowSums(markers_passing_threshold)
  sobj <- AddMetaData(sobj, metadata = celltype_markers_score, col.name = celltype_markers_score_col_name(celltype))
  # normalize the score (divide by number of markers) and square it
  celltype_markers_score_normalized <- (celltype_markers_score / length(markers))^2
  sobj <- AddMetaData(sobj, metadata = celltype_markers_score_normalized, col.name = celltype_markers_score_norm_col_name(celltype))
  
  return(sobj)
}
