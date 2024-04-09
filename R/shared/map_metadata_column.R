# add a metadata column using a lambda mapping function
map_metadata_column <- function(sobj, source, target, lambda) {
  # Apply the lambda function to the source column
  result <- sapply(sobj@meta.data[[source]], lambda)
  
  # Ensure the result has row names matching the Seurat object's cell names
  names(result) <- rownames(sobj@meta.data)
  
  # Use AddMetaData to add the transformed data as a new column
  sobj <- AddMetaData(sobj, metadata = result, col.name = target)
  
  return(sobj)
}

# add a metadata column using a lambda mapping function from multiple source columns
map_metadata_columns <- function(sobj, sources, target, lambda) {
  # Extract the source columns and apply the lambda function to each row
  source_data <- sobj@meta.data[, sources, drop = FALSE] # Ensure it's a dataframe even if single column
  result <- apply(source_data, 1, lambda)
  
  # Ensure the result has row names matching the Seurat object's cell names
  names(result) <- rownames(sobj@meta.data)
  
  # Use AddMetaData to add the transformed data as a new column
  sobj <- AddMetaData(sobj, metadata = result, col.name = target)
  
  return(sobj)
}

# lambda function to extract a prefix before '_'
extract_prefix <- function(x) sub("_.*", "", x)

# lambda function to find the name of the column with the highest score
highest_score_col_name <- function(row) {
  max_col_index <- which.max(row)
  if (max_col_index >= 0) {
    return(names(row)[max_col_index])
  }
  return('unknown')
}
