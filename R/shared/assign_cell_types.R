# WIP: datastructure printing is not nice

# Function to assign cell types to clusters based on marker overlap
assign_cell_types <- function(cluster_marker_lists, celltype_marker_lists) {
  # Initialize a list to hold the cell type assignments for each cluster
  cell_type_assignments <- list()
  
  # Loop over each cluster
  for (cluster in names(cluster_marker_lists)) {
    cluster_markers <- cluster_marker_lists[[cluster]]

    # Initialize a vector to hold matching markers for each cell type
    # overlap_markers <- integer(length(celltype_marker_lists))
    overlap_markers = list()

    # Loop over each cell type
    for (cell_type in names(celltype_marker_lists)) {
      # Calculate the overlap between the cluster's markers and the cell type's markers
      matches = intersect(cluster_markers, celltype_marker_lists[[cell_type]])
      if (length(matches) > 0) {
        overlap_markers[[cell_type]] <- matches
      }
    }

    # Determine the cell types with the highest overlap score for this cluster
    if (length(overlap_markers) > 0) {
      str(overlap_markers)
      overlap_lengths <- sapply(overlap_markers, length)
      max_score <- max(overlap_lengths)
      best_matches <- names(which(overlap_lengths == max_score))
      
      # Assign the cell type with the highest score to the cluster
      # If there are ties, all tied cell types are listed
      if (max_score > 0) {
        cell_type_assignments[[cluster]] <- best_matches
      }
    }
  }
  
  return(cell_type_assignments)
}
