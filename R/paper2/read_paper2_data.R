# Read paper2 dataset (matrix), returns gene matrix
read_paper2_data = function() {
  # read gene data matrix
  gene_matrix <- readRDS("data/paper2/GSE204759_mouse_scRNA_raw_counts.rds")
  
  return(gene_matrix)
}

# read formatted gene metadata
read_paper2_metadata = function() {
  metadata <- read.table("data/paper2/GSE204759_mouse_scRNA_metadata.txt", sep="\t", header=T, quote='')
  rownames(metadata)<-metadata$cellId
  
  return (metadata)
}
