library(Matrix)

# Read paper1 dataset (matrix, gene names, barcodes), returns sorted gene matrix
read_paper1_data = function() {
  # read gene data matrix
  gene_sorted_matrix = readMM('data/paper1/expression/601ae2f4771a5b0d72588bfb/gene_sorted-matrix.mtx.gz')
  
  # read gene names
  genes = read.table('data/paper1/expression/601ae2f4771a5b0d72588bfb/genes.tsv', sep="\t", header=FALSE)
  gene_names = genes[,1]
  head(gene_names)
  
  # read barcodes
  barcodes = read.table('data/paper1/expression/601ae2f4771a5b0d72588bfb/barcodes.tsv', sep="\t", header=FALSE)
  barcode_names = barcodes[,1]
  
  # Set the row names of gene_sorted_matrix to gene identifiers
  rownames(gene_sorted_matrix) = gene_names
  
  # Set the column names of gene_sorted_matrix to barcode names
  colnames(gene_sorted_matrix) = barcode_names
  
  return(gene_sorted_matrix)
}

# read formatted gene metadata
read_paper1_metadata = function() {
  metadata = read.table('data/paper1/metadata/metaData_scDevSC.txt', sep="\t", header=TRUE, quote='')
  metadata = metadata[-1,]
  colnames(metadata)[colnames(metadata) == "seurat_clusters"] <- "seurat_clusters_orig"
  
  return (metadata)
}
