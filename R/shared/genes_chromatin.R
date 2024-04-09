# genes list by Rodrigo
genes_chromatin <- read.csv("input/genes_chromatin_mm39.csv")

# Delete 2nd column & change header
print(names(genes_chromatin))
genes_chromatin <- genes_chromatin[, !names(genes_chromatin) %in% c('X')]
new_column_names <- c('gene', 'ID', 'whole_gene_name', 'function')
names(genes_chromatin) <- new_column_names
View(genes_chromatin)
# check
# print(names(genes_chromatin))
# View(genes_chromatin)
# head(genes_chromatin)