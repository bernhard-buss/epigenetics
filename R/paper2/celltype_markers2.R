library(readr)

# positive marker genes used for classification of cell types
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

# genes list by Rodrigo
genes_chromatin <- read.csv("data/paper2/genes_chromatin_mm39.csv")

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

# genes that match between the 2 lists:
for(celltype_markers in celltype_marker_lists2) {
  matches <- grepl(celltype_markers, genes_chromatin$gene)
              num_matches <- sum(matches)
              cat('Cell Type Markers', celltype_markers, '\n')
              cat('Matches', num_matches, '\n')
}

?grepl

