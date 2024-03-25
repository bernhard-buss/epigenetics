# celltype names using in paper1
celltypes = c(
  'Apical progenitors',
  'Astrocytes',
  'Cajal Retzius cells',
  'CThPN',
#  'Cycling glial cells',
  'Endothelial cells',
  'Immature neurons',
  'Intermediate progenitors',
  'Interneurons',
  'Layer 4',
  'Layer 6b',
  'Microglia',
  'Migrating neurons',
  'NP',
  'Oligodendrocytes',
  'Pericytes',
  'Red blood cells',
  'SCPN',
  'VLMC'
)

# positive marker genes used for classification of cell types
celltype_marker_lists <- list(
  # Neuronal cell markers
  apicalProgenitors = c('Sox2', 'Pax6', 'Hes5', 'Tubb3', 'Neurog2', 'Btg2'),
  intermediateProgenitors = c('Neurog2', 'Tubb3', 'Eomes', 'Neurod6', 'Btg2', 'Pax6', 'Sox2', 'Hes5', 'Neurod1', 'Nrp1', 'Pou3f2'),
  migratingNeurons = c('Neurod2', 'Neurod6', 'Tubb3', 'Nrp1', 'Neurod1', 'Pou3f2'),
  immatureNeurons = c('Neurod2', 'Neurod6', 'Tubb3'),
  cajalretziusNeurons = c('Reln', 'Lhx5', 'Tubb3'),
  interNeurons = c('Tubb3', 'Dlx2', 'Gad1', 'Sox2', 'Gad2'), #Gad2 from paper2
  ependymocytes = c('Hes5', 'Foxj1', 'Sox2', 'Btg2', 'Pax6', 'Tubb3'),
  # cortigo fugal neuron markers
  thalamusProjecting = c('Sox5', 'Bcl11b', 'Tle4', 'Fezf2'),
  subcerebralProjecting = c('Sox5', 'Fezf2', 'Bcl11b', 'Ldb2'),
  nearProjecting = c('Fezf2', 'Bcl11b', 'Fle4', 'Sox5', 'Ldb2'),
  layer6b = c('Bcl11b', 'Tle4', 'Pcp4', 'Nr4a2', 'Sox5'),
  # intra telencephalic neurons
  stellateLayer4 = c('Satb2', 'Ptn', 'Pantr1', 'Lpl', 'Cux1', 'Rorb'),
  cpnLayer23 = c('Pantr1', 'Ptn', 'Satb2', 'Cux1', 'Cux2', 'Lmo4'),
  cpnLayer56 = c('Pantr1', 'Ptn', 'Satb2', 'Lmo4', 'Ldb2', 'Rora'),
  # Glial cell markers
  astrocytes = c('Apoe', 'Slc1a3', 'Cst3', 'Aldh1l1', 'Olig1', 'Aqp4'),
  oligodendrocytes = c('Pdgfra', 'Olig1', 'Cst3', 'Olig2'),
  microglia = c('P2ry12', 'Apoe', 'Cst3', 'Aif1', 'C1qb'),
  cycling = c('Top2a', 'Cst3', 'Slc1a3', 'Pcna', 'Apoe'),
  # vascular cell markers
  endothelial = c('Cldn5', 'Igfbp7', 'Car2', 'Rgs5'),
  vlmc = c('Col3a1', 'Lum', 'Lgals1'),
  pericytes = c('Lgals1', 'Pdgfrb', 'Igfbp7', 'Rgs5', 'Cspg4', 'Col3a1'),
  redBlood = c('Hemgn', 'Car2', 'Lgals1')
)

# top epigenetic modifier genes
dna_methylation_genes <- c("Dnmt1", "Dnmt3a", "Tet1", "Tet2", "Tet3", "Mecp2")
histone_modification_genes <- c("Ezh2", "Hdac1", "Hdac2", "Kmt2a", "Suv39H1", "Ash1l", "Setd1a")
chromatin_loop_related_genes <- c("Ctcf", "Rad21", "Smc3", "Smc1a", "Smc1b", "Stag1", "Stag2", 
                                  "Wapl", "Nipbl", "Mau2", "Znf143", "Ctcfl", "Tcf12")
dna_binding_protein_genes <- c("Yy1", "Sox2", "Pou3f2", "Atf2")

epigenetic_modifiers <- c(dna_methylation_genes, histone_modification_genes, chromatin_loop_related_genes, dna_binding_protein_genes)


