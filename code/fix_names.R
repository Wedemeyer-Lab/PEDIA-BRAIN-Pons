# load packages
library(Signac)
library(Seurat)

# read in data
atlas <- readRDS("output/objects/atlas.RDS")

# Select one celltype column
atlas@meta.data$Cell_Type <- atlas@meta.data$Cell_Type_General
atlas@meta.data$Cell_Type_General <- NULL
atlas@meta.data$Cell_Type_Res0.6 <- NULL

# replace "Oligo" with "OL"
atlas@meta.data$Cell_Type <- gsub("Oligo","OL",atlas@meta.data$Cell_Type)
atlas@meta.data$Subtype <- gsub("Oligo","OL",atlas@meta.data$Subtype)

# replace "Neuron" with "Motor Neuron"
atlas@meta.data$Cell_Type[atlas@meta.data$Cell_Type == "Neuron"] <- "Motor Neuron"

##we customed cell cycle genes (for justification please see the main article)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- subset(g2m.genes, !(g2m.genes %in% c('ANLN','PSRC1','CKAP5', 'TUBB4B')))

##to perform cell cycle/phase calculation
DefaultAssay(atlas) <- "RNA"
atlas <- NormalizeData(atlas)
atlas <- CellCycleScoring(atlas, s.features = s.genes, g2m.features = g2m.genes)
atlas$CC.Difference <- atlas$S.Score - atlas$G2M.Score

# save atlas
saveRDS(atlas,"output/objects/atlas_final.RDS")

