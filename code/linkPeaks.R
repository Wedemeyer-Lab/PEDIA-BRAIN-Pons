#Load packages
library(ggplot2)
library(Seurat)
library(Matrix)
library(Signac)
library(tidyverse)
library(dplyr)
library(readr)
library(openxlsx) 
library(clustree)
library(motifmatchr)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(readxl)
library(harmony)
library(ComplexHeatmap)
library(SummarizedExperiment)
set.seed(6)

options(future.globals.maxSize = 200 * 1024^3)

sobj <- readRDS("output/objects/atlas.RDS")

library(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(sobj) <- "ATAC"

sobj <- RegionStats(
  object = sobj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC"
)


sobj <- LinkPeaks(
  object = sobj,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)



saveRDS(sobj, "output/objects/atlas_peaks_linked.RDS")

meta<- sobj@meta.data
cells_over1k <- rownames(meta[meta$Cell_Type %in%
                                names(table(meta$Cell_Type)[table(meta$Cell_Type) > 1000]),])

sobj <- subset(sobj, cells = cells_over1k)


DefaultAssay(sobj) <- "ATAC"


sobj <- LinkPeaks(
  object = sobj,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)



saveRDS(sobj, "output/objects/atlas_peaks_linked_ct_over1k.RDS")


astros<-readRDS("output/objects/RGL_Astro_wnn_Subtype.RDS")

subtype_meta <- astros@meta.data

sobj@meta.data[rownames(subtype_meta),"Subtype"] <- as.character(subtype_meta$Subtype_updated)

astros <- subset(sobj, Subtype %in% c( "CD44 hi Astrocyte","EAAT2 hi Astrocyte"))

astros <- LinkPeaks(
  object = astros,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)




saveRDS(astros, "output/objects/astros_linked.RDS")

