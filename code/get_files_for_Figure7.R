# import packages
library(Matrix)
library(Seurat)
library(Signac)
library(patchwork)
library(tidyverse)
library(dplyr)
library(readr)
library(magrittr)
library(openxlsx) 
library(harmony)
library(clustree)
library(knitr)
set.seed(6)

# load objects
cancer <-readRDS("output/objects/mapping/cancer_map_to_HD-BAT_dev.group.rds")
ms<- readRDS("output/objects/mapping/ms_map_to_HD-BAT_dev.group.rds")
dropseq <- readRDS("output/objects/mapping/dropseq_map_to_HD-BAT_dev_group.rds")

# create object list
sobj.list <- c(cancer, ms, dropseq)
names(sobj.list) <- c("cancer","ms","dropseq")

# save metadata and umap coordinates for each disease state
for (i in 1:3) {
  sobj<-sobj.list[[i]]
  write.csv(sobj@meta.data,paste0("output/tables/disease_mapping/" ,names(sobj.list)[i],"_meta.csv"))
  write.csv(sobj@reductions$ref.umap@cell.embeddings,paste0("output/tables/disease_mapping",names(sobj.list)[i],"_umap_coords.csv"))
}
