library(Seurat)
library(dplyr)

sobj<-readRDS("output/objects/atlas_peaks_linked.RDS")
sobj@assays$peaksMACS2 <- NULL
gc()


atlas_meta<-sobj@meta.data[,c("sample_id", "Cell_Type","dev.group")]
atlas_meta<-subset(atlas_meta, dev.group != "Fetal")

atlas_meta$condition <- "control"
atlas_meta$dev.group <- NULL


dmgs <- readRDS("output/objects/mapping/dmg_map_to_ATLAS.rds")

dmg_meta<-dmgs@meta.data[,c("sample_name", "Cell_Type", "condition")]
colnames(dmg_meta)[1] <- "sample_id"
dmg_meta$condition <- "tumor"

meta <- rbind(atlas_meta, dmg_meta)

meta<-read.csv("output/tables/mapping/atlas_meta.csv",row.names = 1)
meta$Subtype <- meta$Cell_Type %>% as.character()
OL<-readRDS("output/objects/EarlyOL_OL_wnn_Subtype.RDS")
astro<-readRDS("output/objects/RGL_Astro_wnn_Subtype.RDS")
immune<-readRDS("output/objects/Myeloid_subtype.rds")
OPC <-readRDS("output/objects/OPC_wnn_Subtype.RDS")

astro$Subtype_updated %>% unique()
OL$Subtype_updated  %>% unique()
OPC$Subtype_updated  %>% unique()
immune$Subtype_updated  %>% unique()

meta[rownames(astro@meta.data),"Subtype"] <-astro$Subtype_updated  %>% as.character()
meta[rownames(OL@meta.data),"Subtype"] <-OL$Subtype_updated %>% as.character()
meta[rownames(immune@meta.data),"Subtype"] <-immune$Subtype_updated %>% as.character()
meta[rownames(OPC@meta.data),"Subtype"] <-OPC$Subtype_updated %>% as.character()

meta %>% write.csv("output/atlas_meta.csv")

ms <- readRDS("output/objects/mapping/ms_map_to_HD-BAT.rds")

ms_meta<-ms@meta.data[,c("sample_id", "Subtype")]
ms_meta$condition <- "MS"



atlas_meta<-meta[,c("sample_id", "Subtype","dev.group")]
atlas_meta<-subset(atlas_meta, dev.group != "Fetal")

atlas_meta$condition <- "control"
atlas_meta$dev.group <- NULL

atlas_meta$Subtype[atlas_meta$Subtype =="CD44 hi Astrocyte"]   <-       "CD44-hi Astrocyte" 
atlas_meta$Subtype[atlas_meta$Subtype =="OL like Astrocyte" ]   <-"OL-like Astrocyte" 
atlas_meta$Subtype[atlas_meta$Subtype =="EAAT2 hi Astrocyte"]   <-       "EAAT2-hi Astrocyte"
atlas_meta$Subtype[atlas_meta$Subtype =="Lymphocytes"]   <-       "Lymphocyte"
meta <- rbind(atlas_meta, ms_meta)

write.csv(meta, "output/meta_atlas_ms.csv")

