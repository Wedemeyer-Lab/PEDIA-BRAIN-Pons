library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

cellcolor = c("ROBO2 Ependymal" = "#000000","SLIT2 Ependymal" = "#4d4d4d","Radial Glia" =  "#fcff5d", "Transit Amplifying Cell" = "#DFFF00","PAX3 Neuroblast" = "#7dfc00", "PAX3 Neuron" ="#90EE90", "HOXB3 Neuroblast" = "#0ec434","Motor Neuron" = "#006400","OPC" = "#0096FF", "Fetal COP" = "#00FFFF", "COP" = "#A7C7E7","Early OL" = "#3998f5","OL" = "#1F51FF", "Astrocyte" = "orange", "Mesenchymal" = "#b732cc", "Immune" = "#FF1493", "Endothelial" = "#000080",  "Erythroblast" = "#37294f", "Pericyte" = "#6E260E")

sobj <- readRDS("output/objects/atlas_final.RDS")

sobj <- FindMultiModalNeighbors(sobj, reduction.list = list("rna_harmony", "atac_harmony"), dims.list = list(1:50, 2:50))
sobj <- RunUMAP(sobj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sobj <- FindClusters(sobj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

saveRDS(sobj, file = "output/objects/atlas.wnn.RDS")

p1 <- DimPlot(sobj, reduction = "umap.rna", group.by = "Cell_Type", label = TRUE, label.size = 2.5, repel = TRUE, cols = cellcolor) + ggtitle("RNA")
p2 <- DimPlot(sobj, reduction = "umap.atac", group.by = "Cell_Type", label = TRUE, label.size = 2.5, repel = TRUE, cols = cellcolor) + ggtitle("ATAC")
p3 <- DimPlot(sobj, reduction = "wnn.umap", group.by = "Cell_Type", label = TRUE, label.size = 2.5, repel = TRUE, cols = cellcolor) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
ggsave("output/figures/WNN_UMAP.pdf", height = 8, width = 16, units="in")