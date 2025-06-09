library(GenomeInfoDb)
library(Matrix)
library(Signac)
library(Seurat)
library(patchwork)
library(tidyverse)
library(dplyr)
library(readr)
library(magrittr)
library(readxl)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020) 
library(TFBSTools)
library(gridExtra)
library(openxlsx)
library(harmony)

cellcolor = c("ROBO2 Ependymal" = "#000000","SLIT2 Ependymal" = "#4d4d4d","Radial Glia" =  "#fcff5d", "Transit Amplifying Cell" = "#DFFF00","PAX3 Neuroblast" = "#7dfc00", "PAX3 Neuron" ="#90EE90", "HOXB3 Neuroblast" = "#0ec434","Motor Neuron" = "#006400","OPC" = "#0096FF", "Fetal COP" = "#00FFFF", "COP" = "#A7C7E7","Early OL" = "#3998f5","OL" = "#1F51FF", "Astrocyte" = "orange", "Mesenchymal" = "#b732cc", "Immune" = "#FF1493", "Endothelial" = "#000080",  "Erythroblast" = "#37294f", "Pericyte" = "#6E260E")

subtypecolor <- c("Ependymal" = "#000000","Fetal Ependymal" = "#36454F", "ROBO2 Radial Glia" = "#e2e700","SLIT2 Radial Glia" = "#FFDAB9", 
                  "Transit Amplifying Cell" = "#DFFF00", "PAX3 Neuroblast" = "#7dfc00", "PAX3 Neuron" ="#90EE90", "HOXB3 Neuroblast" = "#0ec434", 
                  "Glut Neuron" = "#228c68", "Neuron" = "#006400","SEMA3E Astrocyte" = "#FFDAB9","ROBO2 Astrocyte" = "#FF4F00","OL-like Astrocyte" = "#C2A162",
                  "Undetermined Astrocyte" = "#BDC3C7","Fetal COP" = "#00FFFF","RG-like OPC" =  "#89CFF0","OPC" = "#0096FF","COP" = "#A7C7E7",
                  "Astrocyte-like OPC"= "#40E0D0", "Early OL" = "#3998f5","RBFOX1 OL" = "#8A2BE2","OPALIN OL" = "#1A3EAB",
                  "Mesenchymal" = "#b732cc", "Fetal Microglia" = "#FFB3BA", "Microglia" = "#FF6F61", "Lymphocyte" = "darkred",
                  "Lymphatic Endothelial" = "#D5006D", "Immune OL" = "#D8BFD8", "Endothelial" = "#000080",
                  "Erythroblast" = "#37294f", "Pericyte" = "#6E260E")
# dev group 
Development = c("#0D0887", "#7E03A8", "#CC4678", "#F89441", "#F0F921")
names(Development) = c("Fetal", "Infant", "Child",
                       "Adolescent","Adult")
postnatal= c("Postnatal" = "#F0F921", "Fetal" =  "#0D0887")


atlas <- readRDS("output/objects/atlas_final.RDS")
DefaultAssay(atlas) <- "SCT"
atlas@assays$peaksMACS2 <- NULL
atlas@assays$ATAC <- NULL
#no axis no legend
atlas$Cell_Type <- factor(atlas$Cell_Type, levels=names(cellcolor))
DimPlot(subset(atlas, postnatal == "Fetal"), label = FALSE, label.size = 3,
        repel = FALSE, group.by = "Cell_Type", cols = cellcolor) + NoLegend()  +ggtitle(NULL)
DimPlot(subset(atlas, postnatal != "Fetal"), label = FALSE, label.size = 3,
        repel = FALSE, group.by = "Cell_Type", cols = cellcolor) + NoLegend()  +ggtitle(NULL) 
# no legend
DimPlot(atlas, label = FALSE, label.size = 3,
        repel = FALSE, group.by = "Cell_Type", cols = cellcolor) + ggtitle(NULL) +NoLegend()
# no axis
DimPlot(atlas, label = FALSE, label.size = 3,
        repel = FALSE, group.by = "Cell_Type", cols = cellcolor) +ggtitle(NULL) + NoAxes()

# axis and legend
DimPlot(atlas, label = FALSE, label.size = 3,
        repel = FALSE, group.by = "Cell_Type", cols = cellcolor) + ggtitle(NULL) 
