#load the libraries
library(Seurat)
library(ggplot2) 
library(tidyverse)
library(devtools)
library(ggpubr)
library(Hmisc)
library(ggpubr)
library(ggrepel)
library(ggfortify)
library(pheatmap)
library(cowplot)
library(patchwork)
library(matrixStats)

# load the single-cell sequencing data
A1.data <- Read10X(data.dir = "/Users/li219/Box/Popovich Lab-Ma Lab collaborative projects/Faith Brennan/CellPhoneDB pilot study/Single cell seq data (sham, 7d, 28 d vehicle and PLX)/1")
A1 <- CreateSeuratObject(A1.data, project = "00 d Control", assay = "RNA", min.cells = 3, min.features = 200)
A2.data <- Read10X(data.dir = "/Users/li219/Box/Popovich Lab-Ma Lab collaborative projects/Faith Brennan/CellPhoneDB pilot study/Single cell seq data (sham, 7d, 28 d vehicle and PLX)/2")
A2 <- CreateSeuratObject(A2.data, project = "00 d PLX5622", assay = "RNA", min.cells = 3, min.features = 200)
A3.data <- Read10X(data.dir = "/Users/li219/Box/Popovich Lab-Ma Lab collaborative projects/Faith Brennan/CellPhoneDB pilot study/Single cell seq data (sham, 7d, 28 d vehicle and PLX)/3")
A3 <- CreateSeuratObject(A3.data, project = "07 d Control", assay = "RNA", min.cells = 3, min.features =200)
A4.data <- Read10X(data.dir = "/Users/li219/Box/Popovich Lab-Ma Lab collaborative projects/Faith Brennan/CellPhoneDB pilot study/Single cell seq data (sham, 7d, 28 d vehicle and PLX)/4")
A4 <- CreateSeuratObject(A4.data, project = "07 d PLX5622", assay = "RNA", min.cells = 3, min.features = 200)
A7.data <- Read10X(data.dir = "/Users/li219/Box/Popovich Lab-Ma Lab collaborative projects/Faith Brennan/CellPhoneDB pilot study/Single cell seq data (sham, 7d, 28 d vehicle and PLX)/7")
A7 <- CreateSeuratObject(A7.data, project = "28 d Control", assay = "RNA", min.cells = 3, min.features =200)
A8.data <- Read10X(data.dir = "/Users/li219/Box/Popovich Lab-Ma Lab collaborative projects/Faith Brennan/CellPhoneDB pilot study/Single cell seq data (sham, 7d, 28 d vehicle and PLX)/8")
A8 <- CreateSeuratObject(A8.data, project = "28 d PLX5622", assay = "RNA", min.cells = 3, min.features =200)

# merge the data
cord.bigthree <- merge(A1, y = c(A2, A3, A4, A7, A8), add.cell.ids = c("00 d Control", "00 d PLX5622", "07 d Control", "07 d PLX5622", "28 d Control", "28 d PLX5622"), project = "all")

# check the data
cord.bigthree
head(colnames(cord.bigthree))
tail(colnames(cord.bigthree))
unique(sapply(X = strsplit(colnames(cord.bigthree), split = "_"), FUN = "[", 1))
table(cord.bigthree$orig.ident)

# perform integration
combinethree.list <- SplitObject(cord.bigthree, split.by = "orig.ident")
combinethree.list <- lapply(X = combinethree.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
combine.anchors <- FindIntegrationAnchors(object.list = combinethree.list, dims = 1:15)
combine.combined <- IntegrateData(anchorset = combine.anchors, dims = 1:15)
DefaultAssay(combine.combined) <- "integrated"
combine.combined <- ScaleData(combine.combined, verbose = FALSE)
combine.combined <- RunPCA(combine.combined, npcs = 30, verbose = FALSE)
combine.combined <- RunUMAP(combine.combined, reduction = "pca", dims = 1:15)
combine.combined <- FindNeighbors(combine.combined, reduction = "pca", dims = 1:15)
combine.combined <- FindClusters(combine.combined, resolution = 0.15)

#Rename cluster identities
combine.combined <- RenameIdents(combine.combined, `13` = "Leptomeningeal cells")
combine.combined <- RenameIdents(combine.combined, `12` = "Pericytes")
combine.combined <- RenameIdents(combine.combined, `11` = "Oligodendrocyte lineage")
combine.combined <- RenameIdents(combine.combined, `10` = "Erythroid cells")
combine.combined <- RenameIdents(combine.combined, `9` = "Neutrophils")
combine.combined <- RenameIdents(combine.combined, `8` = "Intermediate progenitors")
combine.combined <- RenameIdents(combine.combined, `7` = "B cells")
combine.combined <- RenameIdents(combine.combined, `6` = "Ependymal cells")
combine.combined <- RenameIdents(combine.combined, `5` = "T cells")
combine.combined <- RenameIdents(combine.combined, `4` = "Astrocytes") 
combine.combined <- RenameIdents(combine.combined, `3` = "Monocytes")
combine.combined <- RenameIdents(combine.combined, `2` = "Endothelial cells")
combine.combined <- RenameIdents(combine.combined, `1` = "Macrophages")
combine.combined <- RenameIdents(combine.combined, `0` = "Microglia")

table(Idents(combine.combined))

#Export count matrix for cell types of interest for the CellPhone DB package (microglia, macrophages, and astrocytes)
microglia.raw.data <- as.matrix(GetAssayData(combine.combined)[, WhichCells(combine.combined, idents = "Microglia")])
write.csv(microglia.raw.data, "Microglia_raw_data.csv")

macrophage.raw.data <- as.matrix(GetAssayData(combine.combined)[, WhichCells(combine.combined, idents = "Macrophages")])
write.csv(macrophage.raw.data, "Macrophage_raw_data.csv")

astrocyte.raw.data <- as.matrix(GetAssayData(combine.combined)[, WhichCells(combine.combined, idents = "Astrocytes")])
write.csv(astrocyte.raw.data, "Astrocyte_raw_data.csv")

