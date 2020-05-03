#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_KO_Fan")

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:

LG343_79 <- readRDS(file = "LG343_79_singlets_PCA.rds")
LG343_108 <- readRDS(file = "LG343_108_singlets_PCA.rds")
LG343_107 <- readRDS(file = "LG343_107_singlets_PCA.rds")
LG343_126 <- readRDS(file = "LG343_126_singlets_PCA.rds")
LG343_78 <- readRDS(file = "LG343_78_singlets_PCA.rds")
LG343_116 <- readRDS(file = "LG343_116_singlets_PCA.rds")
LG343_80 <- readRDS(file = "LG343_80_singlets_PCA.rds")
LG343_112 <- readRDS(file = "LG343_112_singlets_PCA.rds")

LG343_IKKbKO_combine <- c(LG343_79, LG343_108, LG343_107, LG343_126, LG343_78, LG343_116, LG343_80, LG343_112)

anchors <- FindIntegrationAnchors(object.list = LG343_IKKbKO_combine, dims = 1:30)
LG343_IKKbKO_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

rm(LG343_79, LG343_108, LG343_107, LG343_126, LG343_78, LG343_116, LG343_80, LG343_112, anchors)

#LG343_IKKbKO_integrated[["percent.mt"]] <- PercentageFeatureSet(object = LG343_IKKbKO_integrated, pattern = "^mt-")

pdf("LG343_IKKbKO_merge_VlnPlot.pdf", width=12, height=4)
Idents(LG343_IKKbKO_integrated) <- "orig.ident"
VlnPlot(object = LG343_IKKbKO_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(LG343_IKKbKO_integrated) <- "Condition"
VlnPlot(object = LG343_IKKbKO_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

saveRDS(LG343_IKKbKO_integrated, file = "LG343_IKKbKO_integrated.rds")


