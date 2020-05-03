library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_KO_Fan")

#setwd("~/Desktop/data_analysis/Chao_IKKb")

LG343_IKKbKO_integrated <- readRDS("LG343_IKKbKO_integrated.rds")
DefaultAssay(LG343_IKKbKO_integrated) <- 'integrated'

# LG343_IKKbKO_integrated <- NormalizeData(LG343_IKKbKO_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# LG343_IKKbKO_integrated <- FindVariableFeatures(LG343_IKKbKO_integrated, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(LG343_IKKbKO_integrated)
LG343_IKKbKO_integrated <- ScaleData(LG343_IKKbKO_integrated, features = all.genes)
LG343_IKKbKO_integrated <- RunPCA(LG343_IKKbKO_integrated, features = VariableFeatures(object = LG343_IKKbKO_integrated))
DimPlot(LG343_IKKbKO_integrated, reduction = 'pca')
LG343_IKKbKO_integrated <- FindNeighbors(LG343_IKKbKO_integrated, dims = 1:20)
LG343_IKKbKO_integrated <- FindClusters(LG343_IKKbKO_integrated, resolution = 0.1)
LG343_IKKbKO_integrated <- RunUMAP(LG343_IKKbKO_integrated, dims = 1: 20, perplexity = 20)


DefaultAssay(LG343_IKKbKO_integrated) <- 'RNA'

pdf("LG343_IKKbKO_integrated_umap.pdf", width=4, height=4)
DimPlot(LG343_IKKbKO_integrated, reduction = 'umap', label = T)
dev.off()
pdf("LG343_IKKbKO_integrated_umap_split.pdf", width=20, height=4)
DimPlot(LG343_IKKbKO_integrated, reduction = "umap", split.by = "orig.ident", label = T)
DimPlot(LG343_IKKbKO_integrated, reduction = "umap", split.by = "Condition", label = T)
dev.off()

saveRDS(LG343_IKKbKO_integrated, file = 'LG343_IKKbKO_integrated_PCA_0.1.rds')


