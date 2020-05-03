library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_KO_Fan")

LG343_IKKbKO_integrated <- readRDS("LG343_IKKbKO_integrated_PCA_0.1.rds")

DefaultAssay(LG343_IKKbKO_integrated) <- 'RNA'
LG343_IKKbKO_integrated <- ScaleData(LG343_IKKbKO_integrated)

LG343_markers <- FindAllMarkers(LG343_IKKbKO_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(LG343_markers, "LG343_markers.csv")

