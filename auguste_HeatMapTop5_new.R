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
#LG343_markers <- FindAllMarkers(LG343_IKKbKO_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
LG343_markers <- read.csv(file = "LG343_markers.csv", header=T,row.names =1)
top5 <- LG343_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top5$gene <- as.character(top5$gene)
pdf("HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(LG343_IKKbKO_integrated, features = top5$gene) + NoLegend()
dev.off()