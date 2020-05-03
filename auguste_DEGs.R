library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_KO_Fan")

LG343_IKKbKO_integrated <- readRDS("LG343_IKKbKO_integrated_ready_4_DEGs.rds")
LG343_IKKbKO_integrated$celltype.orig.ident <- paste(Idents(LG343_IKKbKO_integrated), LG343_IKKbKO_integrated$orig.ident, sep = "_")
LG343_IKKbKO_integrated$celltype <- Idents(LG343_IKKbKO_integrated)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_KO_Fan/DEGs")

#pdf("UMAP_0.1.pdf", width=8, height=6)
#DimPlot(LG343_IKKbKO_integrated, label = TRUE)
#dev.off()
#pdf("UMAP_0.1_by_condition.pdf", width=16, height=4)
#DimPlot(LG343_IKKbKO_integrated, split.by = "condition", label = TRUE)
#dev.off()

#Subset Neurons: EN and IN
Cluster_EN_IN <- subset(LG343_IKKbKO_integrated, idents = c("excitatory neurons", "inhibitory neurons"))
#Subset Microglia
Cluster_MG <- subset(LG343_IKKbKO_integrated, idents = "microglia")
Idents(Cluster_MG) <- "orig.ident"
avg.Cluster_MG <- log1p(AverageExpression(Cluster_MG, verbose = FALSE)$RNA)
avg.Cluster_MG$gene <- rownames(avg.Cluster_MG)

#Test Cluster_EN_IN and Cluster_MG
#markers.to.plot <- c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
#                     "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
#                     "Bnc2", "Slc47a1", "Ttr")

#pdf("Verify_Subset_neurons_MG.pdf", width=12, height=8)
#DotPlot(Cluster_EN_IN, features = rev(markers.to.plot), cols = c("orange", "plum", "pink", "grey"), dot.scale = 8, 
#        split.by = "condition") + RotatedAxis()
#DotPlot(Cluster_MG, features = rev(markers.to.plot), cols = c("orange", "plum", "pink", "grey"), dot.scale = 8, 
#        split.by = "condition") + RotatedAxis()
#dev.off()

Idents(Cluster_EN_IN) <- "Condition"
Idents(Cluster_MG) <- "Condition"

#IKKb+/+ vs IKKb-/- DEGs
IKKbCtrl_KO_EN_IN_DEGs <- FindMarkers(Cluster_EN_IN, ident.1 = "IKKbKO; Ctrl", ident.2 = "IKKbKO; Cre/+", logfc.threshold = 0,
                                      test.use = "MAST")
write.csv(IKKbCtrl_KO_EN_IN_DEGs, "IKKbCtrl_KO_EN_IN_DEGs.csv")

IKKbCtrl_KO_MG_DEGs <- FindMarkers(Cluster_MG, ident.1 = "IKKbKO; Ctrl", ident.2 = "IKKbKO; Cre/+", logfc.threshold = 0,
                                      test.use = "MAST")
write.csv(IKKbCtrl_KO_MG_DEGs, "IKKbCtrl_KO_MG_DEGs.csv")

#IKKb+/+ vs IKKb+/+;PS301S+  DEGs
IKKbCtrl_PS301S_EN_IN_DEGs <- FindMarkers(Cluster_EN_IN, ident.1 = "IKKbKO; Ctrl", ident.2 = "IKKbKO; Ctrl; PS19/+", logfc.threshold = 0,
                                     test.use = "MAST")
write.csv(IKKbCtrl_PS301S_EN_IN_DEGs, "IKKbCtrl_PS301S_EN_IN_DEGs.csv")

IKKbCtrl_PS301S_MG_DEGs <- FindMarkers(Cluster_MG, ident.1 = "IKKbKO; Ctrl", ident.2 = "IKKbKO; Ctrl; PS19/+", logfc.threshold = 0,
                                   test.use = "MAST")
write.csv(IKKbCtrl_PS301S_MG_DEGs, "IKKbCtrl_PS301S_MG_DEGs.csv")

#IKKb+/+ vs IKKb-/-;PS301S+  DEGs
IKKbCtrl_KOPS301S_EN_IN_DEGs <- FindMarkers(Cluster_EN_IN, ident.1 = "IKKbKO; Ctrl", ident.2 = "IKKbKO; Cre/+; PS19/+", logfc.threshold = 0,
                                          test.use = "MAST")
write.csv(IKKbCtrl_KOPS301S_EN_IN_DEGs, "IKKbCtrl_KOPS301S_EN_IN_DEGs.csv")

IKKbCtrl_KOPS301S_MG_DEGs <- FindMarkers(Cluster_MG, ident.1 = "IKKbKO; Ctrl", ident.2 = "IKKbKO; Cre/+; PS19/+", logfc.threshold = 0,
                                       test.use = "MAST")
write.csv(IKKbCtrl_KOPS301S_MG_DEGs, "IKKbCtrl_KOPS301S_MG_DEGs.csv")















