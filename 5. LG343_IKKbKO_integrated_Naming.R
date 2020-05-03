library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("~/Desktop/data_analysis/Chao_IKKb")

LG343_IKKbKO_integrated <- readRDS("LG343_IKKbKO_integrated_PCA_0.1.rds")

DefaultAssay(LG343_IKKbKO_integrated) <- 'RNA'

DimPlot(LG343_IKKbKO_integrated, reduction = 'umap', label = T)

#Add marker genes
#excitatory neurons
#sig_EN<-c("SLC17A7", "CAMK2A", "NRGN")
sig_EN<-c("Slc17a7", "Camk2a", "Nrgn")
markers.to.plot <- as.matrix(sig_EN)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
DotPlot(object = LG343_IKKbKO_integrated, features = c("Pla2g7", "Bnc2", "Slc47a1", "Ttr", "Trpm3")) + RotatedAxis()
DotPlot(object = LG343_IKKbKO_integrated, features = c("Map2")) + RotatedAxis()
#inhibitory neurons
#sig_EN<-c("GAD1", "GAD2")
sig_IN<-c("Gad1", "Gad2")
markers.to.plot <- as.matrix(sig_IN)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#oligodendrocytes
#sig_EN<-c("PLP1", "MBP", "MOBP")
sig_OL<-c("Plp1", "Mbp", "Mobp")
markers.to.plot <- as.matrix(sig_OL)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#OPCs
#sig_OPC<-c("SCRG1", "PDGFRA", "OLIG1")
sig_OPC<-c("Scrg1", "Pdgfra", "Olig1")
markers.to.plot <- as.matrix(sig_OPC)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#endothelial cells
#sig_ENC<-c("VTN", "MGP", "IGFBP7")
sig_ENC<-c("Vtn", "Mgp", "Igfbp7")
markers.to.plot <- as.matrix(sig_ENC)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#microglia
#sig_MICROG<-c("C1QA", "TYROBP", "HEXB")
sig_MICROG<-c("Cx3cr1", "P2ry12", "Csf1r")
markers.to.plot <- as.matrix(sig_MICROG)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#astrocytes
#sig_AST<-c("CLU, ALDOC, PLA2G7")
sig_AST<-c("Clu", "Aldoc", "Pla2g7")
markers.to.plot <- as.matrix(sig_AST)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#T lymphocytes
#sig_TL<-c("CD3D", "CD3G", "CCL5")
sig_TL<-c("Cd3d", "Cd3g", "Ccl5")
markers.to.plot <- as.matrix(sig_TL)
DotPlot(object = LG343_IKKbKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

LG343_IKKbKO_integrated <- RenameIdents(LG343_IKKbKO_integrated,
        `0` = "oligodendrocytes", `1`="excitatory neurons", `2`="excitatory neurons", `3`="inhibitory neurons",
        `4`="excitatory neurons", `5`="mixed cell type", `6`="excitatory neurons", `7`="astrocytes",
        `8`="inhibitory neurons", `9`="excitatory neurons", `10`="inhibitory neurons", `11`="microglia",
        `12`="OPCs", `13`="inhibitory neurons", `14`="inhibitory neurons", `15`="excitatory neurons",
        `16`="excitatory neurons", `17`="excitatory neurons", `18`="excitatory neurons", `19`="OPCs",
        `20`="endothelial cells", `21`="vascular cells" , `22`="vascular cells", `23`="ependymal cells"
)
DimPlot(LG343_IKKbKO_integrated, label = TRUE)

DimPlot(LG343_IKKbKO_integrated, split.by = "Condition", label = TRUE)

Idents(LG343_IKKbKO_integrated) <- factor(Idents(LG343_IKKbKO_integrated), levels = c(names(table(LG343_IKKbKO_integrated@active.ident))))
markers.to.plot <- c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
                     "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
                     "Bnc2", "Slc47a1", "Ttr")
DotPlot(LG343_IKKbKO_integrated, features = rev(markers.to.plot), cols = c("orange", "plum", "pink", "grey"), dot.scale = 8, 
        split.by = "Condition") + RotatedAxis()

