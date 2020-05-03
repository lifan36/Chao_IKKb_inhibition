### example Seurat code for initial single cell RNA seq analysis ###
### written by Lay Kodama, Bang Liu, and Li Fan ###
### please reference https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html for details on the Seurat package parameters ###

#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/data_analysis")

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG343_78.counts <- Read10X(data.dir = "/Users/lifan/Box Sync/Li_Fan_Research/10x/ChaoWang_IKKb/LG343_78/filtered_feature_bc_matrix")
LG343_78 <- CreateSeuratObject(counts = LG343_78.counts, project = "LG343_78", min.cells = 3, min.features = 200)
LG343_78[["Condition"]] = c('IKKbKO; Cre/+')

#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG343_78[["percent.mt"]] <- PercentageFeatureSet(object = LG343_78, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG343_78

VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)


#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
all
##An object of class Seurat 
##22387 features across 12315 samples within 1 assay 
##Active assay: RNA (22387 features)
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
##An object of class Seurat 
##22092 features across 7500 samples within 1 assay 
##Active assay: RNA (22092 features)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
ElbowPlot(all)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

DimPlot(all, reduction = "umap", label = T)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG343_78_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
saveRDS(all,"LG343_78_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG343_78_PCA.rds")
all

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10267) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_780", sct = FALSE)
#visualizing clusters and multiplet cells====
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
all <- RunTSNE(object = all, dims = 1:15, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = all, reduction = 'umap')

#calculate nearest neighbors
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
#DimPlot(object = all, reduction = 'tsne')
DimPlot(all, reduction = "umap", label = T)

Idents(object = all) <- "DF.classifications_0.25_0.005_780" #visualizing the singlet vs doublet cells
DimPlot(object = all, reduction = 'umap')

all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
singlets<-saveRDS(singlets,"LG343_78_singlets.rds")
singlets<-readRDS("LG343_78_singlets.rds")

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
ElbowPlot(singlets,ndims=50)
singlets <- RunTSNE(object = singlets, dims = 1:15, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = singlets, reduction = 'umap')

#calculate nearest neighbors
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
DimPlot(object = singlets, reduction = 'umap')

#finding markers of your clusters and annotate cells====
singlets.markers <- FindAllMarkers(object = singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(singlets.markers,"LG343_78_cluster_markers_singlets.csv")

singlets<-saveRDS(singlets,"LG343_78_singlets_PCA.rds")

#annotate cell types based on outside databases - this step is tricky and arbitrary, but crucial to get as correct as possible:
#singlets <- RenameIdents(singlets, `0` = "exN", `1` = "Oligodendrocytes", `2` = "exN", 
#                         `3` = "ihN", `4` = "Astrocytes", `5` = "exN", `6` = "exN", `7` = "exN", '8' = "OPC", '9' = "Microglia", '10' = "Pericyte/EC")

#DimPlot(singlets, reduction = "tsne")

#load in meta data information for your genotypes etc - this step is very specific to your dataset:====
#Column_names_pbmc<-as.data.frame(singlets@meta.data$orig.ident)
#colnames(Column_names_pbmc) <- "Sample"
#sample_info <- read.csv("Sample_key.csv", header=T)
#library_info <- merge(Column_names_pbmc, sample_info, by="Sample")
#write.csv(library_info,"library_info.csv")

#singlets <- AddMetaData(object = singlets, metadata = library_info$Genotype, col.name = "Genotype")
#singlets <- AddMetaData(object =singlets, metadata = library_info$Sample, col.name = "Individual")
#saveRDS(singlets,"singlets_clustered.Rds")
#singlets<-readRDS("singlets_clustered.Rds")

#other functions you can use to analyze your data include:
#proportion of cells by genotype
#plotting tsnes by genotype
#calculating DEGs based on genotype for specific cell types
#subclustering specific cell types for finer clustering information
#please reference the different vignettes on Seurat's lab website: https://satijalab.org/seurat/vignettes.html