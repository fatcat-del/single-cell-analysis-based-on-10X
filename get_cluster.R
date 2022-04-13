library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(harmony)

dir.10x <- "E:/project/data/singlecell/GSE189175/10X/"
pbmc.data <- Read10X(data.dir = dir.10x)
pbmc <- pbmc.data[-57537,]
pbmc.data <- NULL
pbmc <- CreateSeuratObject(counts = pbmc, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 550 )
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc, normalization.method = "CLR")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = F)
ElbowPlot(pbmc, ndims = 50)
pbmc <- RunUMAP(pbmc,, min.dist = 0.5 dims = 1:30, seed.use = 0)
UMAPPlot(object = pbmc,label=T, pt.size=1)

pbmc <- RunTSNE(pbmc,dims = 1:20, seed.use = 0)
TSNEPlot(object = pbmc,label=T, pt.size=1)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

check_gene <- c("PTPRC", "ALB", "KRT18", "KRT19", "AFP", "PECAM1", "CD34", "SOX17", "NGFR")
DotPlot(pbmc, group.by = 'seurat_clusters',features = unique(check_gene)) + RotatedAxis()

DimPlot(pbmc, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE, pt.size = 0.5) + NoLegend()
