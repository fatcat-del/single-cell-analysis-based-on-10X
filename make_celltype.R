library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(SingleR)
library(celldex)

pbmc <- readRDS("E:/project/data/singlecell/Q2_douletfinder/pbmc.clr.umap.rds")
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- pbmc[,pbmc@meta.data$seurat_clusters %in% c(1,3,9,11,14,15,16,18,19)]
pbmc <- NormalizeData(pbmc, normalization.method = "CLR")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = F)
pbmc <- RunUMAP(pbmc, min.dist = 0.5,dims = 1:10, seed.use = 0)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

check_gene <- c("Ptprc", "Thy1")
DotPlot(pbmc, group.by = 'seurat_clusters',features = unique(check_gene)) + RotatedAxis()
DimPlot(pbmc, reduction = 'umap', group.by = 'seurat_clusters',label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(pbmc, features = "Ptprc", label = T)
FeaturePlot(pbmc, features = "Thy1", label = T)

Baron <- CreateSeuratObject(counts = pbmc@assays$RNA@counts, project = "Baron", min.cells = 3)
data("HgProteinCodingGenes")
BaronMatrixProt <- BaronMatrix[rownames(BaronMatrix) %in% HgProteinCodingGenes,]
Baron <- CreateSeuratObject(counts = BaronMatrixProt, project = "Baron", min.cells = 3)
Baron <- NormalizeData(Baron)
Baron <- ScaleData(Baron, features = rownames(Baron)) 
Baron <- RunMCA(Baron)

panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
panglao_pancreas <- panglao %>% filter(organ == "Immune system")
panglao_pancreas <- panglao_pancreas %>%  filter(str_detect(species,"Mm"))

panglao_pancreas <- panglao_pancreas %>%  
  group_by(`cell type`) %>%  
  summarise(geneset = list(`official gene symbol`))
pancreas_gs <- setNames(panglao_pancreas$geneset, panglao_pancreas$`cell type`)

HGT_pancreas_gs <- RunCellHGT(Baron, pathways = pancreas_gs, dims = 1:50, n.features = 200)
pancreas_gs_prediction <- rownames(HGT_pancreas_gs)[apply(HGT_pancreas_gs, 2, which.max)]

pancreas_gs_prediction_signif <- ifelse(apply(HGT_pancreas_gs, 2, max)>2, 
                                        yes = pancreas_gs_prediction, 
                                        "unassigned")
Baron$pancreas_gs_prediction <- pancreas_gs_prediction_signif



ref <- ImmGenData()
pred.pbmc <- SingleR(test = pbmc@assays$RNA@data, ref = ref,labels = ref$label.main, 
                     clusters = pbmc@active.ident, fine.tune = TRUE)
plotScoreHeatmap(pred.pbmc, clusters=pred.pbmc@rownames, fontsize.row = 9,show_colnames = T)
new.cluster.ids <- pred.pbmc$pruned.labels
names(new.cluster.ids) <- levels(pbmc)
scpbmc <- RenameIdents(pbmc,new.cluster.ids)
DimPlot(scpbmc, reduction = "umap",pt.size=0.5,label = TRUE)
