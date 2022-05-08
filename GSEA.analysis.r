library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

genelist_input <- fread(file="all.df.gene.8m.csv", header = T, sep=',', data.table = F)
colnames(genelist_input)[1] <- "Gene"
genename <- as.character(genelist_input[,1])
gene_map <- select(org.Mm.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
write.csv(as.data.frame(gene_map),"ID.transfer.csv",row.names =F)


aaa<-inner_join(gene_map,genelist_input,by = "Gene")
aaa<-na.omit(aaa)
aaa <- aaa[order(aaa$avg_log2FC,decreasing = TRUE),]

geneList <- aaa$avg_log2FC
names(geneList) <- as.character(aaa$ENTREZID)
Go_gseresult <- gseGO(geneList, 'org.Mm.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000,pvalueCutoff = 1)
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1, organism = "mmu")
rectomepath <-enrichPathway(gene=aaa$ENTREZID,pvalueCutoff=0.05, readable=T, qvalueCutoff = 0.05, organism = "mouse")




