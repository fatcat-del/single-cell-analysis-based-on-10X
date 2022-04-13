library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

df <- read.csv("E:/project/data/singlecell/test/net_lr.csv", header = T, row.names = .1)
out.index <- which(df$source == "CD90+ stem cell")
in.index <- which(df$target == "CD90+ stem cell")
outgoing <- df[stem.index,]
incoming <- df[in.index,]

df.gene <- bitr(unique(incoming$receptor), fromType = "SYMBOL",
                  toType = c( "ENTREZID"),
                  OrgDb = org.Mm.eg.db)

df.outgene <- bitr(unique(outgoing$ligand), fromType = "SYMBOL",
                   toType = c( "ENTREZID"),
                   OrgDb = org.Mm.eg.db)

kegg <- enrichKEGG(gene=df.gene$ENTREZID,
                   organism="mmu",
                   pvalueCutoff=0.05,
                   pAdjustMethod="BH",
                   qvalueCutoff=0.1)

barplot(kegg, showCategory=20, font.size=8, title="KEGG pathway")

Go <- enrichGO(gene=df.gene$ENTREZID,
               OrgDb = "org.Mm.eg.db",
               keyType="ENTREZID",
               ont="ALL",
               pvalueCutoff=0.05,
               pAdjustMethod="BH",
               qvalueCutoff=0.1)
barplot(BP,drop =T,showCategory = 5,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scale = 'free')
