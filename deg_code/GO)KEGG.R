options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
library(clusterProfiler)
library(DOSE)
library(DO.db)
library(org.Mm.eg.db)
library(stringr)

#####
E13_ctx_cluster2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/DEG/mfuzz/mfuzz_2.csv", sep=",")
gene <- E13_ctx_cluster2[,1]
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Mm.eg.db)
#Go classification
#Go enrichment
ego_cc<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.2,
                 qvalueCutoff = 0.2)
ego_bp<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.2,
                 qvalueCutoff = 0.2)
dev.new()
barplot(ego_bp,showCategory = 18),title="The GO_BP enrichment analysis of all DEGs ")+ 
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))

write.csv(ego_bp@result,file = "E13_ctx_cluster2_bp.csv")

options(clusterProfiler.download.method = "wininet")
kk<-enrichKEGG(gene      =gene.df$ENTREZID,
               organism = 'mmu',
               pvalueCutoff = 0.05)

write.csv(kk@result,file = "E13_ctx_cluster2_kegg.csv")

p0 <- read.csv("E:/yaner xuenian/NSCPRO/downstream/DEG/E13_P0/E13_P0_DEG_group_group_results.csv", sep=",")
p3 <- read.csv("E:/yaner xuenian/NSCPRO/downstream/DEG/E13_P3/E13_P3_DEG_group_group_results.csv", sep=",")
p6 <- read.csv("E:/yaner xuenian/NSCPRO/downstream/DEG/E13_P6/E13_P6_DEG_group_group_results.csv", sep=",")

p0_1 <-p0[,c(4,14)]
p3_1 <-p3[,c(4,14)]
p6_1 <-p6[,c(4,14)]

p0_1$passage <- "p0"
p3_1$passage <- "p3"
p6_1$passage <- "p6"
E13_volin <- rbind(p0_1,p3_1,p6_1)
write.csv(E13_volin,file = "E13_volin.csv")
