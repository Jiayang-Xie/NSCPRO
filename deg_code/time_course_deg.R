options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("maSigPro")
BiocManager::install("Mfuzz")
library(BiocManager)
library(DESeq2)
library(maSigPro)
library(pheatmap)
library(limma)
library(Mfuzz)
# GO
library(clusterProfiler)
library(DOSE)
library(DO.db)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
# 导入数据  E13_5
p1_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A1_L4_311X11_sorted.count", sep="\t", col.names = c("gene_id","p0_ctx1"))
p1_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A2_L4_312X12_sorted.count", sep="\t", col.names = c("gene_id","p0_ctx2"))
p3_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX1_L3_350X50_sorted.count", sep="\t", col.names = c("gene_id","p3_ctx1"))
p3_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","p3_ctx2"))
p6_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX1_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx1"))
p6_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx2"))
p6_ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx3"))
p1_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B1_L4_313X13_sorted.count", sep="\t", col.names = c("gene_id","p0_ge1"))
p1_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B2_L4_380X80_sorted.count", sep="\t", col.names = c("gene_id","p0_ge2"))
p1_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C2_L4_382X82_sorted.count", sep="\t", col.names = c("gene_id","p0_mid2"))
p1_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C1_L4_381X81_sorted.count", sep="\t", col.names = c("gene_id","p0_mid1"))
p1_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D1_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","p0_hid1"))
p1_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D2_L4_316X16_sorted.count", sep="\t", col.names = c("gene_id","p0_hid2"))
p3_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-GE1_sorted.count", sep="\t", col.names = c("gene_id","p3_ge1"))
p3_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-GE2_sorted.count", sep="\t", col.names = c("gene_id","p3_ge2"))
p3_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","p3_hid1"))
p3_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-HID2_sorted.count", sep="\t", col.names = c("gene_id","p3_hid2"))
p3_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-MID1_sorted.count", sep="\t", col.names = c("gene_id","p3_mid1"))
p3_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-MID2_sorted.count", sep="\t", col.names = c("gene_id","p3_mid2"))
p6_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE1_sorted.count", sep="\t", col.names = c("gene_id","p6_ge1"))
p6_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE2_sorted.count", sep="\t", col.names = c("gene_id","p6_ge2"))
p6_ge3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE_sorted.count", sep="\t", col.names = c("gene_id","p6_ge3"))
p6_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","p6_hid1"))
p6_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID2_sorted.count", sep="\t", col.names = c("gene_id","p6_hid2"))
p6_hid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID_sorted.count", sep="\t", col.names = c("gene_id","p6_hid3"))
p6_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID1_sorted.count", sep="\t", col.names = c("gene_id","p6_mid1"))
p6_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID2_sorted.count", sep="\t", col.names = c("gene_id","p6_mid2"))
p6_mid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID_sorted.count", sep="\t", col.names = c("gene_id","p6_mid3"))


p1_ctx <- merge(p1_ctx1,p1_ctx2,by ="gene_id")
p3_ctx <- merge(p3_ctx1,p3_ctx2,by ="gene_id")
p6_ctx <- merge(merge(p6_ctx1,p6_ctx2,by ="gene_id"),p6_ctx3,by ="gene_id")
p1_ge <- merge(p1_ge1,p1_ge2,by ="gene_id")
p3_ge <- merge(p3_ge1,p3_ge2,by ="gene_id")
p6_ge <- merge(merge(p6_ge1,p6_ge2,by ="gene_id"),p6_ge3,by ="gene_id")
p1_mid <- merge(p1_mid1,p1_mid2,by ="gene_id")
p3_mid <- merge(p3_mid1,p3_mid2,by ="gene_id")
p6_mid <- merge(merge(p6_mid1,p6_mid2,by ="gene_id"),p6_mid3,by ="gene_id")
p1_hid <- merge(p1_hid1,p1_hid2,by ="gene_id")
bp3_hid <- merge(p3_hid1,p3_hid2,by ="gene_id")
p6_hid <- merge(merge(p6_hid1,p6_hid2,by ="gene_id"),p6_hid3,by ="gene_id")
ctx <-merge(merge(p1_ctx,p3_ctx,by ="gene_id"),p6_ctx,by ="gene_id")
ge <-merge(merge(p1_ge,p3_ge,by ="gene_id"),p6_ge,by ="gene_id")
mid <-merge(merge(p1_mid,p3_mid,by ="gene_id"),p6_mid,by ="gene_id")
hid <-merge(merge(p1_hid,p3_hid,by ="gene_id"),p6_hid,by ="gene_id")
rawcount <-merge(merge(merge(ctx,ge,by ="gene_id"),mid,by ="gene_id"),hid,by ="gene_id")
rawcount_1 <- rawcount[-c(1:5),]
ENSEMBL <- gsub("\\.\\d*", "",rawcount_1$gene_id)
row.names(rawcount_1) <-ENSEMBL
rawcount_2 <- rawcount_1[,-1]
rawcount_3 <- rawcount_2[apply(rawcount_2, 1, sum) > 20 , ] #remove the low express gene


#construct the coldata,3 bio replicates each group
condition <- factor(c(rep(c(rep("p0",2),rep("p3",2),rep("p6",3)),4)),levels = c("p0","p3","p6"))
tissue <- factor(rep(c("ctx","ge","mid","hid"),each = 7),levels = c("ctx","ge","mid","hid"))
batch <- factor(c(rep("1",6),2,rep("1",6),2,rep("1",4),2,rep("1",8),2))# batch effect E13
colData <- data.frame(row.names = colnames(rawcount_3),condition,tissue,batch)

# 利用DESeq2包对counts进行归一化
dds<- DESeqDataSetFromMatrix(rawcount_3,colData,design = ~ batch+ condition + tissue  )
dds<- estimateSizeFactors(dds)
normalized_counts <-counts(dds,normalized=T)

## 另一种方法
# n_counts <- estimateSizeFactorsForMatrix(rawcount_3)
# normalized_counts <-counts(n_counts,normalized=T)
# fn<- function(x){
#   return(x/1.269)
# }
# list1<- apply(rawcount_3,2,FUN = fn)#这里需要完善，可以写一个循环

##构建好归一化矩阵后，可以开始进行时序分析。
Time = c(rep(c(rep(c(1,1,3,3,6,6),4)),2))
ctx = c(rep(c(rep(1,6),rep(0,18)),2))
ge = c(rep(c(rep(0,6),rep(1,6),rep(0,12)),2))
mid = c(rep(c(rep(0,12),rep(1,6),rep(0,6)),2))
hid = c(rep(c(rep(0,18),rep(1,6)),2))
eday = c(rep(c(13,15),each = 24))
design<-cbind(Time,eday,ctx,ge,mid,hid)
rownames(edesign) <-colnames(normalized_counts)

## 构建回归模型表达式,采用多项式回归，参数degree设置多项式使用的次数。过高会导致过拟合
d_model<-make.design.matrix(edesign,degree = 2)

## 第一步运用最小二乘法求得每个自变量得系数值，同时，通过F检验评估回归方程得显著性。

fit <- p.vector(normalized_counts,d_model,Q =0.05,MT.adjust = "BH",counts = TRUE,min.obs =20)

## 通过逐步回归法确定最佳得自变量组合
tstep<- T.fit(fit,step.method = "backward",alfa = 0.05)


## 获得显著基因集
sigs = get.siggenes(tstep,rsq = 0.7, vars = "groups")
sigs_summary <- sigs$summary
suma2Venn(sigs$summary[, c(2:4)])
## 可视化，对各个样本中得表达模式进行聚类
cluster_genes_<- see.genes(sigs$sig.genes$ctx,
          show.fit = T,
          dis = d_model$dis,
          cluster.method = "hclust",
          cluster.data = 1,
          k = 9)


# Mfuzz 进行时间序列表达
library(Mfuzz)

## 构建表达矩阵，需要求各组平均值


normalized_counts_ctx <- normalized_counts[,22:28]

Mfuzz_mean <- transform(normalized_counts_ctx,p0 = (p0_hid1+p0_hid2)/2,
                          p3 = (p3_hid1+p3_hid2)/2,
                          p6 = (p6_hid1+p6_hid2+p6_hid3)/3)

Mfuzz_matrix  <-as.matrix( Mfuzz_mean[,c(8:10)])

## 构建ExpressionSet对象，根据标准差去除样本间差异太小的基因
eset <- new("ExpressionSet",exprs = Mfuzz_matrix)
eset <- na.omit(eset)
eset <- filter.std(eset,min.std = 0)

## 标准化
eset <- standardise(eset)

## 聚类
m= mestimate(eset)
cl <- mfuzz(eset,c = 9, m= m)
cl$size

## 可视化
pdf("e13_hid_mfuzz.pdf")
mfuzz.plot2(eset,cl,mfrow = c(3,3),centre = TRUE, x11 = F,centre.lwd=0.2, 
            time.labels = seq(0, 42, 21),xlab="Time")
dev.off()

###
dir.create(path="mfuzz_E13_hid",recursive = TRUE)
for(i in 1:9){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("mfuzz_E13_hid","/mfuzz_",i,".csv"))
}

# go analysis

cluster2 <-names(cl$cluster[unname(cl$cluster)==2])

gene.df<-bitr(cluster2, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Mm.eg.db)


ego_bp<-enrichGO(gene       = gene.df$ENSEMBL,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENSEMBL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.2,
                 qvalueCutoff = 0.2)
barplot(ego_bp,showCategory = 18),title="The GO_BP enrichment analysis of e13_CTX_cluster2.pdf ")+ 
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 25))




