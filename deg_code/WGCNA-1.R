###################### load the relevant packages################################
# install and load the relevant packages
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("factoextra")
BiocManager::install("FactoMineR")
BiocManager::install("hcluster")
BiocManager::install("tidyverse")
BiocManager::install("ggrepel")
BiocManager::install("VennDiagram")
BiocManager::install("ggVennDiagram")
BiocManager::install("ggvenn")
library(VennDiagram)
library(ggVennDiagram)
library(BiocManager)
library(DESeq2)
library(biomaRt)
library(curl)
library(pheatmap)
library(GEOquery)
library(ggplot2)
library("FactoMineR")
library("factoextra") 
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggvenn)
library(magrittr)
library(limma)
library(Mfuzz)
library(dplyr)
###################### load the data and construct the count matrix and design matrix#########
# input the rawcount files to build the matrix
E15_p0_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0A1_L3_307X07_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ctx1"))
E15_p0_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0A2_L3_308X08_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ctx2"))
E15_p0_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0B1_L3_309X09_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ge1"))
E15_p0_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0B2_L4_310X10_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ge2"))
E15_p0_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0C1_L4_310X10_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_mid2"))
E15_p0_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0C2_L4_311X11_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_mid1"))
E15_p0_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0D1_L4_312X12_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_hid1"))
E15_p0_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0D2_L4_313X13_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_hid2"))
E15_p3_ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX_L4_150A50_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ctx3"))
E15_p3_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX2_L4_346X46_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ctx2"))
E15_p3_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX1_L3_343X43_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ctx1"))
E15_p3_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE1_L3_344X44_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ge1"))
E15_p3_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE2_L4_347X47_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ge2"))
E15_p3_ge3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ge3"))
E15_p3_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID1_L4_351X51_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_mid1"))
E15_p3_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID2_L4_352X52_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_mid2"))
E15_p3_mid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID_L4_152A52_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_mid3"))
E15_p3_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID1_L3_345X45_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_hid1"))
E15_p3_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID2_L3_349X49_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_hid2"))
E15_p3_hid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID_L4_153A53_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_hid3"))
E15_p6_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-CTX1_L4_346X46_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ctx1"))
E15_p6_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-CTX2_L4_349X49_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ctx2"))
E15_p6_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-GE1_L4_347X47_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ge1"))
E15_p6_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-GE2_L4_350X50_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ge2"))
E15_p6_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_hid1"))
E15_p6_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-hid_L4_379X79_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_hid2"))
E15_p6_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-mid_L4_378X78_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_mid1"))
E15_p6_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-mid2.count", sep="\t", col.names = c("gene_id","E15_p6_mid2"))
E15_p0_ctx <- merge(E15_p0_ctx1,E15_p0_ctx2,by ="gene_id")
E15_p3_ctx <-  merge(merge(E15_p3_ctx1,E15_p3_ctx2,by ="gene_id"),E15_p3_ctx3,by ="gene_id")
E15_p6_ctx <-merge(E15_p6_ctx1,E15_p6_ctx2,by ="gene_id")
E15_p0_ge <- merge(E15_p0_ge1,E15_p0_ge2,by ="gene_id")
E15_p3_ge <- merge(merge(E15_p3_ge1,E15_p3_ge2,by ="gene_id"),E15_p3_ge3,by ="gene_id")
E15_p6_ge <- merge(E15_p6_ge1,E15_p6_ge2,by ="gene_id")
E15_p0_mid <- merge(E15_p0_mid1,E15_p0_mid2,by ="gene_id")
E15_p3_mid <- merge(merge(E15_p3_mid1,E15_p3_mid2,by ="gene_id"),E15_p3_mid3,by ="gene_id")
E15_p6_mid <- merge(E15_p6_mid1,E15_p6_mid2,by ="gene_id")
E15_p0_hid <- merge(E15_p0_hid1,E15_p0_hid2,by ="gene_id")
E15_p3_hid <-  merge(merge(E15_p3_hid1,E15_p3_hid2,by ="gene_id"),E15_p3_hid3,by ="gene_id")
E15_p6_hid <-merge(E15_p6_hid1,E15_p6_hid2,by ="gene_id")
E15_ctx <-merge(merge(E15_p0_ctx,E15_p3_ctx,by ="gene_id"),E15_p6_ctx,by ="gene_id")
E15_ge <-merge(merge(E15_p0_ge,E15_p3_ge,by ="gene_id"),E15_p6_ge,by ="gene_id")
E15_mid <-merge(merge(E15_p0_mid,E15_p3_mid,by ="gene_id"),E15_p6_mid,by ="gene_id")
E15_hid <-merge(merge(E15_p0_hid,E15_p3_hid,by ="gene_id"),E15_p6_hid,by ="gene_id")
rawcount_E15 <-merge(merge(merge(E15_ctx,E15_ge,by ="gene_id"),E15_mid,by ="gene_id"),E15_hid,by ="gene_id")


E13_p0_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A1_L4_311X11_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ctx1"))
E13_p0_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A2_L4_312X12_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ctx2"))
E13_p3_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX1_L3_350X50_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ctx1"))
E13_p3_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ctx2"))
E13_p6_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ctx1"))
E13_p6_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ctx2"))
E13_p6_ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ctx3"))
E13_p0_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B1_L4_313X13_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ge1"))
E13_p0_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B2_L4_380X80_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ge2"))
E13_p0_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C2_L4_382X82_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_mid2"))
E13_p0_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C1_L4_381X81_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_mid1"))
E13_p0_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D1_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_hid1"))
E13_p0_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D2_L4_316X16_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_hid2"))
E13_p3_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-GE1_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ge1"))
E13_p3_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-GE2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ge2"))
E13_p3_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_hid1"))
E13_p3_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-HID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_hid2"))
E13_p3_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-MID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_mid1"))
E13_p3_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-MID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_mid2"))
E13_p6_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ge1"))
E13_p6_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ge2"))
E13_p6_ge3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ge3"))
E13_p6_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_hid1"))
E13_p6_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_hid2"))
E13_p6_hid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_hid3"))
E13_p6_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_mid1"))
E13_p6_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_mid2"))
E13_p6_mid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_mid3"))
E13_p0_ctx <- merge(E13_p0_ctx1,E13_p0_ctx2,by ="gene_id")
E13_p6_ctx <-  merge(merge(E13_p6_ctx1,E13_p6_ctx2,by ="gene_id"),E13_p6_ctx3,by ="gene_id")
E13_p3_ctx <-merge(E13_p3_ctx1,E13_p3_ctx2,by ="gene_id")
E13_p0_ge <- merge(E13_p0_ge1,E13_p0_ge2,by ="gene_id")
E13_p6_ge <- merge(merge(E13_p6_ge1,E13_p6_ge2,by ="gene_id"),E13_p6_ge3,by ="gene_id")
E13_p3_ge <- merge(E13_p3_ge1,E13_p3_ge2,by ="gene_id")
E13_p0_mid <- merge(E13_p0_mid1,E13_p0_mid2,by ="gene_id")
E13_p6_mid <- merge(merge(E13_p6_mid1,E13_p6_mid2,by ="gene_id"),E13_p6_mid3,by ="gene_id")
E13_p3_mid <- merge(E13_p3_mid1,E13_p3_mid2,by ="gene_id")
E13_p0_hid <- merge(E13_p0_hid1,E13_p0_hid2,by ="gene_id")
E13_p6_hid <-  merge(merge(E13_p6_hid1,E13_p6_hid2,by ="gene_id"),E13_p6_hid3,by ="gene_id")
E13_p3_hid <-merge(E13_p3_hid1,E13_p3_hid2,by ="gene_id")
E13_ctx <-merge(merge(E13_p0_ctx,E13_p3_ctx,by ="gene_id"),E13_p6_ctx,by ="gene_id")
E13_ge <-merge(merge(E13_p0_ge,E13_p3_ge,by ="gene_id"),E13_p6_ge,by ="gene_id")
E13_mid <-merge(merge(E13_p0_mid,E13_p3_mid,by ="gene_id"),E13_p6_mid,by ="gene_id")
E13_hid <-merge(merge(E13_p0_hid,E13_p3_hid,by ="gene_id"),E13_p6_hid,by ="gene_id")
rawcount_E13 <-merge(merge(merge(E13_ctx,E13_ge,by ="gene_id"),E13_mid,by ="gene_id"),E13_hid,by ="gene_id")

rawcount_0 <- merge(rawcount_E13,rawcount_E15,by ="gene_id")
rawcount_1 <- rawcount_0[-c(1:5),]
ENSEMBL <- gsub("\\.\\d*", "",rawcount_1$gene_id)
row.names(rawcount_1) <-ENSEMBL
rawcount_2 <- rawcount_1[,-1]
rawcount_3 <- rawcount_2[apply(rawcount_2, 1, sum) > 20 , ] 
rawcount_3_n <- rawcount_3 
rawcount_3_n$ensembl_gene_id <- rownames(rawcount_3_n)

write.csv(rawcount_2,file ="all_samples_counts.csv")
############################# visualization of the specific genes ##################
##transverse the gene id to name
## get annotation of the all genes(convert the ensemble_id to gene_symbol) for the deg loop annotation
mart <-useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <-row.names(rawcount_3_n)
mms_symbols <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","description","band"),
                     filters = "ensembl_gene_id", values = my_ensembl_gene_id, mart = mart)
rownames(mms_symbols)<- mms_symbols$ensembl_gene_id
rawcount_3_n1 <- merge(rawcount_3_n,mms_symbols[,1:2],by ="ensembl_gene_id")
rawcount_3_n1 <- rawcount_3_n1 %>% distinct(external_gene_name, .keep_all = TRUE)##因为有重复的gene name，要去掉
row.names(rawcount_3_n1) <- rawcount_3_n1$external_gene_name
rawcount_3_n1 <- rawcount_3_n1[,c(-1,-58)]
write.csv(rawcount_2,file ="all_samples_counts_gene_name.csv")
rawcount_2 <- rawcount_2[,-c(7,14,19,28,33,40,47,54)]
rawcount_3 <- rawcount_2[apply(rawcount_2, 1, sum) > 24 , ] #remove the low express gene
##construct the deseq2 component

#construct the coldata,3 bio replicates each group
condition <- factor(c(rep(c(rep(c("p0","p3","p6"),each= 2),"p6"),4),
                      rep(c("p0","p0","p3","p3","p3","p6","p6"),4)),levels = c("p0","p3","p6"))
tissue <- factor(c(rep(c(rep(c("ctx","ge","mid","hid"),each = 7)),2)),levels = c("ctx","ge","mid","hid"))
stage <- factor(c(rep(c("13","15"),each = 28)))
batch <- factor(c(rep("1",6),2,rep("1",6),2,rep("1",4),2,rep("1",8),2,
                  rep(c(1,1,1,1,2,1,1),4)))# batch effect E13
colData <- data.frame(row.names = colnames(rawcount_3_n1),condition,tissue,stage,batch)
write.csv(colData,file ="all_samples_traits.csv")
# p0,p3,p6
E13_ctx =c(1,2,3,4,6,7)
E13_ge =E13_ctx+7
E13_mid =E13_ctx+14
E13_hid =E13_ctx+21
count_loop<- data.frame( E13_ctx,E13_ge,E13_mid,E13_hid)
rawcount_cl <- rawcount_3[,count_loop[,1]]
rawcount_log_cl <-log2(rawcount_cl+1)
colData_cl <- colData[count_loop[,1],c(1,3)]
## tissue
E13_P0 =c(1,2,8,9,15,16,22,23)
E13_P3 =E13_P0+2
E13_P6 =c(5,6,12,13,20,21,26,27)
count_loop<- data.frame( E13_P0,E13_P3,E13_P6)
rawcount_cl <- rawcount_3_n1[,E13_P0]
rawcount_cl <- rawcount_cl[apply(rawcount_cl, 1, sum) > 8 , ]
rawcount_log_cl <-log2(rawcount_cl+1)
colData_cl <- data.frame( colData[E13_P0,2])
colnames (colData_cl)<-"tissue"

## constrcut the regional marker gene list
nsc  <- c("SOX1","FABP7","DLL1","HES5","NOTCH1","NES","SOX2")
cortex <- c("FOXG1","SP8","EMX1","PAX6","SFRP1","LHX2","SIX3","EOMES","NEUROD6","TBR1","MEF2C")
ge <- c("reln","npy","sst","pde1a","sox5","syt2","cplx1","pvalb")
mid_hind <- c("PAX3","DMBX1","DBX1","EN2","HOXA2","HOXB2","HOXB3","PAX7","PHOX2B","WNT1","BARHL1","ATOH1")


dds <- DESeqDataSetFromMatrix(rawcount_cl,colData_cl,design = ~ tissue )
dds<- DESeq(dds)
summary(dds)
dds
tmp <- plotCounts(dds, gene=c("Cdc45","Gnai3"), intgroup="tissue", returnData=T)
  ggplot(tmp,aes(x=tissue, y=count)) + geom_boxplot(aes(fill=tissue)) + scale_y_log10() + ggtitle("Gnai3")


normalized_counts <-counts(dds,normalized=T)
vstdata <-vst(dds,blind = FALSE)
vstdata_1 <- assay(vstdata)
vstdata_mad <- apply(vstdata_1, 1, mean)
vstdata_2 <- vstdata_1[order(vstdata_mad, decreasing=T), ]
vstdata_3 <- vstdata_2[1:5000, ]


vstdata_2 <- t(vstdata_1)

Time = c(rep(c(rep(c(1,1,3,3,6,6),4)),2))
ctx = c(rep(c(rep(1,6),rep(0,18)),2))
ge = c(rep(c(rep(0,6),rep(1,6),rep(0,12)),2))
mid = c(rep(c(rep(0,12),rep(1,6),rep(0,6)),2))
hid = c(rep(c(rep(0,18),rep(1,6)),2))
eday = c(rep(c(13,15),each = 24))
design<-cbind(Time,eday,ctx,ge,mid,hid)
rownames(design) <-colnames(vstdata_1)

write.table(vstdata_1,"wgcna_rnaseq_exp_1.txt",sep = " ",quote = F)
write.table(design,"wgcna_rnaseq_trait_1.txt",sep = " ",quote = F)
write.csv(vstdata_3,file ="wgcna_rnaseq_exp__5000_mean_2.csv")
write.csv(design,file = "wgcna_rnaseq_trait_1.csv")

all_exp = read.csv(file = "E:/yaner xuenian/NSCPRO/downstream/DEG/wgcna_rnaseq_exp_1.csv",header = T)

### p0 heatmap
rawcount_E13 <-merge(merge(merge(E13_p0_ctx,E13_p0_ge,by ="gene_id"),E13_p0_mid,by ="gene_id"),E13_p0_hid,by ="gene_id")
rawcount_E15 <-merge(merge(merge(E15_p0_ctx,E15_p0_ge,by ="gene_id"),E15_p0_mid,by ="gene_id"),E15_p0_hid,by ="gene_id")
p0_rawcount <- merge(rawcount_E13,rawcount_E15,by = "gene_id")
rawcount_0 <- merge(rawcount_E13,rawcount_E15,by ="gene_id")
rawcount_1 <- rawcount_0[-c(1:5),]
ENSEMBL <- gsub("\\.\\d*", "",rawcount_1$gene_id)
row.names(rawcount_1) <-ENSEMBL
rawcount_2 <- rawcount_1[,-1]
rawcount_3 <- rawcount_2[apply(rawcount_2, 1, sum) > 20 , ] #remove the low express gene
rawcount_log <-log2(rawcount_3+1)
rawcount_log<- na.omit(rawcount_log)
##mad transvert
mad_counts_log <- apply(rawcount_log, 1, mad)
mad_counts_log <- rawcount_log[order(mad_counts_log, decreasing=T), ]#for count-heatmap

##

condition <- factor(c(rep(c(rep(c("p0","p3","p6"),each= 2)),8)),levels = c("p0","p3","p6"))
tissue <- factor(c(rep(c(rep(c("ctx","ge","mid","hid"),each = 2)),2)),levels = c("ctx","ge","mid","hid"))
stage <- factor(c(rep(c("E13","E15"),each = 8)))
colData <- data.frame(row.names = colnames(rawcount_3),tissue,stage)


dds <- DESeqDataSetFromMatrix(rawcount_3,colData,design = ~ stage+ tissue  )
dds<- DESeq(dds)

normalized_counts <-counts(dds,normalized=T)
normalized_counts<- na.omit(normalized_counts)

pheatmap(mad_counts_log[1:1000,],
         cluster_rows = T, 
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         fontsize_col = 12,
         cellwidth = 40,
         annotation_col = colData,
         filename = "P0_E13_E15_mad_1000.png"
         #annotation_col =df
)