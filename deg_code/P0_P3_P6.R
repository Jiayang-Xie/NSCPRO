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
###################### load the data and construct the count matrix and design matrix#########
# input the rawcount files to build the matrix
p0_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0A1_L3_307X07_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ctx1"))
p0_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0A2_L3_308X08_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ctx2"))
p0_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0B1_L3_309X09_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ge1"))
p0_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0B2_L4_310X10_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_ge2"))
p0_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0C1_L4_310X10_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_mid2"))
p0_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0C2_L4_311X11_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_mid1"))
p0_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0D1_L4_312X12_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_hid1"))
p0_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0D2_L4_313X13_sorted.count", sep="\t", col.names = c("gene_id","E15_p0_hid2"))
p3_ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX_L4_150A50_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ctx3"))
p3_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX2_L4_346X46_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ctx2"))
p3_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX1_L3_343X43_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ctx1"))
p3_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE1_L3_344X44_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ge1"))
p3_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE2_L4_347X47_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ge2"))
p3_ge3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_ge3"))
p3_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID1_L4_351X51_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_mid1"))
p3_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID2_L4_352X52_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_mid2"))
p3_mid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID_L4_152A52_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_mid3"))
p3_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID1_L3_345X45_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_hid1"))
p3_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID2_L3_349X49_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_hid2"))
p3_hid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID_L4_153A53_sorted.count", sep="\t", col.names = c("gene_id","E15_p3_hid3"))
p6_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-CTX1_L4_346X46_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ctx1"))
p6_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-CTX2_L4_349X49_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ctx2"))
p6_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-GE1_L4_347X47_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ge1"))
p6_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-GE2_L4_350X50_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_ge2"))
p6_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_hid1"))
p6_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-hid_L4_379X79_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_hid2"))
p6_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-mid_L4_378X78_sorted.count", sep="\t", col.names = c("gene_id","E15_p6_mid1"))
p6_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-mid2.count", sep="\t", col.names = c("gene_id","E15_p6_mid2"))
###  E15
p0_ctx <- merge(p0_ctx1,p0_ctx2,by ="gene_id")
p3_ctx <-  merge(merge(p3_ctx1,p3_ctx2,by ="gene_id"),p3_ctx3,by ="gene_id")
p6_ctx <-merge(p6_ctx1,p6_ctx2,by ="gene_id")
p0_ge <- merge(p0_ge1,p0_ge2,by ="gene_id")
p3_ge <- merge(merge(p3_ge1,p3_ge2,by ="gene_id"),p3_ge3,by ="gene_id")
p6_ge <- merge(p6_ge1,p6_ge2,by ="gene_id")
p0_mid <- merge(p0_mid1,p0_mid2,by ="gene_id")
p3_mid <- merge(merge(p3_mid1,p3_mid2,by ="gene_id"),p3_mid3,by ="gene_id")
p6_mid <- merge(p6_mid1,p6_mid2,by ="gene_id")
p0_hid <- merge(p0_hid1,p0_hid2,by ="gene_id")
p3_hid <-  merge(merge(p3_hid1,p3_hid2,by ="gene_id"),p3_hid3,by ="gene_id")
p6_hid <-merge(p6_hid1,p6_hid2,by ="gene_id")
ctx <-merge(merge(p0_ctx,p3_ctx,by ="gene_id"),p6_ctx,by ="gene_id")
ge <-merge(merge(p0_ge,p3_ge,by ="gene_id"),p6_ge,by ="gene_id")
mid <-merge(merge(p0_mid,p3_mid,by ="gene_id"),p6_mid,by ="gene_id")
hid <-merge(merge(p0_hid,p3_hid,by ="gene_id"),p6_hid,by ="gene_id")
rawcount_E15 <-merge(merge(merge(ctx,ge,by ="gene_id"),mid,by ="gene_id"),hid,by ="gene_id")


p0_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A1_L4_311X11_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ctx1"))
p0_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A2_L4_312X12_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ctx2"))
p3_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX1_L3_350X50_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ctx1"))
p3_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ctx2"))
p6_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ctx1"))
p6_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ctx2"))
p6_ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ctx3"))
p0_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B1_L4_313X13_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ge1"))
p0_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B2_L4_380X80_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_ge2"))
p0_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C2_L4_382X82_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_mid2"))
p0_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C1_L4_381X81_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_mid1"))
p0_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D1_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_hid1"))
p0_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D2_L4_316X16_sorted.count", sep="\t", col.names = c("gene_id","E13_p0_hid2"))
p3_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-GE1_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ge1"))
p3_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-GE2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_ge2"))
p3_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_hid1"))
p3_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-HID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_hid2"))
p3_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-MID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_mid1"))
p3_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-MID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p3_mid2"))
p6_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ge1"))
p6_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ge2"))
p6_ge3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-GE_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_ge3"))
p6_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_hid1"))
p6_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_hid2"))
p6_hid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-HID_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_hid3"))
p6_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID1_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_mid1"))
p6_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID2_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_mid2"))
p6_mid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-MID_sorted.count", sep="\t", col.names = c("gene_id","E13_p6_mid3"))
###  E13
p0_ctx <- merge(p0_ctx1,p0_ctx2,by ="gene_id")
p6_ctx <-  merge(merge(p6_ctx1,p6_ctx2,by ="gene_id"),p6_ctx3,by ="gene_id")
p3_ctx <-merge(p3_ctx1,p3_ctx2,by ="gene_id")
p0_ge <- merge(p0_ge1,p0_ge2,by ="gene_id")
p6_ge <- merge(merge(p6_ge1,p6_ge2,by ="gene_id"),p6_ge3,by ="gene_id")
p3_ge <- merge(p3_ge1,p3_ge2,by ="gene_id")
p0_mid <- merge(p0_mid1,p0_mid2,by ="gene_id")
p6_mid <- merge(merge(p6_mid1,p6_mid2,by ="gene_id"),p6_mid3,by ="gene_id")
p3_mid <- merge(p3_mid1,p3_mid2,by ="gene_id")
p0_hid <- merge(p0_hid1,p0_hid2,by ="gene_id")
p6_hid <-  merge(merge(p6_hid1,p6_hid2,by ="gene_id"),p6_hid3,by ="gene_id")
p3_hid <-merge(p3_hid1,p3_hid2,by ="gene_id")
ctx <-merge(merge(p0_ctx,p3_ctx,by ="gene_id"),p6_ctx,by ="gene_id")
ge <-merge(merge(p0_ge,p3_ge,by ="gene_id"),p6_ge,by ="gene_id")
mid <-merge(merge(p0_mid,p3_mid,by ="gene_id"),p6_mid,by ="gene_id")
hid <-merge(merge(p0_hid,p3_hid,by ="gene_id"),p6_hid,by ="gene_id")
rawcount_E13 <-merge(merge(merge(ctx,ge,by ="gene_id"),mid,by ="gene_id"),hid,by ="gene_id")

rawcount_0 <- merge(rawcount_E13,rawcount_E15,by ="gene_id")
rawcount_1 <- rawcount_E13[-c(1:5),]
ENSEMBL <- gsub("\\.\\d*", "",rawcount_1$gene_id)
row.names(rawcount_1) <-ENSEMBL
rawcount_2 <- rawcount_1[,-1]
rawcount_2 <- rawcount_2[,-c(7,14,19,28)]
rawcount_3 <- rawcount_2[apply(rawcount_2, 1, sum) > 10 , ] #remove the low express gene
#construct the coldata,3 bio replicates each group
condition <- factor(c(rep(c(rep("p0",2),rep("p3",2),rep("p6",2)),4)),levels = c("p0","p3","p6"))
tissue <- factor(rep(c("ctx","ge","mid","hid"),each = 6),levels = c("ctx","ge","mid","hid"))
batch <- factor(c(rep(c(rep("1",4),2,rep("1",2)),4)))# batch effect
#batch <- factor(c(rep("1",6),2,rep("1",6),2,rep("1",4),2,rep("1",8),2))# batch effect E13
colData <- data.frame(row.names = colnames(rawcount_3),condition,tissue)







## get annotation of the all genes(convert the ensemble_id to gene_symbol) for the deg loop annotation
mart <-useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <-row.names(rawcount_2)
mms_symbols <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","description","band"),
                     filters = "ensembl_gene_id", values = my_ensembl_gene_id, mart = mart)
rownames(mms_symbols)<- mms_symbols$ensembl_gene_id

#construct the coldata,3 bio replicates each group
condition <- factor(c(rep(c(rep("p0",2),rep("p3",2),rep("p6",3)),4)),levels = c("p0","p3","p6"))
tissue <- factor(rep(c("ctx","ge","mid","hid"),each = 6),levels = c("ctx","ge","mid","hid"))
batch <- factor(c(rep(c(rep("1",4),2,rep("1",2)),4)))# batch effect
#batch <- factor(c(rep("1",6),2,rep("1",6),2,rep("1",4),2,rep("1",8),2))# batch effect E13
colData <- data.frame(row.names = colnames(rawcount_3),condition,tissue,batch)

###################### construct the loop for the analysis  #############################################
E13_ctx =c(1,2,3,4,6,7)
E13_ge =E13_ctx+7
E13_mid =E13_ctx+14
E13_hid =E13_ctx+21
count_loop<- data.frame( E13_ctx,E13_ge,E13_mid,E13_hid)
rawcount_cl <- rawcount_3[,count_loop[,1]]
rawcount_log_cl <-log2(rawcount_cl+1)
colData_cl <- colData[count_loop[,1],c(1,3)]

#remove the batch2 version

#start the DEG workflow
##construct the DESeqDataSet(dds) Object
dds <- DESeqDataSetFromMatrix(rawcount_cl,colData_cl,design = ~condition)
dds<- DESeq(dds)

normalized_counts <- counts(dds,normalized = T) #提取标准化数据
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
normalized_counts_log <- apply(rawcount_log_cl, 1, mad)
normalized_counts_log <- rawcount_log_cl[order(normalized_counts_log, decreasing=T), ]#for count-heatmap
#### the three major plot to evaluate the sample
### 01 PCA by rlog
dev.new()
rld_plot <- rlogTransformation(dds,blind = F)
pca_plot<-plotPCA(rld_plot,intgroup = ("condition"),returnData = T)
plotPCA(rld_plot,intgroup = ("condition"))#package function
png("E15_ctx_PCA.png")


E15_P6_PCA<- plotPCA(rld_plot,intgroup = ("condition"))
print(E15_P6_PCA)
dev.off()

### 02 sample distance heatmap  
sampleDist_data <-dist(t(assay(rld_plot)))
sampleDist_data_matrix <- as.matrix(sampleDist_data)
rownames(sampleDist_data_matrix)<-paste(pca_plot$name,sep = "-")
colnames(sampleDist_data_matrix)<- paste(pca_plot$name,sep = "-")
colors<- colorRampPalette(rev(brewer.pal(9,"GnBu")))(100)
pheatmap_sample<- pheatmap(sampleDist_data_matrix,
                           clustering_distance_rows = sampleDist_data,
                           clustering_distance_cols = sampleDist_data,
                           border_color = NA,
                           cellwidth = 70,
                           fontsize_row = 15,
                           fontsize_col = 15,
                           #height = 100,
                           #width = 100,
                           filename = "E15_ctx_dis.png",
                           col =colors)  
### 03 sample expression heatmap
### count matrix 
#slect <-order(rowMeans(counts(dds,normalized =TRUE)),decreasing = TRUE)[1:1000] #按照每行counts平均数排序，取指定数量
#df <- as.data.frame(colData(dds)[,c("tissue")])#图例注释信息
#ntd <- normTransform(dds) # log2(n+1)
pheatmap(normalized_counts_log[1:1000,],
         cluster_rows = T, 
         show_rownames = FALSE,
         cluster_cols = FALSE,
         fontsize_col = 12,
         cellwidth = 40,
         filename = "E15_ctx_count_heat.png"
         #annotation_col =df
)


###################### the DEG analysis and results extraction#########################
loop_tissue<-data.frame(p3_vs_p0=c("p3","p0"), ##P1,P2,P3更改
                        p6_vs_p3=c("p6","p3"), 
                        p6_vs_p0=c("p6","p0"))
loop_deg <- data.frame(p3_vs_p0=c(1:4),##P1,P2,P3更改
                       p6_vs_p3=c(3:6), 
                       p6_vs_p0=c(1,2,5,6))
top10sig <-c() #for the multi_deg plot
diff_multi<-c()#for the multi_deg plot
gene_list <-list() # for the vennplot
#gene_list <-list()#for the extraction of the intersection of the multi_group deg
for (i in 1:3) {##P1,P2,P3更改
  res <-results(dds,contrast = c("condition",loop_tissue[,i]))# the contrast should be set in multigroup comparison
  res = res[order(res$pvalue),]
  summary(res)
  
  table(res$padj <0.05) # the distribution of the p-value
  res_data <- as.data.frame(res)
  res_data <- na.omit(res_data)
  write.csv(res_data,file = paste0("E15_ctx","_",colnames(loop_tissue[i]),"_all_results",".csv"))
  
  ### vocanal plot of the DEG
  logFC_cutoff <- with(res_data, mean(abs(log2FoldChange)) + 2* sd(log2FoldChange))# logFC_cutoff
  res_data$change <- as.factor(ifelse(res_data$padj <0.05 &   abs(res_data$log2FoldChange) > 1,
                                      ifelse(res_data$log2FoldChange > 1,'up','down'),'normal'))
  pig_title <- paste0("E15_ctx"," ",colnames(loop_tissue[i]),"(up:",
                      nrow(res_data[res_data$change == 'up',])," ","down:", 
                      nrow(res_data[res_data$change == 'down',]),")")
  g = ggplot(data = res_data,
             aes(x= log2FoldChange,y= -log10(padj), color = change)) + 
    geom_point(alpha =0.4 ,size =1.75) +
    theme_set(theme_set(theme_bw(base_size = 20))) + 
    xlab("log2 fold change") +ylab("-log10 padj") +
    ggtitle(pig_title)+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black') +
    geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
    scale_color_manual(values = c('blue','black','red')) ##corresponding to the levels(res_data$change)
  print(g)
  ggsave(g,filename = paste0("E15_ctx","_",colnames(loop_tissue[i]),"_deg_volcano.png"))
  
  ###extract the DEGs for the further analysis
  diff_gene <- subset(res_data,padj <0.05 & abs (log2FoldChange) >1)
  dim(diff_gene)
  write.csv(diff_gene, file = paste0("E15_ctx","_",colnames(loop_tissue[i]),"_DEG_results",".csv"))
  
  ## heatmap of the DEG
  choose_gene = rownames(diff_gene)
  choose_matrix <-normalized_counts_log[choose_gene,loop_deg[,i]] 
  #choose_matrix <-t(scale(t(choose_matrix)))
  choose_matrix <- na.omit(choose_matrix)
  pheatmap(choose_matrix,
           cluster_rows = T, 
           show_rownames = FALSE,
           cluster_cols = FALSE,
           #fontsize_col = 12,
           cellwidth = 40,
           filename = paste0("E15_ctx","_",colnames(loop_deg[i]),"_deg_heatmap.png")
           #annotation_col =df
  )
  ## get annotation of the deg genes(convert the ensemble_id to gene_symbol)
  
  ## merge the mms_symbols and the diff_gene
  diff_gene_id <- cbind(rownames(diff_gene),diff_gene)
  colnames(diff_gene_id)[1] <- c("ensembl_gene_id")
  mms_symbols_deg <-mms_symbols[choose_gene,]
  diff_gene_ano <- merge(diff_gene_id,mms_symbols,by = "ensembl_gene_id")
  write.csv(diff_gene_ano, file = paste0("E15_ctx","_",colnames(loop_tissue[i]),"_DEG_all_anno_results",".csv"))
  
  ### construct the multi_group DEG matrix for the plot
  diff_gene_ano$group<- i
  diff_gene_ano$label <- ifelse(diff_gene_ano$padj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
  top10sig0 <- filter(diff_gene_ano) %>% distinct(external_gene_name,.keep_all = T) %>% top_n(10,abs(log2FoldChange))
  diff_gene_ano$size <- case_when(!(diff_gene_ano$external_gene_name %in% top10sig0$external_gene_name)~ 1,
                                  diff_gene_ano$external_gene_name %in% top10sig0$external_gene_name ~ 2)
  top10sig <-rbind(top10sig,top10sig0)
  diff_multi<-rbind(diff_multi,diff_gene_ano)
  
  ### construct the list for the extraction of the intersection of the deg
  #gene_list[[i]] <- diff_gene_ano$ensembl_gene_id
  gene_list[[i]] <-diff_gene_ano %>% filter(change == "up") %$% .[,1]
  gene_list[[i+3]] <-diff_gene_ano %>% filter(change == "down") %$% .[,1]
}

diff_multi_no_10 <- filter(diff_multi,size == 1)
write.csv(diff_multi, file = "E15_ctx_DEG_group_group_results.csv")

### multigroup deg plot


#根据图p中log2FC区间确定背景柱长度：
dfbar<-data.frame(x=c(1,2,3), ##P1,P2,P3更改
                  y=c(13,10,12))
dfbar1<-data.frame(x=c(1,2,3),##P1,P2,P3更改
                   y=c(-14,-5,-15))


#添加X轴的group色块和标签：
label_up <-c()
label_down <-c()
for (a in 1:3) {##P1,P2,P3更改
  label_up<-c(label_up,dim(filter(diff_multi,group == a, change == "up"))[1])
  label_down<-c(label_down,dim(filter(diff_multi,group == a, change == "down"))[1])
}
dfcol_up<-data.frame(x=c(1:3),y=0,label=label_up)##P1,P2,P3更改
dfcol_down<-data.frame(x=c(1:3),y=0,label=label_down)##P1,P2,P3更改
dfcol_name <-data.frame(label=c("p3_vs_p0","p6_vs_p3","p6_vs_p0"),##P1,P2,P3更改
                        y=0,x=c(1:3))



dfcol<-data.frame(x=c(1:3),y=0,label=c(1:3))
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F")##P1,P2,P3更改
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+  #geom_col the background
  geom_jitter(data = diff_multi_no_10,
              aes(x = group, y = log2FoldChange, color = label),
              size = 0.5,
              width =0.4)+
  geom_jitter(data = top10sig,
              aes(x = group, y = log2FoldChange, color = label),
              size = 0.75,
              width =0.4)+ #geom_jitter the data of deg
  geom_tile(data = dfcol,
            aes(x=x,y=y),
            height=1.8,
            color = "black",
            fill = mycol,
            alpha = 0.6,
            show.legend = F)+ # add the color block
  geom_text_repel(data=top10sig,
                  aes(x=group,y=log2FoldChange,label=external_gene_name),
                  size = 2,
                  force = 0.2,
                  point.padding = NA,
                  segment.size = NA,
                  #arrow = arrow(length = unit(0.008, "npc"),
                  #type = "open", ends = "last"),
                  max.overlaps = 20) + # geom_text_repel add the top10 gene name
  labs(x="group",y="log2FoldChange")+
  geom_text(data=dfcol_name,
            aes(x=x,y=y,label=label),
            size =3,
            color ="black")+
  geom_text(data=dfcol_up,
            aes(x=x,y=y,label=label),
            nudge_y = 0.65,
            size =3,
            color ="white")+
  geom_text(data=dfcol_down,
            aes(x=x,y=y,label=label),
            nudge_y = -0.65,
            size =3,
            color ="white")+ # add the x/y lab and the group number
  theme_minimal()+     # 主题美化
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 10)
  )

print(p1)

ggsave(p1,filename = "E15_ctx_multi_deg_volcano.png")


########################### vennplot extract the intersection of the DEGs####################
#ctx_ge<- read.table("E:/yaner xuenian/NSCPRO/downstream/DEG/E15_P1/each_group_results/E15_P1_ctx_vs_ge_DEG_results.csv", sep=",", header = T)
#up_ctx_ge <- ctx_ge %>% filter(change == "up") %$% .[,1]


gene_list_venn <- list(list(gene_list[[1]],gene_list[[2]],gene_list[[3]]),      ##P1,P2,P3更改
                       list(gene_list[[1]],gene_list[[5]],gene_list[[3]]),
                       list(gene_list[[4]],gene_list[[5]],gene_list[[6]]),
                       list(gene_list[[4]],gene_list[[2]],gene_list[[6]]))

gene_ggvenn_name<-data.frame (up_up=c('p3-p0_up','p6-p3_up','p6-p0_up'),
                              up_down=c('p3-p0_up','p6-p3_down','p6-p0_up'),  ##P1,P2,P3更改
                              down_down= c('p3-p0_down','p6-p3_down','p6-p0_down'), 
                              down_up=c('p3-p0_down','p6-p3_up','p6-p0_down'))
gene_list_all <-list()
for (c in 1:4) {
  if(c%%2 == 0){
    p3<- ggVennDiagram(gene_list_venn[[c]],
                       category.names =gene_ggvenn_name[,c],
                       set_size = 8,
                       label_color = "black",
                       label_size = 8
    ) + scale_fill_gradient(low = "#43C6AC", high = "#F8FFAE")+
      scale_color_brewer(palette = "Greens")+
      scale_x_continuous(expand = expansion(mult = .2))
    print(p3 )
    ggsave(p3,filename = paste0("E15_ctx_",colnames(gene_ggvenn_name[c]),"_intersection_vennplot.png"))
    inter<- get.venn.partitions(gene_list_venn[[c]])
    for (d in 1:nrow(inter)) {
      inter[d,'values'] <- paste(inter[[d,'..values..']], collapse = ', ')
    }
    write.csv(inter[,-5],file = paste0("E15_ctx_",colnames(gene_ggvenn_name[c]),"_intersection.csv")) 
  }else{
    p2<- ggVennDiagram(gene_list_venn[[c]],
                       category.names = gene_ggvenn_name[,c],
                       set_size = 8,
                       label_color = "black",
                       label_size = 8
    ) + scale_fill_gradient2(low = "#c6ffdd",mid ="#fbd786" , high = "#f7797d")+
      scale_color_brewer(palette = "YlOrRd")+
      scale_x_continuous(expand = expansion(mult = .2))
    print(p2 )
    ggsave(p2,filename = paste0("E15_ctx_",colnames(gene_ggvenn_name[c]),"_intersection_vennplot.png"))
    inter<- get.venn.partitions(gene_list_venn[[c]])
    for (d in 1:nrow(inter)) {
      inter[d,'values'] <- paste(inter[[d,'..values..']], collapse = ', ')
    }
    write.csv(inter[,-5],file = paste0("E15_ctx_",colnames(gene_ggvenn_name[c]),"_intersection.csv"))
    
    
  }
  inter_gene <- data.frame(unlist(strsplit(inter[1,7], ', ')))
  colnames(inter_gene) <- c(inter[1,4])
  gene_list_all[c] <-inter_gene[1]
  # inter<- get.venn.partitions(gene_list_venn[[c]])
}

###
#inter_all <- do.call(cbind, lapply(lapply(inter$..values..,unlist ),`length<-`,max(lengths(inter$..values..))))# the all results of intersetcion



############################ GO and KEGG and GSEA enrichment analysis##############################