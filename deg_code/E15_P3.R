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
###################### load the data and construct the count matrix and design matrix#########
# input the rawcount files to build the matrix
p0_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0A1_L3_307X07_sorted.count", sep="\t", col.names = c("gene_id","p0_ctx1"))
p0_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0A2_L3_308X08_sorted.count", sep="\t", col.names = c("gene_id","p0_ctx2"))
p0_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0B1_L3_309X09_sorted.count", sep="\t", col.names = c("gene_id","p0_ge1"))
p0_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0B2_L4_310X10_sorted.count", sep="\t", col.names = c("gene_id","p0_ge2"))
p0_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0C1_L4_310X10_sorted.count", sep="\t", col.names = c("gene_id","p0_mid2"))
p0_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0C2_L4_311X11_sorted.count", sep="\t", col.names = c("gene_id","p0_mid1"))
p0_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0D1_L4_312X12_sorted.count", sep="\t", col.names = c("gene_id","p0_hid1"))
p0_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/E15-5-P0D2_L4_313X13_sorted.count", sep="\t", col.names = c("gene_id","p0_hid2"))
p3_ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX_L4_150A50_sorted.count", sep="\t", col.names = c("gene_id","p3_ctx3"))
p3_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX2_L4_346X46_sorted.count", sep="\t", col.names = c("gene_id","p3_ctx2"))
p3_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-CTX1_L3_343X43_sorted.count", sep="\t", col.names = c("gene_id","p3_ctx1"))
p3_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE1_L3_344X44_sorted.count", sep="\t", col.names = c("gene_id","p3_ge1"))
p3_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE2_L4_347X47_sorted.count", sep="\t", col.names = c("gene_id","p3_ge2"))
p3_ge3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-GE_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","p3_ge3"))
p3_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID1_L4_351X51_sorted.count", sep="\t", col.names = c("gene_id","p3_mid1"))
p3_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID2_L4_352X52_sorted.count", sep="\t", col.names = c("gene_id","p3_mid2"))
p3_mid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-MID_L4_152A52_sorted.count", sep="\t", col.names = c("gene_id","p3_mid3"))
p3_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID1_L3_345X45_sorted.count", sep="\t", col.names = c("gene_id","p3_hid1"))
p3_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID2_L3_349X49_sorted.count", sep="\t", col.names = c("gene_id","p3_hid2"))
p3_hid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P3-E15-5-HID_L4_153A53_sorted.count", sep="\t", col.names = c("gene_id","p3_hid3"))
p6_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-CTX1_L4_346X46_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx1"))
p6_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-CTX2_L4_349X49_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx2"))
p6_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-GE1_L4_347X47_sorted.count", sep="\t", col.names = c("gene_id","p6_ge1"))
p6_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-GE2_L4_350X50_sorted.count", sep="\t", col.names = c("gene_id","p6_ge2"))
p6_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-HID1_sorted.count", sep="\t", col.names = c("gene_id","p6_hid1"))
p6_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-hid_L4_379X79_sorted.count", sep="\t", col.names = c("gene_id","p6_hid2"))
p6_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-mid_L4_378X78_sorted.count", sep="\t", col.names = c("gene_id","p6_mid1"))
p6_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5/P6-E15-5-mid2.count", sep="\t", col.names = c("gene_id","p6_mid2"))

p0_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A1_L4_311X11_sorted.count", sep="\t", col.names = c("gene_id","p0_ctx1"))
p0_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0A2_L4_312X12_sorted.count", sep="\t", col.names = c("gene_id","p0_ctx2"))
p3_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX1_L3_350X50_sorted.count", sep="\t", col.names = c("gene_id","p3_ctx1"))
p3_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P3-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","p3_ctx2"))
p6_ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX1_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx1"))
p6_ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX2_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx2"))
p6_ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/P6-E13-5-CTX_sorted.count", sep="\t", col.names = c("gene_id","p6_ctx3"))
p0_ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B1_L4_313X13_sorted.count", sep="\t", col.names = c("gene_id","p0_ge1"))
p0_ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0B2_L4_380X80_sorted.count", sep="\t", col.names = c("gene_id","p0_ge2"))
p0_mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C2_L4_382X82_sorted.count", sep="\t", col.names = c("gene_id","p0_mid2"))
p0_mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0C1_L4_381X81_sorted.count", sep="\t", col.names = c("gene_id","p0_mid1"))
p0_hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D1_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","p0_hid1"))
p0_hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E13_5/E13-5-P0D2_L4_316X16_sorted.count", sep="\t", col.names = c("gene_id","p0_hid2"))
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
rawcount <-merge(merge(merge(ctx,ge,by ="gene_id"),mid,by ="gene_id"),hid,by ="gene_id")
rawcount_1 <- rawcount[-c(1:5),]
ENSEMBL <- gsub("\\.\\d*", "",rawcount_1$gene_id)
row.names(rawcount_1) <-ENSEMBL
rawcount_2 <- rawcount_1[,-1]
rawcount_3 <- rawcount_2[apply(rawcount_2, 1, sum) > 10 , ] #remove the low express gene
rawcount_log <-log2(rawcount_3+1)

## get annotation of the all genes(convert the ensemble_id to gene_symbol) for the deg loop annotation
mart <-useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <-row.names(rawcount_2)
mms_symbols <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","description","band"),
                     filters = "ensembl_gene_id", values = my_ensembl_gene_id, mart = mart)
rownames(mms_symbols)<- mms_symbols$ensembl_gene_id

#construct the coldata,3 bio replicates each group
condition <- factor(c(rep(c(rep("p0",2),rep("p3",2),rep("p6",3)),4)),levels = c("p0","p3","p6"))
tissue <- factor(rep(c("ctx","ge","mid","hid"),each = 7),levels = c("ctx","ge","mid","hid"))
batch <- factor(c(rep(c(rep("1",4),2,rep("1",2)),4)))# batch effect
colData <- data.frame(row.names = colnames(rawcount_3),condition,tissue)

###################### construct the loop for the analysis  #############################################
E13_P0 =c(1,2,8,9,15,16,22,23)
E13_P3 =E13_P0+2
E13_P6 =c(5,6,12,13,20,21,26,27)
count_loop<- data.frame( E13_P0,E13_P3,E13_P6)
rawcount_cl <- rawcount_3[,E13_P3]
rawcount_log_cl <-log2(rawcount_cl+1)
colData_cl <- colData[E13_P3,c(2,2)]

#remove the batch2 version

#start the DEG workflow
##construct the DESeqDataSet(dds) Object
dds <- DESeqDataSetFromMatrix(rawcount_cl,colData_cl,design = ~ tissue)
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
pca_plot<-plotPCA(rld_plot,intgroup = ("tissue"),returnData = T)
plotPCA(rld_plot,intgroup = ("tissue"))#package function
png("E13_P0_PCA.png")


E15_P6_PCA<- plotPCA(rld_plot,intgroup = ("tissue"))
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
                           filename = "E13_P0_dis.png",
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
         filename = "E13_P0_count_heat.png"
         #annotation_col =df
         )


###################### the DEG analysis and results extraction#########################
loop_tissue<-data.frame(ctx_vs_ge=c("ctx","ge"), ##P1,P2,P3更改
                        ctx_vs_mid=c("ctx","mid"), 
                        ctx_vs_hid=c("ctx","hid"),
                        ge_vs_mid=c("ge","mid"),
                        ge_vs_hid=c("ge","hid"),
                        mid_vs_hid=c("mid","hid"))
loop_deg <- data.frame(ctx_vs_ge=c(1:4), ##P1,P2,P3更改
                       ctx_vs_mid=c(1,2,5,6), 
                       ctx_vs_hid=c(1,2,7,8),
                       ge_vs_mid=c(3:6),
                       ge_vs_hid=c(3,4,7,8),
                       mid_vs_hid=c(5:8))
top10sig <-c() #for the multi_deg plot
diff_multi<-c()#for the multi_deg plot
gene_list <-list() # for the vennplot
#gene_list <-list()#for the extraction of the intersection of the multi_group deg
for (i in 1:6) {##P1,P2,P3更改
res <-results(dds,contrast = c("tissue",loop_tissue[,i]))# the contrast should be set in multigroup comparison
res = res[order(res$pvalue),]
summary(res)

table(res$padj <0.05) # the distribution of the p-value
res_data <- as.data.frame(res)
res_data <- na.omit(res_data)
write.csv(res_data,file = paste0("E13_P3","_",colnames(loop_tissue[i]),"_all_results",".csv"))

### vocanal plot of the DEG
logFC_cutoff <- with(res_data, mean(abs(log2FoldChange)) + 2* sd(log2FoldChange))# logFC_cutoff
res_data$change <- as.factor(ifelse(res_data$padj <0.05 &   abs(res_data$log2FoldChange) > 1,
                                    ifelse(res_data$log2FoldChange > 1,'up','down'),'normal'))
pig_title <- paste0("E13_P3"," ",colnames(loop_tissue[i]),"(up:",
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
ggsave(g,filename = paste0("E13_P3","_",colnames(loop_tissue[i]),"_deg_volcano.png"))

###extract the DEGs for the further analysis
diff_gene <- subset(res_data,padj <0.05 & abs (log2FoldChange) >1)
dim(diff_gene)
write.csv(diff_gene, file = paste0("E13_P3","_",colnames(loop_tissue[i]),"_DEG_results",".csv"))

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
         filename = paste0("E13_P3","_",colnames(loop_deg[i]),"_deg_heatmap.png")
         #annotation_col =df
         )
## get annotation of the deg genes(convert the ensemble_id to gene_symbol)

## merge the mms_symbols and the diff_gene
diff_gene_id <- cbind(rownames(diff_gene),diff_gene)
colnames(diff_gene_id)[1] <- c("ensembl_gene_id")
mms_symbols_deg <-mms_symbols[choose_gene,]
diff_gene_ano <- merge(diff_gene_id,mms_symbols,by = "ensembl_gene_id")
write.csv(diff_gene_ano, file = paste0("E13_P3","_",colnames(loop_tissue[i]),"_DEG_all_anno_results",".csv"))

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
gene_list[[i+6]] <-diff_gene_ano %>% filter(change == "down") %$% .[,1]
}

diff_multi_no_10 <- filter(diff_multi,size == 1)
write.csv(diff_multi, file = "E13_P3_DEG_group_group_results.csv")

### multigroup deg plot


#根据图p中log2FC区间确定背景柱长度：
dfbar<-data.frame(x=c(1,2,3,4,5,6), ##P1,P2,P3更改
                  y=c(5,9,11,11,12,11))
dfbar1<-data.frame(x=c(1,2,3,4,5,6),##P1,P2,P3更改
                   y=c(-11,-11,-14,-10,-15,-4))


#添加X轴的group色块和标签：
label_up <-c()
label_down <-c()
for (a in 1:6) {##P1,P2,P3更改
  label_up<-c(label_up,dim(filter(diff_multi,group == a, change == "up"))[1])
  label_down<-c(label_down,dim(filter(diff_multi,group == a, change == "down"))[1])
}
dfcol_up<-data.frame(x=c(1:6),y=0,label=label_up)##P1,P2,P3更改
dfcol_down<-data.frame(x=c(1:6),y=0,label=label_down)##P1,P2,P3更改
dfcol_name <-data.frame(label=c("ctx-ge","ctx-mid","ctx-hid","ge-mid","ge-hid","mid-hid"),##P1,P2,P3更改
                        y=0,x=c(1:6))

dfcol<-data.frame(x=c(1:6),y=0,label=c(1:6))
mycol <- c("#E64B357F","#4DBBD57F","#00A0877F","#F39B7F7F","#8491B47F","#91D1C27F")##P1,P2,P3更改
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

ggsave(p1,filename = "E13_P0_multi_deg_volcano.png")


########################### vennplot extract the intersection of the DEGs####################
#ctx_ge<- read.table("E:/yaner xuenian/NSCPRO/downstream/DEG/E15_P1/each_group_results/E15_P1_ctx_vs_ge_DEG_results.csv", sep=",", header = T)
#up_ctx_ge <- ctx_ge %>% filter(change == "up") %$% .[,1]


gene_list_venn <- list(list(gene_list[[1]],gene_list[[2]],gene_list[[3]]),      ##P1,P2,P3更改
                       list(gene_list[[7]],gene_list[[8]],gene_list[[9]]),
                       list(gene_list[[4]],gene_list[[5]],gene_list[[7]]),
                       list(gene_list[[10]],gene_list[[11]],gene_list[[1]]),
                       list(gene_list[[6]],gene_list[[8]],gene_list[[10]]),
                       list(gene_list[[12]],gene_list[[2]],gene_list[[4]]),
                       list(gene_list[[9]],gene_list[[11]],gene_list[[12]]),
                       list(gene_list[[3]],gene_list[[5]],gene_list[[6]])
                       )
gene_ggvenn_name<-data.frame (ctx=c('ge','mid','hid'),ctx=c('ge','mid','hid'),  ##P1,P2,P3更改
                              ge= c('mid','hid',"ctx"), ge=c('mid','hid',"ctx"),
                              mid=c('hid',"ctx","ge"), mid=c('hid',"ctx","ge"),
                              hid= c('ctx',"ge","mid"), hid=c('ctx',"ge","mid"))
gene_list_all <-list()
for (c in 1:8) {
  if(c%%2 == 0){
    p3<- ggVennDiagram(gene_list_venn[[c]],
                       category.names =gene_ggvenn_name[,c],
                       set_size = 10,
                       label_color = "black",
                       label_size = 8
    ) + scale_fill_gradient(low = "#43C6AC", high = "#F8FFAE")+
      scale_color_brewer(palette = "Greens")
    print(p3 )
    ggsave(p3,filename = paste0("E13_P3_down_",colnames(gene_ggvenn_name[c]),"_intersection_vennplot.png"))
    inter<- get.venn.partitions(gene_list_venn[[c]])
    for (d in 1:nrow(inter)) {
      inter[d,'values'] <- paste(inter[[d,'..values..']], collapse = ', ')
    }
    write.csv(inter[,-5],file = paste0("E13_P3_down_",colnames(gene_ggvenn_name[c]),"_intersection.csv")) 
  }else{
    p2<- ggVennDiagram(gene_list_venn[[c]],
                       category.names = gene_ggvenn_name[,c],
                       set_size = 10,
                       label_color = "black",
                       label_size = 8
    ) + scale_fill_gradient2(low = "#c6ffdd",mid ="#fbd786" , high = "#f7797d")+
      scale_color_brewer(palette = "YlOrRd")
    print(p2 )
    ggsave(p2,filename = paste0("E13_P3_up_",colnames(gene_ggvenn_name[c]),"_intersection_vennplot.png"))
    inter<- get.venn.partitions(gene_list_venn[[c]])
    for (d in 1:nrow(inter)) {
      inter[d,'values'] <- paste(inter[[d,'..values..']], collapse = ', ')
    }
    write.csv(inter[,-5],file = paste0("E13_P3_up_",colnames(gene_ggvenn_name[c]),"_intersection.csv"))
    

    }
  inter_gene <- data.frame(unlist(strsplit(inter[1,7], ', ')))
  colnames(inter_gene) <- c(inter[1,4])
  gene_list_all[c] <-inter_gene[1]
 # inter<- get.venn.partitions(gene_list_venn[[c]])
}

###
#inter_all <- do.call(cbind, lapply(lapply(inter$..values..,unlist ),`length<-`,max(lengths(inter$..values..))))# the all results of intersetcion



############################ GO and KEGG and GSEA enrichment analysis##############################
BiocManager::install('clusterProfiler')
BiocManager::install('DOSE')
BiocManager::install('org.Mm.eg.db')
BiocManager::install('stringr')
library(clusterProfiler)
library(DOSE)
library(DO.db)
library(org.Mm.eg.db)
library(stringr)

## construct the dataframe that contained the entrez id for the further analysis
gene_enrich_up<- list(gene_list_all[[1]],gene_list_all[[3]],gene_list_all[[5]],gene_list_all[[7]])

gene_list_entrezid <-list()
for (e in 1:8) {
  gene_df <-bitr(gene_list_all[[e]],fromType = 'ENSEMBL',
                 toType = c('SYMBOL','ENTREZID'),
                 OrgDb = org.Mm.eg.db)
  gene_list_entrezid[e] <-gene_df[3]
}
gene_enrich_up<- list(gene_list_entrezid[[1]],gene_list_entrezid[[3]],
                      gene_list_entrezid[[5]],gene_list_entrezid[[7]])
names(gene_enrich_up)<-c('ctx','ge','mid','hid')
gene_df <-bitr(gene_list_all[[1]],fromType = 'ENSEMBL',
               toType = c('SYMBOL','ENTREZID'),
               OrgDb = org.Mm.eg.db)
ego<- compareCluster(geneClusters = gene_enrich_up,fun = enrichKEGG)
ego <- setReadable(ego, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
table(ego@compareClusterResult$Cluster)
dotplot(ego_cc)

## go classification 

ego_cc <- enrichGO(gene = gene_df$ENSEMBL,
                   OrgDb = org.Mm.eg.db,
                   keyType = 'ENSEMBL',
                   ont = 'CC',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.2,
                   qvalueCutoff = 0.2)
ego_bp <- enrichGO(gene = gene_df_ABC$ENSEMBL,
                   OrgDb = org.Mm.eg.db,
                   keyType = 'ENSEMBL',
                   ont = 'BP',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)
ego_MF <- enrichGO(gene = gene_df_ABC$ENSEMBL,
                   OrgDb = org.Mm.eg.db,
                   keyType = 'ENSEMBL',
                   ont = 'MF',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)
head(ego_cc@result)
ego_cc@
barplot(ego_cc,showCategory = 10) +title = ('the GO_BP enrichment analysis of E15_P3_hid')+
  scale_size(range = c(5,6))+
  scale_x_discrete(labels =function(ego_bp) strwrap(ego_bp,width = 50))
barplot(ego_cc,showCategory = 10,title = 'the GO_CC enrichment analysis of E15_P3_hid')
barplot(ego_MF,showCategory = 10,title = 'the GO_MF enrichment analysis of E15_P3_hid')

## kegg enrichment

options(clusterProfiler.download.method = "wininet")
kk<- enrichKEGG(gene = gene_df_ABC$ENTREZID,
                organism = 'mmu',
                pvalueCutoff = 0.05)
dim(kk)
kk[1:30]
barplot(kk, showCategory = 25,title = 'The KEGG enrichment analysis of E15_P3_hid')+
  scale_size(range = c(2,12))+
  scale_x_discrete(labels = function(kk) str_wrap(kk,width = 25))
dotplot(kk,showCategory = 25,title = 'The KEGG enrichment analysis of E15_P3_hid')+
  scale_size(range = c(2,12))+
  scale_x_discrete(labels = function(kk) str_wrap(kk,width = 25))

## Gene Set Enrichment Analysis (GSEA)
###  获取按照log2FC 大小排序的基因列表
gsea_list <- ctx_vs_mid$log2FoldChange
names(gsea_list) <- ctx_vs_mid[,2]
gsea_list<- sort(gsea_list,decreasing = TRUE)
### GSEA analysis
gsemf <- gseGO(gsea_list,
               OrgDb = org.Mm.eg.db,
               keyType = 'ENSEMBL',
               ont = 'BP',
               pvalueCutoff = 0.5)
gseaplot(gsemf,geneSetID = 'GO:0001819')
head(gsemf)
