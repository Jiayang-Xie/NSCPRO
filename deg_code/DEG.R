# install and load the relevant packages
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
library(BiocManager)
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("biomaRt")# gene id convert
BiocManager::install("curl")
BiocManager::install("limma")#removeBatchEffect
BiocManager::install("sva")#removeBatchEffect
BiocManager::install("GEOquery")
library(DESeq2)
library(biomaRt)
library(curl)
library(pheatmap)
library(GEOquery)
library(ggplot2)
# input the rawcount files to build the matrix
ctx3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-CTX_L4_150A50_sorted.count", sep="\t", col.names = c("gene_id","ctx3"))
ctx2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-CTX2_L4_346X46_sorted.count", sep="\t", col.names = c("gene_id","ctx2"))
ctx1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-CTX1_L3_343X43_sorted.count", sep="\t", col.names = c("gene_id","ctx1"))
ge1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-GE1_L3_344X44_sorted.count", sep="\t", col.names = c("gene_id","ge1"))
ge2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-GE2_L4_347X47_sorted.count", sep="\t", col.names = c("gene_id","ge2"))
ge3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-GE_L4_314X14_sorted.count", sep="\t", col.names = c("gene_id","ge3"))
mid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-MID1_L4_351X51_sorted.count", sep="\t", col.names = c("gene_id","mid1"))
mid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-MID2_L4_352X52_sorted.count", sep="\t", col.names = c("gene_id","mid2"))
mid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-MID_L4_152A52_sorted.count", sep="\t", col.names = c("gene_id","mid3"))
hid1 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-HID1_L3_345X45_sorted.count", sep="\t", col.names = c("gene_id","hid1"))
hid2 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-HID2_L3_349X49_sorted.count", sep="\t", col.names = c("gene_id","hid2"))
hid3 <- read.table("E:/yaner xuenian/NSCPRO/downstream/E15_5-P3/P3-E15-5-HID_L4_153A53_sorted.count", sep="\t", col.names = c("gene_id","hid3"))
ctx <- merge(merge(ctx1,ctx2,by ="gene_id"),ctx3,by ="gene_id")
ge <- merge(merge(ge1,ge2,by ="gene_id"),ge3,by ="gene_id")
mid <- merge(merge(mid1,mid2,by ="gene_id"),mid3,by ="gene_id")
hid <- merge(merge(hid1,hid2,by ="gene_id"),hid3,by ="gene_id")
rawcount<-merge(merge(merge(ctx,ge,by ="gene_id"),mid,by ="gene_id"), hid,by ="gene_id")
rawcount_1 <- rawcount[-c(1:5),]
ENSEMBL <- gsub("\\.\\d*", "",rawcount_1$gene_id)
row.names(rawcount_1) <-ENSEMBL
rawcount_2 <- rawcount_1[,-1]
rawcount_3 <- rawcount_2[apply(rawcount_2, 1, sum) > 1 , ] #remove the low express gene
rawcount_log <-log2(rawcount_3+1)

#construct the coldata,3 bio replicates each group
condition <- factor(c(rep("ctx",3),rep("ge",3),rep("mid",3),rep("hid",3)),levels = c("ctx","ge","mid","hid"))
batch <- factor(c(rep("bat1",2),rep("bat2",1),rep("bat1",2),rep("bat2",1),rep("bat1",2),rep("bat2",1),
                  rep("bat1",2),rep("bat2",1)))# batch effect
colData <- data.frame(row.names = colnames(rawcount_3),condition,batch)

#start the DEG workflow
##construct the DESeqDataSet(dds) Object

dds <- DESeqDataSetFromMatrix(rawcount_3,colData,design = ~ batch + condition)
dds<- DESeq(dds)
dds

## the visualization of data(evaluate the samples before extract the results)
### Dispersion plot
plotDispEsts(dds, main="Dispersion plot")

###PCA (principal components analysis),before the PCA;the count matrix should be feature normalized by rlog or vst 
#### 归一化 blind，转换时是否忽视实验设计。blind=T，不考虑实验设计，用于样品质量保证（sample quality assurance，QA）
#### blind=F，考虑实验设计，用于downstream analysis。
vstdata <-vst(dds,blind = FALSE) 
assay(vstdata) <- limma::removeBatchEffect(assay(vstdata),vstdata$batch) #去批次效应
plotPCA(vstdata,intgroup = ("condition"))#package function
vst_plot<-plotPCA(vstdata,intgroup = ("condition"),returnData =T)
plot(vst_plot[,1:2],pch=19,col= c(rep("red",3),rep("blue",3),rep("green",3),rep("yellow",3)),cex=1.5,
     main = "The PCA analysis of P3-E15.5",xlab ="PC1: 44% variance",ylab="PC2: 22% variance",
     col.main ="black",xlim=c(-15,15))
text(vst_plot[,1],vst_plot[,2],row.names(vst_plot),cex=1, font = 1,adj = 1)
legend("topright",inset=.05,title="group",c("ctx","ge","mid","hid"),col = c("red","blue","green","yellow"),
       pch = 19,cex = 0.75,bty = "n")
###PCA by rlog
rld_plot <- rlogTransformation(dds)
plotPCA(rld_plot,intgroup = ("condition"))

### lfshrink and MAplot
library(apeglm)  
resultsNames(dds_norm)  #看一下要shrink的维度;shrink数据更加紧凑,少了一项stat，但并未改变padj，但改变了foldchange
res_shrink <- lfcShrink(dds_norm, coef="condition_treat2_vs_treat1", type="apeglm") #最推荐apeglm算法;根据resultsNames(dds)的第5个维度，coef=5，也可直接""指定;apeglm不allow contrast，所以要指定coef
pdf("MAplot.pdf", width = 6, height = 6) 
plotMA(res_shrink, ylim=c(-10,10), alpha=0.1, main="MA plot: ")
dev.off()


### sample to sample distance heatmap
library(RColorBrewer)
sampleDist_vstdata <-dist(t(assay(vstdata))) #vst transformation
sampleDist_vstdata_matrix <- as.matrix(sampleDist_vstdata)
rownames(sampleDist_vstdata_matrix)<-paste(vstdata$condition,sep = "-")
colnames(sampleDist_vstdata_matrix)<- NULL
colors<- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDist_vstdata_matrix,
         clustering_distance_rows = sampleDist_vstdata,
         clustering_distance_cols = sampleDist_vstdata,
         col =colors)

### count matrix 
slect <-order(rowMeans(counts(dds,normalized =TRUE)),decreasing = TRUE)[1:1000] #按照每行counts平均数排序，取指定数量
df <- as.data.frame(colData(dds)[,c("condition","sizeFactor","batch")])#图例注释信息
ntd <- normTransform(dds) # log2(n+1)
pheatmap(assay(ntd)[slect,],cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col =df)

## the overview of the result
res <-results(dds,contrast = c("condition","ge","hid"))# the contrast should be set in multigroup comparison
res = res[order(res$pvalue),]
head(res)
summary(res)

table(res$padj <0.05) # the distribution of the p-value
res_data <- as.data.frame(res)
res_data <- na.omit(res_data)
write.csv(res_data,file = "all_results_ge_vs_hid")

## vocanal plot of the DEG
logFC_cutoff <- with(res_data, mean(abs(log2FoldChange)) + 2* sd(log2FoldChange))# logFC_cutoff
res_data$change <- as.factor(ifelse(res_data$pvalue <0.05 &   abs(res_data$log2FoldChange) > logFC_cutoff,
                                    ifelse(res_data$log2FoldChange > logFC_cutoff,'up','down'),'not'))
pig_title <- paste0('cut_off for logFC is ',round(logFC_cutoff,3),'\nThe number of up gene is',
                     nrow(res_data[res_data$change == 'up',]), '\nThe number of down gene is', 
                     nrow(res_data[res_data$change == 'down',]))
g = ggplot(data = res_data,
           aes(x= log2FoldChange,y= -log10(pvalue), color = change)) + 
           geom_point(alpha =0.4 ,size =1.75) +
           theme_set(theme_set(theme_bw(base_size = 20))) + 
           xlab("log2 fold change") +ylab("-log10 pvalue") +
           ggtitle(pig_title) + 
           theme(plot.title = element_text(size = 15,hjust = 0.5)) +
           scale_color_manual(values = c('blue','black','red')) ##corresponding to the levels(res_data$change)
print(g)
ggsave(g,filename = 'volcano.png')

## extract the DEGs for the further analysis
diff_gene <- subset(res,padj <0.05 & abs (log2FoldChange) >1)
dim(diff_gene)
head(diff_gene,200)
summary(diff_gene)
diff_gene_data <- as.data.frame(diff_gene)
write.csv(diff_gene, file = "E15_P3_DEG_ge_vs_hid.csv")

## heatmap of the DEG
choose_gene = head(rownames(diff_gene),100)
choose_matrix <-rawcount_log[choose_gene,] 
choose_matrix <-t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename = 'E15_P3_DEG_ge_vs_hid.png')



## get annotation of the genes(convert the ensemble_id to gene_symbol)
mart <-useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <-row.names(diff_gene_data)
mms_symbols <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","description","band"),
                     filters = "ensembl_gene_id", values = my_ensembl_gene_id, mart = mart)
## merge the mms_symbols and the diff_gene
diff_gene_id <- cbind(rownames(diff_gene_data),diff_gene_data)
colnames(diff_gene_id)[1] <- c("ensembl_gene_id")
diff_gene_ano <- merge(diff_gene_id,mms_symbols,by = "ensembl_gene_id")
write.csv(diff_gene_ano, file = "DEG_all_ge_vs_hid.csv")

# GO analysis
