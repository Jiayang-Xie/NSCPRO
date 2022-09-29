options(repos="https://mirror.lzu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("GenomicAlignments", type="binary")
BiocManager::install("tximportData")
BiocManager::install("stringr")
BiocManager::install("tidyverse")
library(BiocManager)
library(DESeq2)
library(biomaRt)
library(curl)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(installr)
library(tximportData)
library(tximport)
library(stringr)
library(tidyverse)

rawcount_1<- read.csv(file = "E:/Sequence/downstream analysis/DEG/all_samples_counts_gene_name.csv")
codldata <- read.csv(file = "E:/Sequence/downstream analysis/DEG/all_samples_traits.csv")
rawcount_1 <- rawcount_1[apply(rawcount_1, 1, sum) > 10 , ] 
test <- rawcount_1$E13_p0_ctx1
colSums(rawcount_1[,2:57])


E13_P0 =c(1,2,8,9,15,16,22,23)
E13_P3 =E13_P0+2
E13_P6 =c(5,6,12,13,20,21,26,27)
count_loop<- data.frame( E13_P0,E13_P3,E13_P6)
rawcount_cl <- rawcount_3[,E13_P3]
rawcount_log_cl <-log2(rawcount_cl+1)
colData_cl <- colData[E13_P3,c(2,2)]






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


##φ(*￣0￣)ヾ(≧▽≦*)o(～￣▽￣)～o(*^＠^*)o`(*>﹏<*)′🤣
txdb <-makeTxDbFromGFF(file = "E:/Sequence/downstream analysis/DEG/gencode.vM30.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
df <- select(txdb, keys = k,"GENEID",  "TXNAME")
write.csv(df, file = "tx2gene_vM30_anno.csv")



############ 1 hisat2 + featurecount(gencode.vM25.annotation.gtf)   ###########
E13_P0 <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/E15_P0.txt")
E13_P3 <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/E15_P3.txt")
E13_P6 <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/E15_P6.txt")
E15_P0 <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/E15_P0.txt")
E15_P3 <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/E15_P3.txt")
E15_P6 <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/E15_P6.txt")
coldata <- read.csv("E:/Sequence/downstream analysis/DEG/all_samples_traits.csv",stringsAsFactors = T,row.names = 1)
coldata_1 <-( coldata[-56,])
coldata_1$stage <-as.factor(coldata_1$stage) 
coldata_1$batch <-as.factor(coldata_1$batch)
coldata_1$condition <- as.factor(coldata_1$condition)
rownames(coldata_1) <-colnames(normalized_counts_t)

raw_count <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/all_raw_count.txt")
raw_count_1 <- raw_count[,-c(10,19,32,41,54)]
colnames(raw_count_1) <- c("geneid",
                           "E13_P0_CTX1","E13_P0_CTX2","E13_P0_GE1","E13_P0_GE2","E13_P0_MID1","E13_P0_MID2","E13_P0_HID1","E13_P0_HID2",
                           "E13_P3_CTX1","E13_P3_CTX2","E13_P3_GE1","E13_P3_GE2","E13_P3_HID1","E13_P3_HID2","E13_P3_MID1","E13_P3_MID2",
                           "E13_P6_CTX1","E13_P6_CTX2","E13_P6_CTX3","E13_P6_GE1","E13_P6_GE2","E13_P6_GE3",
                           "E13_P6_HID1","E13_P6_HID2","E13_P6_HID3","E13_P6_MID1","E13_P6_MID2","E13_P6_MID3",
                           "E15_P0_CTX1","E15_P0_CTX2","E15_P0_GE1","E15_P0_GE2","E15_P0_MID1","E15_P0_MID2","E15_P0_HID1","E15_P0_HID2",
                           "E15_P3_CTX1","E15_P3_CTX2","E15_P3_CTX3","E15_P3_GE1","E15_P3_GE2","E15_P3_GE3",
                           "E15_P3_HID1","E15_P3_HID2","E15_P3_HID3","E15_P3_MID1","E15_P3_MID2","E15_P3_MID3",
                           "E16_P6_CTX1","E15_P6_CTX2","E15_P6_GE1","E15_P6_GE2","E15_P6_HID1","E15_P6_HID2","E15_P6_MID1")
##### 继续计算
rawcount_2 <- raw_count_1[-1,]
ENSEMBL <- gsub("\\.\\d*", "",rawcount_2$geneid)
rownames(rawcount_2)<- ENSEMBL
rawcount_2 <- rawcount_2[,-1]

rawcount_3 <- as.data.frame(lapply(rawcount_2, as.integer))
colSums(rawcount_3)
rownames(rawcount_3)<-rownames(rawcount_2)

######################## 转换id(all)#####################
mart <-useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <-row.names(rawcount_2)
mms_symbols <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","description","band"),
                     filters = "ensembl_gene_id", values = ENSEMBL, mart = mart)
write.csv(mms_symbols,file ="mms_symbols_gene_name.csv")
rownames(mms_symbols)<- mms_symbols$ensembl_gene_id
rawcount_3$ensembl_gene_id <-rownames(rawcount_3)
rawcount_3_n1 <- merge(rawcount_3,mms_symbols[,1:2],by ="ensembl_gene_id")
rawcount_3_n1 <- rawcount_3_n1 %>% distinct(external_gene_name, .keep_all = TRUE)##因为有重复的gene name，要去掉
row.names(rawcount_3_n1) <- rawcount_3_n1$external_gene_name
rawcount_3_n1 <- rawcount_3_n1[,c(-1,-57)]   ### all samples with no filter 
write.csv(rawcount_3_n1,file ="all_samples_counts_gene_name.csv")





################ heatmap for all  ######################
rawcount_3 <- rawcount_3[,-56]
rawcount_3 <- rawcount_3[apply(rawcount_3[,-56],1,sum) > 20,] #remove the low express gene
rawcount_log <-log2(rawcount_3+1)

mad_counts_log <- apply(rawcount_log, 1, mad)
mad_counts_log <- rawcount_log[order(mad_counts_log, decreasing=T), ]#for count-heatmap
mad_counts_log_p0 <-mad_counts_log[,c(1:8,29:36)]
rawcount_log<- na.omit(rawcount_log)



pheatmap(mad_counts_log_p0[1:500,],
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize_col = 12,
         cellwidth = 40,
         annotation_col = coldata_1,
         filename = "P0_mad_500_1.png"
         #annotation_col =df
)


##################### 合并单个不同来源count文件 ###############
# E15_P6_mid2 <- read.table("E:/Sequence/downstream analysis/DEG/featurecounts/E15_P6_mid2.sf")
# E15_P6_mid2 <-E15_P6_mid2[-1,c(1,5)]
# colnames(E15_P6_mid2)<- c("geneid","E15_P6_MID2")
# ENSEMBL_mid <- gsub("\\.\\d*", "",E15_P6_mid2$geneid)
# 
# E15_P6_mid2$gene_id <- ENSEMBL_mid
# E15_P6_mid2 <- E15_P6_mid2[,-1]
# raw_count_1 <- raw_count_1[-1,]
# ENSEMBL <- gsub("\\.\\d*", "",raw_count_1$geneid)
# 
# rawcount_2 <- raw_count_1[,-1]
# rawcount_2$gene_id <- ENSEMBL
# rawcount_2<- merge(rawcount_2,E15_P6_mid2, by="gene_id")
# rawcount_2_1 <- rawcount_2[,-1]




################ the normalized exp by DEseq2 ####################


dds <- DESeqDataSetFromMatrix(rawcount_3,coldata_1,design = ~ stage+ tissue +condition )
dds<- DESeq(dds)

normalized_counts <-counts(dds,normalized=T)
normalized_counts<- na.omit(normalized_counts)
vstdata <-vst(dds,blind = T) 
vstdata_2 <- assay(vstdata)
vstdata_2 <- as.data.frame(vstdata_2)
write.csv(vstdata_2,file = "all_samples_by_vst.csv")
vstdata_1 <- as.data.frame(vstdata_1) 

write.csv(normalized_counts,file ="all_samples_counts_by_deseq2.csv")
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- normalized_counts[,c(1:12,15,16,13,14,17:22,26:28,23:25,29:42,46:48,43:45,49:52,55,53,54)]
normalized_counts<-read.csv("E:/Sequence/downstream analysis/DEG/all_samples_counts_by_deseq2.csv", header = T,row.names = 1)
normalized_counts_t<-read.csv("E:/Sequence/downstream analysis/DEG/all_samples_counts_by_deseq2.csv", header = T,row.names = 1)

############# the gene list of regional markers ######################


nsc <- c("SXO1","FABP7","DLL1","HES5","NOTCH1","NES","SOX2")
nsc <- str_to_title(nsc)
ctx <-c("Foxg1","Emx1","Emx2","Pax6","Lhx2","Six3","Tbr1","Tbr2","Gli3","Nr2f1","Otx2")
ctx <- str_to_title(ctx)
ge <- c("Dlx2","LHX6","GSX2","NKX2.1","Couptfii","GAD2","GAD1","Olig2","Isl1")
ge <- str_to_title(ge)
mid <- c("En1","pax2","otx2","En2","Sim1","Pax5","Lmx1a","Pax2","Lhx1","Lmx1b","Wnt1",)
mid <- str_to_title(mid)
hid <- c("Gbx2","FGF8","EN2","KROX20","IRX3","HOXB4","HOXB2","pax2","pax3","lhx9","pitx1","nrg1","pax8")
hid <- str_to_title(hid)
mid_hid <- c("PAX3","DMBX1","DBX1","EN2","HOXA2","HOXB2","HOXB3","PAX7","PHOX2B","WNT1","BARHL1","ATOH1")

mid_hid <- str_to_title(mid_hid)
neuron <- c("stmn2","gap43","gng3","tubb3","mapt","tmem130","gria2","Syp","dcx","l1cam","nrxn3","eno2","pcdha5",
            "elavl2","elavl3","map2")
neuron <- str_to_title(neuron)
N_progenitor <- c("Sox1","Nes","Vim","Sox2","Hes5","Hes1","Notch1","Notch2","gpr98","sparc","slc1a3",
                  "fabp7","bmp7","anxa2")
N_progenitor <- str_to_title(N_progenitor)
regional_marker <- c(ctx,ge,mid,hid,nsc,N_progenitor,neuron)
regional_marker_name <- list(ctx,ge,mid,hid,nsc,N_progenitor,neuron)


###### heatmap of markers in brain region  ###############
counts_region_all <- normalized_counts[which(rownames(normalized_counts)%in%regional_marker),]
counts_region_all_log <-log2(counts_region_all+1)

counts_ctx <- normalized_counts[which(rownames(normalized_counts)%in%ctx),]
counts_ctx_log <-log2(counts_ctx+1)
counts_ge <- normalized_counts[which(rownames(normalized_counts)%in%ge),]
counts_ge_log <-log2(counts_ge+1)
counts_mid <- normalized_counts[which(rownames(normalized_counts)%in%mid),]
counts_mid_log <-log2(counts_mid+1)
counts_hid <- normalized_counts[which(rownames(normalized_counts)%in%hid),]
counts_hid_log <-log2(counts_hid+1)
counts_nsc <- normalized_counts[which(rownames(normalized_counts)%in%nsc),]
counts_nsc_log <-log2(counts_nsc+1)
counts_neuron <- normalized_counts[which(rownames(normalized_counts)%in%neuron),]
counts_progenitor <- normalized_counts[which(rownames(normalized_counts)%in%N_progenitor),]
counts_progenitor_log <-log2(counts_progenitor+1)

counts_E_name <- c("counts_ctx","counts_ge","counts_mid","counts_hid","counts_nsc","counts_progenitor")
counts_E_log_name <- c("counts_ctx_log","counts_ge_log","counts_mid_log","counts_hid_log",
                       "counts_nsc_log","counts_progenitor_log")

E13_P0 =c(1:8)
E13_P3 =c(9:16)
E13_P6 =c(17:28)
E15_P0 =c(29:36)
E15_P3 =c(37:48)
E15_P6 =c(49:55)
count_loop<- list(E13_P0,E13_P3,E13_P6,E15_P0,E15_P3,E15_P6)
count_loop_name <- c("E13_P0","E13_P3","E13_P6","E15_P0","E15_P3","E15_P6")
regional_name <- c("ctx","ge","mid","hid","nsc","N_progenitor","neuron")
bk <- c(seq(0,7.9,by =0.01),seq(8,16,by=0.01))

for (a in 1:7) {
  counts_region <- normalized_counts[which(rownames(normalized_counts)%in%regional_marker_name[[a]]),] 
  counts_region_log <- log2(counts_region+1)
  # for (b in 1:6) {
  #   counts_E_P <- counts_region[,count_loop[[b]]]
  # 
  #   pheatmap(counts_E_P,
  #            cluster_rows = T,
  #            cluster_cols = F,
  #            show_rownames = T,
  #            show_colnames = T,
  #            fontsize_col = 12,
  #            cellwidth = 40,
  #            annotation_col = coldata_1,
  #            filename = paste0("E:/Sequence/downstream analysis/DEG/counts_genes/",
  #                              regional_name[a],"_", count_loop_name[b],"counts_heatmap.png")
  #   )
  # }
  # 
  for (b in 1:6) {
    counts_E_P_log <- counts_region_log[,count_loop[[b]]]
    
    pheatmap(counts_E_P_log,
             scale = "none",
             cluster_rows = F,
             cluster_cols = F,
             show_rownames = T,
             show_colnames = T,
             fontsize_col = 12,
             cellwidth = 40,
             color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                       colorRampPalette(colors = c("white","red"))(length(bk)/2)),
             breaks = bk,
             #annotation_col = coldata_1,
             filename = paste0("E:/Sequence/downstream analysis/DEG/count_heatmap/",
                               regional_name[a],"_", count_loop_name[b],"_counts_log_heatmap.png")
    ) 
  
  }
}

###### line_box of markers in passage ############

E13_ctx =c(1,2,9,10,17,18,19)
E13_ge =c(3,4,11,12,20,21,22)
E13_mid =c(5,6,13,14,23,24,25)
E13_hid =c(7,8,15,16,26,27,28)
E15_ctx =c(29,30,37,38,39,49,50)
E15_ge =c(31,32,40,41,42,51,52)
E15_mid =c(33,34,43,44,45,53)
E15_hid =c(35,36,46,47,48,54,55)
passage <- list(E13_ctx,E13_ge,E13_mid,E13_hid,E15_ctx,E15_ge,E15_mid,E15_hid)
passage_name <- c("E13_ctx","E13_ge","E13_mid","E13_hid","E15_ctx","E15_ge","E15_mid","E15_hid")
topbar <- function(x){      
  return(mean(x)+sd(x)/sqrt(length(x))) #误差采用了mean+-sem
}
bottombar <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
for (c in 1:8) {
  
  if(c<5) {
    counts_region <- normalized_counts[which(rownames(normalized_counts)%in%regional_marker_name[[c]]),]
    counts_region_log <- log2(counts_region+1)
  } else {
    counts_region <- normalized_counts[which(rownames(normalized_counts)%in%regional_marker_name[[c-4]]),] 
    counts_region_log <- log2(counts_region+1)
  }
  counts_region_E <- counts_region_log [,passage[[c]]]
  counts_region_E_t <- as.data.frame( t(counts_region_E))
count_line <- data.frame()
for (i in colnames(counts_region_E_t)) {
  temp <- as.data.frame(counts_region_E_t[,i])
  colnames(temp) <- "Expression"
  temp$gene <- i
  temp$time <- coldata_1[rownames(counts_region_E_t),1]
  count_line <- rbind(count_line,temp)
}
count_line$time <-as.numeric(count_line$time)
count_line[,3][count_line[,3] == 1] = 0
count_line[,3][count_line[,3] == 2] = 21
count_line[,3][count_line[,3] == 3] = 42
passage_line <- ggplot(count_line,aes(time,Expression, color = gene,group =gene))+
  stat_summary(geom = 'line',fun='mean',cex=1)+
  stat_summary(geom = 'errorbar',
               fun.min = bottombar,fun.max = topbar,
               width=1,cex=0.5,aes(color=gene))+
  stat_summary(geom = 'point',fun='mean',aes(fill=gene),
               size=3,pch=21,color='black')+
  theme_classic(base_size = 15)
  #theme(legend.position = 'none')
png (filename = paste0("E:/Sequence/downstream analysis/DEG/line_plot/",
                       passage_name[c],"_passage_line.png"), width = 1600,height = 1200,res =300)
print(passage_line)
dev.off()

}

###### multibox plot of marker genes
counts_ctx_E13_P0 <- counts_ctx [,1:8]
counts_ctx_E13_P0_t <- as.data.frame( t(counts_ctx_E13_P0))
count_box <- data.frame()
for (i in colnames(counts_ctx_E13_P0_t)) {
  temp <- as.data.frame(counts_ctx_E13_P0_t[,i])
  colnames(temp) <- "counts"
  temp$gene <- i
  temp$group <- coldata_1[rownames(counts_ctx_E13_P0_t),"tissue"]
  count_box <- rbind(count_box,temp)
}
multi_box <- ggplot(count_box,aes(x=gene, y = counts)) +
                    geom_boxplot(aes(fill = group),position = position_dodge(0.5),width = 0.6)
multi_box


##############  PCA plot #################
vstdata <-vst(dds,blind = FALSE) 
assay(vstdata) <- limma::removeBatchEffect(assay(vstdata),vstdata$batch) #去批次效应
plotPCA(vstdata_2,intgroup = ("condition"))#package function
vst_plot<-plotPCA(vstdata,intgroup = ("condition"),returnData =T)
plot(vst_plot[,1:2],pch=19,col= c(rep("red",3),rep("blue",3),rep("green",3),rep("yellow",3)),cex=1.5,
     main = "The PCA analysis of P3-E15.5",xlab ="PC1: 44% variance",ylab="PC2: 22% variance",
     col.main ="black",xlim=c(-15,15))
text(vst_plot[,1],vst_plot[,2],row.names(vst_plot),cex=1, font = 1,adj = 1)
legend("topright",inset=.05,title="group",c("ctx","ge","mid","hid"),col = c("red","blue","green","yellow"),
       pch = 19,cex = 0.75,bty = "n")


library(FactoMineR)
library(factoextra)
library(rvest)
library(tidyverse)

dev.off()

vstdata_2_t <- as.data.frame(t(vstdata_2))
vst_pca <- PCA(vstdata_2_t,scale.unit = F,graph = T)
fviz_pca_ind(vst_pca,repel = T)

####### sample to sample distance heatmap #############
library(RColorBrewer)
sampleDist_vstdata <-dist(vstdata_2_t) #vst transformation
sampleDist_vstdata_matrix <- as.matrix(sampleDist_vstdata)
rownames(sampleDist_vstdata_matrix)<-paste(colnames(vstdata),sep = "-")
colnames(sampleDist_vstdata_matrix)<- paste(colnames(vstdata),sep = "-")
colors<- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDist_vstdata_matrix,
         clustering_distance_rows = sampleDist_vstdata,
         clustering_distance_cols = sampleDist_vstdata,
         col =colors)


#######  go enrichment analysis  ###########
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

module_gene <- read.csv(file = "E:/Sequence/downstream analysis/WGCNA/wgcna/geneInfo.csv",header = T)
rownames(module_gene) <- module_gene$X
module_gene <- module_gene[,-1]
module_turquoise <- module_gene[module_gene$moduleColor == "turquoise",1]

gene_module <- bitr(module_turquoise,fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)

ego_cc <- enrichGO(gene = gene_module$ENSEMBL,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENSEMBL",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
ego_bp <- enrichGO(gene = gene_module$ENSEMBL,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENSEMBL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

barplot(ego_cc,showCategory =18 ,title = " The GO_CC enrichment analysis of turquoise module") +
  scale_size(range = c(2,12)) +
  scale_x_discrete(labels = function(ego_cc) str_wrap(ego_cc,width = 25))

barplot(ego_cc,showCategory =20)
write.csv(ego_cc@result,file = "The GO_CC enrichment analysis of turquoise module.csv")

barplot(ego_bp,showCategory =20)
write.csv(ego_bp@result,file = "The GO_BP enrichment analysis of turquoise module.csv")

kk <- enrichKEGG(gene = gene_module$ENTREZID,
                 organism = "mmu",
                 pvalueCutoff = 0.05)

kk[1:30]
barplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")
 
dotplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")+
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))








