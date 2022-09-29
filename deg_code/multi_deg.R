rm(list=ls())
options(stringsAsFactors = F)
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
library(BiocManager)
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("githubinstall")
BiocManager::install("openxlsx")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(githubinstall)
library(devtools)
library(openxlsx)
devtools::install_github("BioinformaticsFMRP/TCGAbiolinksGUI.data")
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
query <- GDCquery(project = "TCGA-SARC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")
GDCdownload(query)
SARC_DATA<-GDCprepare(query,save=FALSE)
save(SARC_DATA,file='sarc.Rdata')
summary(SARC_DATA)
sarc_data<- assay(SARC_DATA)
dim(sarc_data)
clinic <- read.xlsx('NIHMS912614-supplement-10.xlsx',sheet = 2,startRow = 2)
colnames(sarc_data)<- substr(colnames(sarc_data),start = 1,stop = 15)
sarc_choose<- sarc_data[,clinic$TCGA.barcode]
save(sarc_choose,clinic,file='ex_pd.Rdata')
library(edgeR)
#####过滤掉cpm<2的counts
tmp<- cpm(sarc_choose)
tmp1<- tmp[rowSums(tmp<2)>1,]
sarc_choose<- sarc_choose[rownames(tmp1),]
identical(colnames(sarc_choose),clinic$TCGA.barcode)
####用TMM方法进行normalization
dge <- DGEList(counts=sarc_choose)
dge <- calcNormFactors(dge,method = 'TMM')
####构建对应的group_list
load(file = "ex_pd.Rdata")
subtype <- unique(clinic$short.histo)
group_list <- c()
for( i in 1:length(subtype)){
  tmp<- ifelse(clinic$short.histo==subtype[i],subtype[i],'other')
  group_list<- cbind(group_list,tmp)
  colnames(group_list)[i] <- subtype[i]
}
#####差异分析和PCA 
deg_re <- matrix()
  design <- model.matrix(~0+factor(group_list[,1]))
  colnames(design)=levels(factor(group_list[,1]))
  rownames(design)=colnames(sarc_choose)
  data_v <- voom(dge,design)
  ###PCA
  library("FactoMineR")
  library("factoextra") 
  ####
  constrsts<- paste0(colnames(group_list)[1],"-other")
  contrast.matrix<-makeContrasts(contrasts=constrsts,levels = design)
  ######limma三部曲,只需要归一化后的数据、实验矩阵、比较矩阵的输入
  fit <- lmFit(data_v,design)
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit1, 0.01)
  ####上面得到了p值等统计的结果，topTable对p值校验，对基因排序
  tT_tmp <- topTable(fit2, adjust="fdr", number=nrow(fit2))
  tT_tmp<- subset(tT_tmp, select=c('adj.P.Val',"P.Value","logFC"))
  colnames(tT_tmp)<- paste0(colnames(group_list)[1],colnames(tT_tmp))
  deg_re <- cbind(deg_re,tT_tmp)
deg_re <- deg_re[,-1]
save(deg_re,file='deg.Rdata')


