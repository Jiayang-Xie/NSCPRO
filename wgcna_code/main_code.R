options(repos="https://mirror.lzu.edu.cn/CRAN/")
library(BiocManager)
BiocManager::install("WGCNA")
BiocManager::install("sva")
library(WGCNA)
library(DESeq2)
library(stringr)
library(dplyr)
library(ggplot2)

devtools::install_github("zhangyuqing/sva-devel")

######  0.data preparation for preparation (the expression data and the phenotypic data)

rawcount <-read.csv(file = "E:/Sequence/downstream analysis/DEG/all_samples_counts_by_deseq2.csv",row.names = 1)
coldata <- read.csv(file = "E:/Sequence/downstream analysis/DEG/all_samples_traits.csv",row.names = 1)
coldata <- coldata[-56,]
rownames(coldata) <- colnames(rawcount)
coldata$ctx  <- c(0)
coldata$ge  <- c(0)
coldata$mid  <- c(0)
coldata$hid  <- c(0)
coldata[,5][coldata[,2] == "ctx"] = 1
coldata[,6][coldata[,2] == "ge"] = 1
coldata[,7][coldata[,2] == "mid"] = 1
coldata[,8][coldata[,2] == "hid"] = 1
colnames(coldata)[1] <- "passage"
coldata[,1][coldata[,1] == 0] = 1
coldata[,1][coldata[,1] == 21] = 2
coldata[,1][coldata[,1] == 42] = 3
coldata <- coldata[,-c(2,4)]
write.csv(coldata,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_trait_wgcna.csv")

coldata1 <- read.csv(file = "E:/Sequence/downstream analysis/DEG/all_samples_traits.csv",row.names = 1,stringsAsFactors = T)
coldata1 <- coldata1[-56,]
coldata1$condition <- as.factor(coldata1$condition)
coldata1$stage <- as.factor(coldata1$stage)
rownames(coldata1) <- colnames(rawcount_1)
write.csv(coldata1,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_trait_wgcna.csv")
rawcount_1<- read.csv(file = "E:/Sequence/downstream analysis/DEG/all_samples_counts_gene_name.csv",row.names = 1)

### remove the batch effect

batch <- coldata1$batch
E13_P6 <- as.matrix(rawcount_1[,17:28])
E13_P6_adjust <- ComBat_seq(E13_P6,batch = coldata1[17:28,4],group = coldata1[17:28,2])
E13_P6_martix <- as.matrix(E13_P6)


dds <- DESeqDataSetFromMatrix(rawcount_1,coldata1,design = ~ stage+ tissue + condition )
dds<- DESeq(dds)
normalized_counts <-counts(dds,normalized=T)
normalized_counts<- na.omit(normalized_counts)
vstdata <-vst(dds,blind = FALSE)
vstdata_1 <- assay(vstdata)
write.csv(vstdata_1,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_count_genename_vst.csv")
vstdata_mad <- apply(vstdata_1, 1, mad)
vstdata_2 <- vstdata_1[order(vstdata_mad, decreasing=T), ]
vstdata_3 <- vstdata_2[1:5000, ]
write.csv(vstdata_2,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_count_genename_vst_mad.csv")
write.csv(vstdata_3,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_count_genename_vst_mad_5000.csv")

### extract the matrix without the batch effect
vstdata_4 <-vstdata_1[,-c(19,22,25,26,39,42,45,48)] 
coldata2 <- coldata[-c(19,22,25,26,39,42,45,48),]
write.csv(vstdata_4,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_count_genename_vst_nobatcheffect.csv")
vstdata_mad <- apply(vstdata_4, 1, mad)
vstdata_2 <- vstdata_4[order(vstdata_mad, decreasing=T), ]
vstdata_3 <- vstdata_2[1:5000, ]
write.csv(vstdata_2,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_count_genename_vst_mad_nobatcheffect.csv")
write.csv(vstdata_3,file = "E:/Sequence/downstream analysis/WGCNA/wgcna/E13_E15_sample_count_genename_vst_mad_5000_nobatcheffect.csv")

### extract the E13 matrix
vstdata_5 <- vstdata_1[,c(1:18,20,21,23,24,27,28)]
vstdata_mad <- apply(vstdata_5, 1, mad)
vstdata_2 <- vstdata_5[order(vstdata_mad, decreasing=T), ]
vstdata_3 <- vstdata_2[1:5000,]
coldata_5 <- coldata2[1:24,-2]


### 筛选中位数绝对偏差MAD前75%
vstdata_3 = vstdata_5[which(vstdata_mad > 
                              max(quantile(vstdata_mad,probs = seq(0,1,0.25))[2],0.01)),]



#### start the WGCNA
### 1 data input, wash, and preprocessing
## 1.1 检查缺失值和识别离群值（异常值）
dataExp <- t(vstdata_3)
gsg <- goodSamplesGenes(dataExp,verbose =3 )
gsg$allOK

# 1.1.1 如果结果为FALSE ，需要删除缺失值
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dataExp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dataExp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  dataExp = dataExp[gsg$goodSamples, gsg$goodGenes]
}

## 1.2 聚类所有样本，观察是否有缺失值
sampleTree <- hclust(dist(dataExp),method = "average")
sizeGrWindow(12,9)#试图
par(cex =0.6);
par(mar = c(0,4,2,0))
plot(sampleTree,main = "sample clustering to detect outliers",sub ="",xlab ='',
     cex.lab = 1.5,cex.axis = 1.5,cex.main =2)
## 1.3 cmdscale 降维处理 可视化
 sampleClu <- cmdscale(dist(dataExp),k = 2)
tibble(x = sampleClu[,1],
       y = sampleClu[,2],
      class = )


## 1.4 combine the expression data with the phenotypic data and reconstituted the sampletree

sampleTree2 = hclust(dist(dataExp),method = "average")
traitcolors = numbers2colors(coldata2,signed = F)# create a color representation for the expression and trait
plotDendroAndColors(sampleTree2,traitcolors,
                    groupLabels = names(coldata2),
                    main = "sample dendrogram and trait heatmap")
save(dataExp,coldata,file = "E11_E15_01WGCNA_datainput.RData")

### 2 construct the co-expression network
## 2.1 参数设置
getwd()
workingDir = "E:/Sequence/downstream analysis/WGCNA/wgcna"
setwd(workingDir)
enableWGCNAThreads() #开启多线程

## 2.2 挑选合适的软阈值
# 2.2.1 picksoftthreshold
powers = c(c(1:10),seq(from = 12, to =20,by =2))
sft = pickSoftThreshold(dataExp, networkType = "signed", powerVector = powers, verbose = 3)
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9

# 2.2.2 无标度拓扑拟合指数

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold(power)",
     ylab ="Scale Free Topology Model Fit, signed R^2",
     type="n",
     main = "Scale independence");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,cex = cex1, col = "red");
abline(h =0.9,col="red")

power = sft$powerEstimate

#2.2.3 平均连接度

plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab ="Soft Threshold(power)",
     ylab= "Mean Connectivity",type ="n",
     main = paste("Mean connectivity")
     )
text(sft$fitIndices[,1],sft$fitIndices[,5],labels = powers,cex = cex1,col = "red")

# 2.2.4 如无合适软阈值，按以下条件选择
power = sft$powerEstimate
nSamples = 55
type = "signed"

if(is.na(power)){
  power = ifelse(nSamples <20, ifelse(type == "unsigned",9,18),
                 ifelse(nSamples < 30, ifelse(type == "unsigned",8,16),
                        ifelse(nSamples < 40, ifelse(type =="unsigned",7,14),
                               ifelse(type == "unsigned",6,12))
                        )
                 )
}

## 2.3 构建网络和模块检测（一步法）
cor <- WGCNA::cor
net = blockwiseModules(dataExp,power = power,
                       TOMType = "signed", minModuleSize = 30, 
                       reassignThreshold = 0,mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = F,
                       saveTOMs = T,
                       saveTOMFileBase = "E11_E15_TOM")

# minModuleSize 模块中最少的基因数
# mergeCutHeight ：模块合并阈值，阈值越大，模块越少（重要）
# saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM"保存TOM矩阵，名字为"femaleMouseTOM"

##2.4 模块标识的层次聚类树状图

sizeGrWindow(12,9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = F, hang =0.03,
                    addGuide = T, guideHang = 0.05)
## 2.5 保存分配模块和模块包含的基因信息
moduleLabels = net$colors
mergedColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs,moduleLabels,mergedColors,geneTree,
     file = "E13_15_02networkconstruction_auto.RData")

### 3 模块与表型数据关联并识别重要基因
## 3.0 参数设置和载入之前的结果
getwd()
workingDir = "."
options(stringsAsFactors = F)
lnames = load(file = "") # 01 01WGCNA_datainput
lnames = load(file = "") # 02 tworkconstruction_auto

## 3.1 模块-表型数据关联
nGenes = ncol(dataExp)
nSamples = nrow(dataExp)

#重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(dataExp,mergedColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs,coldata2,use = "p")
moduleTraitPvalue =corPvalueStudent(moduleTraitCor,nSamples)

#通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和p值
textMatrix = paste(signif(moduleTraitCor,2),"\n",
                   "(",signif(moduleTraitPvalue,3),")",sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

# 用热图的形式展示相关系数

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(coldata2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relations"))

## 3.2 基因与表型数据的关系、重要模块：基因显著性和模块成员
# 基因的显著性GS定义为基因与性状的相关性

passage = as.data.frame(coldata2$passage)
names(passage) = "passage"
modNames = substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(dataExp,MEs,use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership) = paste("MM",modNames,sep = "")
names(MMPvalue) = paste("p.MM",modNames,sep = "")
geneTraitSignificance = as.data.frame(cor(dataExp,passage,use = "p")) # the correlation with passage
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) = paste("GS.",names(passage),sep = "")


## 3.3 模块内分析：鉴定具有高GS和高MM的基因
# GS：所有基因表达谱与这个模块的eigengene的相关性(cor) 
# MM: 基因和表型性状之间相关性的绝对值
module = "turquoise"
column = match(module, modNames)
moduleGenes = mergedColors ==module
sizeGrWindow(7,7)
par(mfrow =c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab =paste( "Module Membership in", module,"module"),
                   ylab = "gene significance for passage",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2, cex.axis = 1.2, col = module)

## 3.4 输出网络分析结果

geneInfo0 = data.frame(substanceBXH = names(as.data.frame(dataExp)),
                                            moduleColor = mergedColors,
                                            geneTraitSignificance,
                                            GSPvalue)

# sort by the significance level with passage
modOrder = order(-abs(cor(MEs,passage,use = "p")))
# 添加模块成员的信息
for (mod in 1:ncol(geneModuleMembership)){
   oldnames = names(geneInfo0)
   geneInfo0  = data.frame(geneInfo0,geneModuleMembership[,modOrder[mod]],
                           MMPvalue[,modOrder[mod]])
   names(geneInfo0) =c(oldnames,paste("MM.",modNames[modOrder[mod]],sep = ""),
                       paste("P.MM",modNames[modOrder[mod]],sep = ""))
}
geneOrder = order(geneInfo0$moduleColor,-abs(geneInfo0$GS.passage))
geneInfo = geneInfo0[geneOrder,]
write.csv(geneInfo,file = "geneInfo.csv")

### 4 网络交互分析（go注释等）
## 4.1 输出基因列表
dataExp_module = as.data.frame(mergedColors,colnames(dataExp))
write.csv(dataExp_module,file = "dataExp_module.csv")



### 5 网络可视化
## 5.1 参数设置

nGenes = ncol(dataExp)
nSamples = nrow(dataExp) 

## 5.2 可视化基因网络
# 5.2.1 计算TOM矩阵

dissTOM = 1-TOMsimilarityFromExpr(dataExp,power = power)
plotTOM = dissTOM^7
diag(plotTOM) = NA
sizeGrWindow(9,9) 
TOMplot(plotTOM,geneTree,mergedColors,main = "Network heatmap plot, all genes")

# 5.2.2 选取基因子集进行构建树状图

nSelect = 400
set.seed(10) # 随机数种子
select_num = sample(nGenes,size = nSelect)
selectTOM = dissTOM[select_num,select_num]
selectTree = hclust(as.dist(selectTOM),method = "average")
selectColors = mergedColors[select_num]
sizeGrWindow(9,9)
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss,selectTree,selectColors,main = "Network heatmap plot,selected genes")

# 改变热图的深色背景为白色背景
library(gplots)
myheatcol = colorpanel(250,"red","orange","lemonchiffon")
TOMplot(plotDiss,selectTree,selectColors,main = "Network heatmap plot,selected genes",Colors = myheatcol)

# 5.2.3 可视化表征基因网络
# 重新计算模块的eigengenes
MEs = moduleEigengenes(dataExp,mergedColors)$eigengenes
# 提取感兴趣的表型数据
passage = as.data.frame(coldata2$passage)
names(passage) = "passage"
# 加入到相应的模块
MET =orderMEs(cbind(MEs,passage))
#画图
sizeGrWindow(5,7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET,"",marDendro = c(0,4,1,2),marHeatmap = c(3,4,1,2),cex.lab = 0.8, 
                      xLabelsAngle = 90)
# 拆分聚类图和热图
sizeGrWindow(6,6)
par(cex=1.0)
plotEigengeneNetworks(MET,"Eigengene dendrogram",marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET,"Eigengene adjacency heatmap",marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE,xLabelsAngle = 90)

### 6 将网络导出到网络可视化软件
## 6.1 输出到VisANt 软件所学的数据
TOM = TOMsimilarityFromExpr(dataExp,power = power)
module = "turquoise"
probes = names(dataExp)
inModule = (mergedColors == module)
modProbes = probes[inModule]
modTOM = TOM[inModule,inModule]
dimnames(modTOM) = list(modProbes,modProbes)
vis = exportNetworkToVisANT( modTOM,
                             file = paste("VisANTInput-",module,".txt",sep = ""),
                             weighted = TRUE,
                             threshold = 0,
                             )
# 6.1.1 控制输出hubgene 个数
nTOP = 30 
IMConn = softConnectivity(dataExp[,modProbes])
top = (rank(-IMConn) <= nTOP)
vis = exportNetworkToVisANT( modTOM[top,top],
                             file = paste("VisANTInput-",module,"-top30.txt",sep = ""),
                             weighted  = TRUE,
                             threshold = 0,
)

## 6.2 输出到Cytoscape
TOM = TOMsimilarityFromExpr(dataExp,power = power)
# 选择需要的模块
module = "turquoise"
probes =names(dataExp)
inModule = is.finite(match(mergedColors,module))
modProbes= probes[inModule]
modGenes = 
# 选择相关的TOM矩阵
modTOM =TOM[inModule,inModule]
dimnames(modTOM) = list(modProbes,modProbes)
# export the network into edge and node list files Cytoscape can read

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-",
                                                paste(module,collapse = "-"),
                                                ".txt",sep = ""),
                               nodeFile = paste("CytoscapeInput-nodes-",
                                                paste(module,collapse = "-"),
                                                ".txt",sep = ""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]
)




### addition 




