options(repos="https://mirror.lzu.edu.cn/CRAN/")
library(BiocManager)
BiocManager::install("FactoMineR")
BiocManager::install("sva")

## PCA analysis with VST or rlog transformation

dataexpr <- rawcount_1
dataexpr_nor <- normalized_counts
