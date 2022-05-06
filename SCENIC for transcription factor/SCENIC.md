- https://m.sciencenet.cn/home.php?mod=space&uid=118204&do=blog&id=1208136   # 官方中文
- https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247498091&idx=1&sn=5c157da6e889f2134ef9342cee560197&chksm=9b4bb9d0ac3c30c60adf722ae7acc54bb65a030163c0f43150a3cf874f46a390bda4f189aaeb&scene=178&cur_album_id=1909628995961159685&ascene=0&devicetype=android-30&version=28001257&nettype=WIFI&abtest_cookie=AAACAA%3D%3D&lang=zh_CN&exportkey=AR729ICKOPWdGRTIas05ESA%3D&pass_ticket=JNlhFADLUv85a7%2BD9vOnCFNP7%2FYDDuNdeyQslGY8d3HJN5UK7Sn3Y1bsozDXA75z&wx_header=3 # jimmy教程


library(SCENIC)
library("doRNG")
library("doMC")
library("DT")
library(NMF)
library(ComplexHeatmap)
library(R2HTML)
library(Rtsne)
library(zoo)
library(mixtools)
library(rbokeh)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(doParallel)

# Load data
loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
library(SCopeLoomR)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)
## 保存分组信息
saveRDS(cellInfo, file="int/cellInfo.Rds")

# 下载数据
- https://resources.aertslab.org/cistarget/
if(F){
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
}
# 保证 cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", nCores=1)
## 添加分组信息，一定用这种方法才能在热图中添加
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Co-expression network
## 过滤基因
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

# SCENIC
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)


# 分组转录因子
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
