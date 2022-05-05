- https://mp.weixin.qq.com/s?__biz=MzUzMTEwODk0Ng==&mid=2247493674&idx=1&sn=c5f5272b9448787542dc52d1410e91e7&scene=21#wechat_redirect
- https://www.jianshu.com/p/c15db7c9d3c5

# 读取GEO数据
dat1=read.table("GSE60361_C1-3005-Expression.txt.gz",header = TRUE,row.names=1)
dat2=dat1[!duplicated(dat1$cell_id),]
rownames(dat2)=dat2$cell_id
dat2$cell_id=NULL
mouseBrainExprMatrix=as.matrix(dat2)
## 随机取5000个基因
exprMatrix <- mouseBrainExprMatrix[sample(rownames(mouseBrainExprMatrix), 5000),] 

# 开始做AUCell
library(AUCell)
## 按照基因表达排序
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE) 
## 创建基因集
library(GSEABase)
gmtFile <- paste(file.path(system.file('examples', package='AUCell')), "geneSignatures.gmt", sep="/")
gmtFile
geneSets <- getGmt(gmtFile)
geneSets
names(geneSets) 
geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, 
                            newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep="")
)
geneSets
names(geneSets)
### 创建随机基因集
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
### 创建Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), 
                           setName="HK-like (100g)"))
geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)
# 计算AUC值
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_AUC
names(geneSets)
## 绘制AUC图
par(mfrow=c(3,3))
test1=AUCell_exploreThresholds(cells_AUC)

