- https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247497365&idx=1&sn=051651c84403a9ce6a7411b639a32220&scene=21#wechat_redirect
- https://www.jianshu.com/p/6e1d71db4220
# 读取示例基因
```r
library(RcisTarget)

geneList1 <- read.table(file.path(system.file('examples', package='RcisTarget'), "hypoxiaGeneSet.txt"), 
                        stringsAsFactors=FALSE)[,1]
geneLists <- list(hypoxia=geneList1 ) 
geneLists
```

# 得到转录因子对应的motif。Motifs - Version 8 (mc8nr): 20003 motifs # Motifs - Version 9 (mc9nr): 24453 motifs
```r
data("motifAnnotations_mgi")
motifAnnotations_hgnc 
```
# 把每个基因的TSS 周围10kb的序列用来预测motif排序，同一motif可以对应多个TF
```r
motifRankings <- importRankings("cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.feather")
```
# 整合差异基因的motif,可分为以下3steps
```r
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)
```

##1 计算AUC值
motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
### 验证AUC的分布
auc <- getAUC(motifs_AUC)[1,]
hist(auc, main="hypoxia", xlab="AUC histogram",
     breaks=100, col="#ff000050", border="darkred")
nes3 <- (3*sd(auc)) + mean(auc) # 取的3倍SD，p<0.01
abline(v=nes3, col="red")

##2 注释motif（转换成TF），并加NES得分
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
motifAnnot=motifAnnotations_hgnc)

##3 添加富集的基因，加参数method=“iCisTarget”更为精确
motifEnrichmentTable_wGenes2 <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=geneLists
                                                   )
### 富集曲线图
selectedMotifs <-sample(motifEnrichmentTable$motif, 3)
par(mfrow=c(2,2))
getSignificantGenes(geneLists, 
                    motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="aprox")

### motif AGCT图
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)

resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]

library(DT)
datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))


# 画网络图
signifMotifNames <- motifEnrichmentTable$motif[1:5]
incidenceMatrix <- getSignificantGenes(geneLists, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix
library(igraph)
mugh<- graph_from_incidence_matrix(incidenceMatrix, directed = F)

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]


length(layouts)
par(mfrow=c(3,5), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(mugh)) 
  plot(mugh, edge.arrow.mode=0, layout=l, main=layout) }







