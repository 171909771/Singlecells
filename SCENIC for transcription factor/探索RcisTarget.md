- https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247497365&idx=1&sn=051651c84403a9ce6a7411b639a32220&scene=21#wechat_redirect
- https://www.jianshu.com/p/6e1d71db4220
# 读取示例基因
```r
library(RcisTarget)

geneList1 <- read.table(file.path(system.file('examples', package='RcisTarget'), "hypoxiaGeneSet.txt"), 
                        stringsAsFactors=FALSE)[,1]
geneLists <- list(hypoxia=geneList1 ) 
geneLists
```r

# 得到转录因子对应的motif。Motifs - Version 8 (mc8nr): 20003 motifs # Motifs - Version 9 (mc9nr): 24453 motifs
data("motifAnnotations_mgi")
motifAnnotations_hgnc 
# 把每个基因的TSS 周围10kb的序列用来预测motif排序，同一motif可以对应多个TF
motifRankings <- importRankings("cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.feather")
# 整合差异基因的motif
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)

# cisTarget包含了
motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
auc <- getAUC(motifs_AUC)[1,]
hist(auc, main="hypoxia", xlab="AUC histogram",
     breaks=100, col="#ff000050", border="darkred")
nes3 <- (3*sd(auc)) + mean(auc) #谨慎，取的3倍sd，p<0.01
abline(v=nes3, col="red")




