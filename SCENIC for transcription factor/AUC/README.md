## 官方介绍
- https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html

## 非官方介绍
- http://www.bio-info-trainee.com/7823.html

AUC的运算
--------
- 1.对所有基因按表达排序
- 2.用排序的基因与基因集进行匹配后提取（只保留匹配的基因，作为x轴），出现一个基因y轴就+1
- 3.取提排序前5%的基因进行分析，形成下图
- 4.计算AUC值（还没有研究）
![image](https://user-images.githubusercontent.com/41554601/202910492-042d535f-3c3d-4878-a642-89899d0027f5.png)

取阈值的方法
----
*AUCell_exploreThresholds函数进行，加载了下面的代码
- https://github.com/aertslab/AUCell/blob/bf2b9e19b5b7a5736b12d118fa88382e8101146e/R/priv_auc.assignmnetThreshold_v6.R

**思路： 对模型进行binomal分布拟合，分析后选取以下最佳阈值：**
- minimumDens (plot in Blue)： 分析核密度曲线，取第一个上升的拐点
```r
densCurve <- density(auc, adjust=densAdjust, cut=0)  # 构建和密度
  maximumsDens <- NULL
  inflPoints <- diff(sign(diff(densCurve$y)))  # 很巧的方法来取密度曲线的拐点
  maximumsDens <- which(inflPoints==-2)   # 取波峰
  globalMax <- maximumsDens[which.max(densCurve$y[maximumsDens])]
  minimumDens <- which(inflPoints==2)    # 取波谷
  smallMin <- NULL
  if(!skipSmallDens)
    smallMin <- data.table::last(minimumDens[which(minimumDens < globalMax)]) #1prev to max
  minimumDens <- c(smallMin,
        minimumDens[which(minimumDens > globalMax)]) # all after maximum

  # Density-based threshold (V4):
  # First minimum after the biggest maximum   (adjust=2)    就是取第一个波谷卫阈值，并且后面的波峰需要大于最高波峰的5%
  densTrh <- NULL
  if(length(minimumDens)>0) # && (!skipMinimumDens))
  {
    densTrh <- densCurve$x[min(minimumDens)]
    # Commented on V6
    # Only keep if it is a real inflextion point
    # (i.e. next max at least 5% of the global max)
    if(length(maximumsDens)>0)
    {
      nextMaxs <- maximumsDens[which(densCurve$x[maximumsDens] > densTrh)]
      if((max(densCurve$y[nextMaxs])/max(densCurve$y))<.05)
      {
        densTrh <- NULL
        # print(gSetName)
      }
    }
  }
```

- Global_k1 (plot in Grey): 对整个数据进行正态分布分析，取0.01的BH校正值对应x轴
```r
  meanAUC <- mean(auc)
  sdAUC <- sd(auc)
  maybeNormalDistr <- !suppressWarnings(     # 正态分布检验
    ks.test(auc, rnorm(max(100,length(auc)),mean=meanAUC, sd = sdAUC),
            alternative = "less")$p.value < .01)
  if(maybeNormalDistr){
    commentMsg <- paste0(commentMsg,
            "The AUC might follow a normal distribution (random gene-set?). ")
    skipGlobal <- FALSE

    # aucThrs["outlierOfGlobal"] <- meanAUC + 2*sdAUC
    aucThrs["outlierOfGlobal"] <- qnorm(1-(thrP/nCells), mean=meanAUC, sd=sdAUC)   # 关键步骤
```

- tenPercentOfMax：取auc最大值的前10%, 注意code中中文注释的2个条件
```r
  histogram <- hist(c(0, auc/max(auc)), breaks=100, plot=FALSE)$count
  if((sum(histogram[1:5]) / sum(histogram)) >= notPopPercent*.75) {   # hist出的数据中前5%的频率大于总的计数的（.75*.75），即太多0了
    skipGlobal <- FALSE
    skipRed <- TRUE
    skipSmallDens <- TRUE
  }
  if((sum(histogram[1:10]) / sum(histogram)) >= notPopPercent*.50) {   # hist出的数据中前5%的频率大于总的计数的（.75*.5）才能用。
    skipSmallDens <- TRUE
    skipGlobal <- FALSE
    # skipRed <- TRUE ?
    aucThrs["tenPercentOfMax"] <- max(auc)*.10
  }
```


**bimodal拟合后的做分布，注释k值是拟合的正态分布个数，是L_k2核R_k3的前处理**
```r
 na <- capture.output(distrs[["k2"]] <-
        tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=2, verb=FALSE),
        # With fast, if there are many zeroes, it fails quite often
        error = function(e) {
          return(NULL)
        }))

 na <- capture.output(distrs[["k3"]] <-
        tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=3, verb=FALSE),
        error = function(e) {
          return(NULL)
        }))
```

- L_k2 (plot in Red): 左分布，2个混合分布中，取左分布（均值最小）的右侧0.01的BH校正值对应x轴
```r
  if(!is.null(distrs[["k2"]]))
  {
    k2_L <- which.min(distrs[["k2"]][["mu"]]) # (sometimes the indexes are shifted)
    aucThrs["L_k2"] <- qnorm(1-(thrP/nCells),
                             mean=distrs[["k2"]][["mu"]][k2_L],
                             sd=distrs[["k2"]][["sigma"]][k2_L])
  }
```
- R_k3 (plot in Pink): 右分布，3个混合分布中，取右分布(均值最大)的左侧0.01的BH校正值对应x轴
```r
  if(!is.null(distrs[["k3"]]))
  {
    k3_R <- which.max(distrs[["k3"]][["mu"]]) # R: right distribution
    k3_R_threshold <- qnorm(thrP,
                            mean=distrs[["k3"]][["mu"]][k3_R],
                            sd=distrs[["k3"]][["sigma"]][k3_R])
    if(k3_R_threshold > 0) aucThrs["R_k3"] <- k3_R_threshold
  }
```






