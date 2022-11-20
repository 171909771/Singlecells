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
### AUCell_exploreThresholds函数进行，加载了下面的代码
- https://github.com/aertslab/AUCell/blob/bf2b9e19b5b7a5736b12d118fa88382e8101146e/R/priv_auc.assignmnetThreshold_v6.R

**思路： 对模型进行binomal分布拟合，分析后选取以下最佳阈值：**
- minimumDens： 分析核密度曲线，取第一个上升的拐点
```r
ensCurve <- density(auc, adjust=densAdjust, cut=0)
  maximumsDens <- NULL
  inflPoints <- diff(sign(diff(densCurve$y)))
  maximumsDens <- which(inflPoints==-2)
  globalMax <- maximumsDens[which.max(densCurve$y[maximumsDens])]
  minimumDens <- which(inflPoints==2)
  smallMin <- NULL
  if(!skipSmallDens)
    smallMin <- data.table::last(minimumDens[which(minimumDens < globalMax)]) #1prev to max
  minimumDens <- c(smallMin,
        minimumDens[which(minimumDens > globalMax)]) # all after maximum

  # Density-based threshold (V4):
  # First minimum after the biggest maximum   (adjust=2)
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
