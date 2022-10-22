不同分布的算法
------
```
# 离散分布 https://stats.oarc.ucla.edu/other/dae/
WCL2=zeroinfl(y1~x1,dist = 'negbin')  # 零膨胀分布 ,注意：一定要检查，数值不能全为负  library(pscl)
WCL=glm.nb(y1~x1)    # 负二项分布  library(MASS)
# 连续分布
WCL1=glm(y1~x1,family=tweedie(var.power=1.5,link.power=0))   # tweedie分布:解释 - https://freakonometrics.hypotheses.org/56493 library(statmod)
y=cpglm(y1~x,link = 'log') #tweedie分布,自动寻找var.power  library(cplm, quietly = TRUE)

```

线性回归中的link function的推导
------
- https://www.r-bloggers.com/2018/10/generalized-linear-models-understanding-the-link-function/
