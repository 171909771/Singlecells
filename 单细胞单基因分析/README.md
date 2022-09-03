不同分布的算法
------
```
WCL2=zeroinfl(y1~x1,dist = 'negbin')  # 零膨胀分布 ,注意：一定要检查，数值不能全为负  library(pscl)
WCL1=glm(y1~x1,family=tweedie(var.power=1.5,link.power=0))   # tweedie分布:解释 - https://freakonometrics.hypotheses.org/56493 library(statmod)
y=cpglm(y1~x,link = 'log') #tweedie分布,自动寻找var.power  library(cplm, quietly = TRUE)
WCL=glm.nb(y1~x1)    # 负二项分布  library(MASS)
```
