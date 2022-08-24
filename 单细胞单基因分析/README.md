不同分布的算法
------
```
library(MASS)
library(tweedie)
library(statmod)
library(pscl)
WCL2=zeroinfl(y1~x1,dist = 'negbin')  # 零膨胀分布 ,注意：一定要检查，数值不能全为负
WCL1=glm(y1~x1,family=tweedie(var.power=1.5,link.power=0))   # tweedie分布:解释 - https://freakonometrics.hypotheses.org/56493
WCL=glm.nb(y1~x1)    # 负二项分布
```
