不同分布的算法
------
```
library(MASS)
library(tweedie)
library(statmod)
library(pscl)
WCL2=zeroinfl(y1~x1,dist = 'negbin')  # 零膨胀分布
WCL1=glm(y1~x1,family=tweedie(var.power=1.5,link.power=0))   # tweedie分布
WCL=glm.nb(y1~x1)    # 负二项分布
```
