## 算不同细胞群、基因群对方差的的贡献值
#### 按照细胞文库大小，等分5份
```
cellcluster=as.numeric(cut(cellsum,quantile(cellsum,seq(0,1,1/5))))
cellcluster[is.na(cellcluster)]=5
cluster=lapply(1:5, function(x) a1[,cellcluster==x]) #
```
#### 用counts算方差百分比
```
var=lapply(cluster, function(x) top6[,colnames(x)]) #基因也分了6份，一第6份基因为例
sdeach=sapply(var,function(x) sd(apply(x,1,function(y) (mean(y))^2)))   #
sdper=sapply(sdeach,function(x) x/sum(sdeach))  
```
#### 用log-normalization算方差百分比
```
logvar=lapply(cluster, function(x) logtop6[,colnames(x)]) #
logsdeach=sapply(logvar,function(x) sd(apply(x,1,function(y) (mean(y))^2)))   #
logsdper=sapply(logsdeach,function(x) x/sum(logsdeach))  #
```
##### 画图
```
dat3=data.frame("value"=c(sdper,logsdper),"condition"=rep(c("h1","h2","h3","h4","h5"),2),"name"=c(rep(c("con","log"),each=5)))
library(ggplot2)
ggplot(dat3, aes(fill=condition, y=value, x=name)) + geom_bar(position="fill", stat="identity")
```
![image](https://user-images.githubusercontent.com/41554601/183298387-02444557-a7a5-4577-b792-36d0069b10d1.png)
