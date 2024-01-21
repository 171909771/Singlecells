### 整理单细胞数据集
### 每个数据集整理为以下格式
![image](https://github.com/171909771/Singlecells/assets/41554601/cf403827-13b6-4bb4-b1f0-613ea67e9b34)

```
fs=list.files('./','^GSM')
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]

filename1="GSE227651/"
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste0(filename1, str_split(y[1],'_',simplify = T)[,1])
  #为每个样本创建子文件夹
  dir.create(folder,recursive = T)
  #移动文件
  file.copy(y[1],folder)
  file.copy(y[2],folder)
  file.copy(y[3],folder)
  #重命名文件，注意改
  file.rename(paste0(folder,'/',y[1]),paste0(folder,'/',str_split(y[1],'-',simplify = T)[,2]))
  file.rename(paste0(folder,'/',y[2]),paste0(folder,'/',"features.tsv.gz"))
  file.rename(paste0(folder,'/',y[3]),paste0(folder,'/',str_split(y[3],'-',simplify = T)[,2]))
})

```

