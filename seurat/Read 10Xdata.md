# if just read a matrix from GEO 3files (barcode, matrix and feature)
```r
library(Matrix)
matrix_dir = "~/filtered_feature_bc_matrix/hg19/"   ##根据实际文件夹进行修改
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
```
# using read10x() to read single sample
```r
library(Seurat)
getwd() # get current path
data_dir <- 'D:/test/GSE174574_RAW/1'
list.files(data_dir)
data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = data)
```
# read multiple samples
## error codes needed to modify, the purpose is that each sample has an own folder containing 3 files (barcode, features and matrix),
```r
if(F){
fs=list.files('./','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]

lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste0("GSE174574/", str_split(y[1],'_',simplify = T)[,1])
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE174574/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE174574/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE174574/",y[3]),file.path(folder,"matrix.mtx.gz"))
})}
```

## after creating a  set folder containing each sample folder
```
library(Seurat)
samples=list.files("GSE174574/")
samples
dir <- file.path('./GSE174574',samples)
names(dir) <- samples
```
## 合并方法1 全部合并
```r
counts <- Read10X(data.dir = dir)
scRNA1 = CreateSeuratObject(counts, min.cells=3)
dim(scRNA1)   #查看基因数和细胞总数
table(scRNA1@meta.data$orig.ident)  #查看每个样本的细胞数
```

## 合并方法2 分别建立后合并
```r
scRNAlist <- list()
for(i in 1:length(dir)){
  print(i)
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=3)
}
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]]))
dim(scRNA2)   #查看基因数和细胞总数
table(scRNA2@meta.data$orig.ident)  #查看每个样本的细胞数
```
# 最简单的读取，单个，多个都适用
## 创建文件夹，文件夹里面放这个数据集
```r
library(Seurat)
samples=list.files("MCAO/")
samples
dir <- file.path('./MCAO',samples)
names(dir) <- samples
```

# 读取h5文件
data2=Read10X_h5("filtered_feature_bc_matrix.h5")
