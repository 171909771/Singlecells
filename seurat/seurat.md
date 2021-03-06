https://www.jianshu.com/p/a4a556578a02



# 读取数据
```r
library(pheatmap)
library(Seurat)
if(F){
        library(Seurat)
        # https://satijalab.org/seurat/mca.html
        # 构建 Seurat 需要的对象
        # 其中 min.cells 和 min.genes 两个参数是经验值
        sce <- CreateSeuratObject(counts = counts, 
                                  meta.data =meta,
                                  min.cells = 5, #每个基因最低表达的细胞量
                                  min.features  = 200, #每个细胞最低表达的基因量
                                  project = "sce")
        # 参考：https://github.com/satijalab/seurat/issues/668
}
```
# 正式开始！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
#  load data into seurat 
## 过滤细胞和基因，第一种方法
```{r}
counts <- Read10X(data.dir = "./")
sce <- CreateSeuratObject(counts = counts, 
                          meta.data =meta,
                          min.cells = 5, #每个基因最低表达的细胞量
                          min.features  = 200, #每个细胞最低表达的基因量
                          project = "sce")
```

## 过滤细胞和基因，第二种方法
### step1 先看每个细胞中的参数比例
```r
pbmc <- PercentageFeatureSet(sce, pattern = "^Ercc", col.name = "percent.Ercc") 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","percent.Ercc"), ncol = 3)
```
### step2 based on the above pic, set cutoff: nUMI >2000; ERCC<50%
```{r}
sce <- subset(sce, subset = nFeature_RNA > 2000 & percent.ERCC < 50)
```
### step3 clean genes; average < 1，可以不用，SCTransform可以过滤基因
```{r}
counts1 <- GetAssayData(object = sce, slot = "counts")
test1=apply(counts1,1, mean)
filtered_counts <- counts1[test1>1, ]
sce <- CreateSeuratObject(filtered_counts, meta.data = sce@meta.data)
```

if(F){
# QC and clean data
## https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
## 2 pics to remove cells
### pic 1
```{r}
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.ERCC")
plot2
```
### pic 2
```{r}
qcdata=sce@meta.data
library(ggplot2)
qcdata %>% 
        ggplot(aes(x=nCount_RNA , y=nFeature_RNA, color=percent.ERCC)) + 
        geom_point() + 
        scale_colour_gradient(low = "gray90", high = "black") +
        stat_smooth(method=lm) +
        scale_x_log10() + 
        scale_y_log10() + 
        theme_classic() +
        geom_vline(xintercept = 500) +
        geom_hline(yintercept = 250) +
        facet_wrap(~plate)
```
}

# 下面三步现在可以用SCTransform来做,用这种方法不能选择簇JackStraw命令
- http://events.jianshu.io/p/854f7493b293
## 上面没有做percent，这里需要做
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt") 
## run sctransform 
### 用这个函数可以用更高维度的PCA
- https://www.biostars.org/p/423306/  #上面的话的解答
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

## step1 normalization
```{r}
scenorm <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
```

if(F){
### get overlaping genes,note: this filtered method is different from that paper
```{r}
gs1=intersect(gs,rownames(scenorm@assays[["RNA"]]@counts))
#just do a heatmap
pheatmap(scenorm@assays[["RNA"]]@data[gs1,])
```
}
## step2 find highly variable genes
```{r}
scenormvc <- FindVariableFeatures(scenorm, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(scenormvc)
top10 <- head(VariableFeatures(scenormvc), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

## step3 scale data, may reduce avariable effect
```{r}
all.genes <- rownames(scenormvc)
scaledata1 <- ScaleData(scenormvc, features = all.genes)

if(F) #remove avariable effect
{
scaledata2 <- ScaleData(scenormvc, vars.to.regress = "percent.ERCC")
}

pheatmap(scaledata1@assays[["RNA"]]@scale.data[gs1,])
```

# if know the dim, skip below 2 steps

# linear dimensional reduction:PCA
```{r}
PCAdata1 <- RunPCA(scaledata1, features = VariableFeatures(object = scaledata1))
## observe pca result
### method1 text
print(PCAdata1[["pca"]], dims = 1:5, nfeatures = 5)
### method2 
VizDimLoadings(PCAdata1, dims = 1:2, reduction = "pca")
### method3
DimPlot(PCAdata1, reduction = "pca",group.by = "g")
### method4
DimHeatmap(PCAdata1, dims = 1:4, cells = 691, balanced = TRUE,ncol=2)
```
## 如果有批次效应，可以参考CCA，但是生物学差异也被抹平
- http://www.360doc.com/content/21/0805/00/76149697_989558460.shtml

## 用SCTransform就可以不用确定dim，后面直接用高PCA值（30，50），效果一样
```{r}
if(F){
ddim1 <- JackStraw(PCAdata1, num.replicate = 100)
ddim2 <- ScoreJackStraw(ddim1, dims = 1:20)
##choose dim by gap,method1
JackStrawPlot(ddim2, dims = 1:20)
##choose dim by gap,method2
ElbowPlot(ddim2)
}
##set cluster
# 分辨率用多少 https://www.jianshu.com/p/ef28c8e723be

cluster1 <- FindNeighbors(ddim2, dims = 1:20)
cluster2 <- FindClusters(cluster1, resolution = seq(0.5,1.2,by=0.1))
library("clustree")
clustree(cluster2)
# 可以把这个参数加入到seurat数据中
Idents(object = pbmc) <- "RNA_snn_res.0.3"

head(Idents(cluster2))
```

# non-linear dimensional reduction
```{r}
## umap
umap1 <- RunUMAP(cluster2, dims = 1:10)
DimPlot(umap1, reduction = "umap")
## tsne
tsne1 <- RunTSNE(cluster2, dims = 1:10)
DimPlot(tsne1, reduction = "tsne")
```
# 看分群后效果，以microglia为例
```r
genes_to_check = c("Tmem119", "Hexb", "Cx3cr1", "Sparc",
                    "P2ry12", "Cst3","Aif1")
DotPlot(tsne1, group.by = 'RNA_snn_res.0.9',
        features = unique(genes_to_check)) + RotatedAxis()
```


# Finding differentially expressed features 
```{r}
## find cluster 2 DEGs from all clusters,  min.pct regarded as minimum percentage in each cluster
cluster2.markers <- FindMarkers(tsne1, ident.1 = 2, min.pct = 0.25)
## find cluster 5 DEGs from cluster 0 and 3
cluster5.markers <- FindMarkers(tsne1, ident.1 = 3, ident.2 = c(0, 3), min.pct = 0.25)
## find all cluster DEGs
pbmc.markers <- FindAllMarkers(tsne1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
## find cluster 2 DEGs from all clusters,  test.use = "roc" means adding a AUC 
cluster0.markers <- FindMarkers(tsne1, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
## summarize top2 gene in each cluster
library(dplyr)
pbmc.markers %>% group_by(cluster) %>%slice_max(n = 2, order_by = avg_log2FC) -> top2markers
top2markers
#visual
## single gene was grouped by cluster
VlnPlot(tsne1, features = c("Lrrc15"), slot = "counts", log = TRUE)
## maping in tSNE
FeaturePlot(tsne1, features = top2markers$gene)
## 
pbmc.markers %>% group_by(cluster) %>%slice_max(n = 10, order_by = avg_log2FC) -> top10markers
### no right information
DoHeatmap(tsne1, features = top10markers$gene) + NoLegend()
```

# dividing a cluster interested
```{r}
Cluster2.only <- subset(x = tsne1, idents = c('2'))
```

# add label in tsne
## https://satijalab.org/seurat/articles/pbmc3k_tutorial.html






