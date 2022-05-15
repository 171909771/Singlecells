- https://www.jianshu.com/p/1d9189db05b0 #示例
- http://events.jianshu.io/p/b662ffa45120  # 解读箭头
- https://www.cnblogs.com/raisok/p/12425258.html
library(devtools)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(hdf5r)
library(dplyr )

ldat <- ReadVelocity(file = "SCG71.loom")
bm <- as.Seurat(x = ldat)


bm <-SCTransform(object = bm, assay = "spliced") %>% RunPCA(verbose = FALSE) %>% 
  FindNeighbors(dims = 1:20) %>% FindClusters() %>% RunUMAP(dims = 1:20)
# 降维聚类分群都是基于"spliced"这个assay做的
DimPlot(bm)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)


ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)

## arrow.scale为箭头密度
show.velocity.on.embedding.cor(emb = bm@reductions[["umap"]]@cell.embeddings, vel = Tool(object = bm, 
                                                                                         slot = "RunVelocity"), n = 300, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
