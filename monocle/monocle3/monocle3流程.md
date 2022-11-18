- https://cole-trapnell-lab.github.io/monocle3/docs/installation/
### 载入数据
```r
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
dat1 <- readRDS('prunnedfinalcluster.rds')
dat1 <- subset(dat1,tmp=='Microglia')
dat1 <- subset(dat1,celltype.group=='MCAO')
```
### 载入monocle矩阵及数据前处理
```r
DefaultAssay(object = dat1 ) <- "RNA"
expr_matrix <- as.matrix(dat1@assays[["RNA"]]@counts)

pd <- dat1@meta.data
fd <- data.frame(gene_short_name= rownames(expr_matrix) , row.names = rownames(expr_matrix))
cds <- new_cell_data_set(expr_matrix, 
                          cell_metadata  = pd, 
                          gene_metadata  = fd)

## 数据前处理
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
```
### 分群
```r
cds <- cluster_cells(cds)
### 绘图
plot_cells(cds, color_cells_by = 'partition',group_label_size = 5)
```
![image](https://user-images.githubusercontent.com/41554601/202502404-4300d985-891e-43f6-8daa-50adcb630af6.png)


### 绘制路径
```r
## use_partition参数：F表示整张图一个路径；close_loop：F表示不闭环.参考下面的说明
## https://github.com/cole-trapnell-lab/monocle3/issues/328
cds <- learn_graph(cds,use_partition = F,close_loop =F)
## 合并多个partition
## https://github.com/cole-trapnell-lab/monocle3/issues/328

### 绘图
plot_cells(cds,
           color_cells_by = "partition",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```
![image](https://user-images.githubusercontent.com/41554601/202502445-fc2fc2d3-ab7a-4214-b8b5-fa4af3999ba2.png)

### 决定root nodes
```r
## 利用meta信息中的小cluster（time_bin）作为标识，[, "seurat_clusters"]可以目标信息
get_earliest_principal_node <- function(cds, time_bin="31"){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

## 绘图
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```
![image](https://user-images.githubusercontent.com/41554601/202502546-a2af039a-57eb-4fdb-9859-6eb4c241841b.png)

##### 用连续数值确定轨迹点
```r
## 取每个细胞中P2的表达矩阵
tmp <-cds@assays@data@listData[["counts"]]['P2ry12',]
tmp <- data.frame(tmp)
tmp$cellname <- rownames(tmp)
## 取每个细胞中的路径点值
closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- data.frame(closest_vertex)
closest_vertex$cellname <- rownames(closest_vertex)
## 融合2个矩阵
tmpmerg <- merge(tmp,closest_vertex,by='cellname')
## 计算在每个轨迹点的平均值，并取最大值的点
tmpmerg <- data.frame(tmpmerg)
rownames(tmpmerg) <- tmpmerg$cellname
tmpmerg$cellname <- NULL
test <- tmpmerg %>% group_by(closest_vertex) %>% summarise('meannumber'=mean(tmp))
test.numb <- test$closest_vertex[which.max(test$meannumber)]
root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(test.numb)]
cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
```

### 划分modules
```r
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=cds@clusters@listData[["UMAP"]][["clusters"]])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
```
![image](https://user-images.githubusercontent.com/41554601/202597203-b1135bc5-48c2-4f7c-ba79-42ebd1a655e8.png)

### 颜色改变
```r
library(tidyverse)
## 改变颜色
## https://github.com/cole-trapnell-lab/monocle3/issues/181
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(60,44)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)+scale_color_gradient(low="green", high="red")
```
![image](https://user-images.githubusercontent.com/41554601/202597190-76e34a57-8540-4441-a2b2-a0a10efcc2ce.png)

### 手动选择目标cluster后绘制module
```r
## 手动选cluster
cds_subset <- choose_cells(cds)
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
## 确定展现的基因
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])
## 每个基因出轨迹图
plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
```
