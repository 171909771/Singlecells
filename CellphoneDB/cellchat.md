### 参考
- https://www.jianshu.com/p/b3d26ac51c5a  原版，不能复制
- https://developpaper.com/cell-interaction-in-single-cell-analysis-3-cellchat/ 英文版，可以复制


```r
Sys.setenv(RETICULATE_PYTHON="/home/chan87/miniconda3/bin/python3") 
library(Seurat) 
library(SeuratData) 
library(tidyverse) 
library(CellChat) 
library(NMF) 
library(ggalluvial) 
library(patchwork) 
library(ggplot2) 
library(svglite) 
options(stringsAsFactors = FALSE) 


- https://www.jianshu.com/p/5d6fd4561bc0  准备pbmc.rds文件

pbmc3k.final <- readRDS("pbmc.rds") 

cellchat <- createCellChat(object=pbmc3k.final,group.by = "cell_type") 
#cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "cell_type") cellchat 
summary(cellchat) 
str(cellchat) 
levels(cellchat@idents) 
#cellchat <- setIdent(cellchat, ident.use = "cell_type") 
groupSize <- as.numeric(table(cellchat@idents))


# 导入cellchat数据库
CellChatDB <- CellChatDB.human 
#导入小鼠是CellChatDB <- CellChatDB.mouse
str(CellChatDB)
#查看数据库信息 #包含interaction、complex、cofactor和geneInfo这4个dataframe
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4] 
head(CellChatDB$cofactor) 
head(CellChatDB$complex)
head(CellChatDB$geneInfo) 
#dev.new() #下一步不出图的时候运行 
showDatabaseCategory(CellChatDB)

## 用特定的数据库进行分析
unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种 
#选择"Secreted Signaling"进行后续细胞互作分析 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis 
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB 
# set the used database in the object 
cellchat@DB <- CellChatDB.use # set the used database in the object

## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个） 
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4) 
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB 
#上一步运行的结果储存在cellchat@LR$LRsig 
cellchat <- projectData(cellchat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project

## 计算
#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。 # Filter out the cell-cell communication if there are only few number of cells in certain cell groups cellchat <- filterCommunication(cellchat, min.cells = 10) 
df.net <- subsetCommunication(cellchat) 
write.csv(df.net, "net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat) 
df.netp <- subsetCommunication(cellchat, slot.name = "netP") 


#### 展示图片
## 总体的互作关系
cellchat <- aggregateNet(cellchat)
#Calculate the number of each cell
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")

## 每种细胞的互作关系

mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
par("mar")
par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

## 每种细胞的互作强度
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


## 单信号通路作图
cellchat@netP$pathways # 查看细胞通路
pathways.show <- c("TGFb")

###  层次图
levels(cellchat@idents)    # show all celltype
vertex.Receiver = c(1,2,6) # define a numeric vector giving the index of the cell type as targets
par(mar=c(5.1,4.1,4.1,2.1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver  = vertex.Receiver,layout="hierarchy" )

### 圈圈图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",  vertex.receiver  = vertex.Receiver)

### 玹图: 不能用，报错
par(mfrow=c(1,1))
vertex.Receiver = c(2,5,6)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord",vertex.receiver  = vertex.Receiver)

### 热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

## 在指定通路中贡献值最大的受体-配体对
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TGFb = extractEnrichedLR (cellchat, signaling = pathways.show, geneLR.return = F) 


### 分层图
LR.show <- pairLR.TGFb[1,] 
vertex.receiver = c(1,2,4,6) # 不画分层就可以不用
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver, layout="hierarchy")

### 圈圈图
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#### 自动出所有图，没有试过
if(F){
#Access all the signaling pathways showing significant Communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex. Receiver = C (1,2,4,6) # this step is not necessary without drawing a hierarchy diagram
dir. Create ("all_paths_com_circle") # create a folder to save batch drawing results
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
            vertex. receiver = vertex. Receiver, layout = "circle") # draw network diagram
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}
setwd("../")
}

## 最重要的气泡图
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#Receptor cells and ligand cells need to be specified
netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), 
                     targets.use = c(1,2,4,6), remove.isolate = FALSE)

### 指定通路作图
netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), targets.use = c(1,2,4,6), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)


### 指定通路和细胞
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","TGFb"))
netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(1,4,2),  # 按照levels(cellchat@idents)的细胞标记，source是ligand，target是recepter
                 pairLR.use = pairLR.use, remove.isolate = TRUE)

### 查看表达量
plotGeneExpression(cellchat, signaling = "TGFb")
plotGeneExpression(cellchat, signaling = "TGFb", type = "dot",color.use="red")
```

