
## arrow.scale为箭头密度
show.velocity.on.embedding.cor(emb = bm@reductions[["umap"]]@cell.embeddings, vel = Tool(object = bm, 
                                                                                         slot = "RunVelocity"), n = 300, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
