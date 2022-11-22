---
title: "alra"
author: "longlong"
date: '2022-11-22'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
- https://github.com/KlugerLab/ALRA/blob/master/alra_test.R

## step1 归一化，类似seurat的Lognormalization
```{r log,eval=F}
normalize_data(A)
```

## step2 计算K值
$$
M=U\Sigma V^*
$$
- U是m×m階酉矩陣；Σ是m×n階非負實數對角矩陣；而V*，即V的共軛轉置

Σ相当于是PCA，计算K值就相当于是计算第几个PCA

randomized.svd函数得到3个值，其中*d*就是Σ

- K值是我们要求的维度
- thresh是标准差
- noise_start是Σ方差停止波动的预估起始值
- q是迭代次数
- 后面的参数是高科技，用来加速（看原文怎么启动）
```{r K-value,eval=F}
choose_k <- function (A_norm,K=100, thresh=6, noise_start=80,q=2,use.mkl=F, mkl.seed =-1) {
  #  Heuristic for choosing rank k for the low rank approximation based on
  #  statistics of the spacings between consecutive singular values. Finds
  #  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
  #  is significantly different than spacings in the tail of the singular values.
  # 
  #
  # Args:
  #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
  #   K: Number of singular values to compute. Must be less than the
  #   smallest dimension of the matrix.
  #   thresh: Number of standard deviations away from the ``noise'' singular
  #   values which you consider to be signal
  #   noise_start : Index for which all smaller singular values are
  #   considered noise
  #   q : Number of additional power iterations
  #   use.mkl : Use the Intel MKL based implementation of SVD. Needs to be
  #             installed from https://github.com/KlugerLab/rpca-mkl.
  #   mkl.seed : Only relevant if use.mkl=T. Set the seed for the random
  #   generator for the Intel MKL implementation of SVD. Any number <0 will
  #   use the current timestamp. If use.mkl=F, set the seed using
  #   set.seed() function as usual.
  # Returns:
  #   A list with three items
  #       1) Chosen k
  #       2) P values of each possible k 
  #       3) Singular values of the matrix A_norm
  
  if (K > min(dim(A_norm))) {
    stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
  }
  if (noise_start >K-5) {
    stop("There need to be at least 5 singular values considered noise.\n")
  }
  noise_svals <- noise_start:K
  if (!use.mkl) {
    rsvd_out <- randomized.svd(A_norm,K,q=q)
  }else {
    rsvd_out <- randomized.svd(A_norm,K,q=q, method='rsvd-mkl', mkl.seed=mkl.seed)
  }
  diffs <- rsvd_out$d[1:(length(rsvd_out$d)-1)] - rsvd_out$d[2:length(rsvd_out$d)]  # 1
  mu <- mean(diffs[noise_svals-1])  # 2
  sigma <- sd(diffs[noise_svals-1])   # 3
  num_of_sds <- (diffs-mu)/sigma    # 4
  k <- max (which(num_of_sds > thresh))   #  5
  return (list( k=k,num_of_sds = num_of_sds,d=rsvd_out$d))
}
```
1. 得到相邻Σ的差值

2. 求预估Σ停止波动范围的均值

3. 求预估Σ停止波动范围的标准差

4. 将相邻Σ的差值Z-score

5. 取大于6个标准差的末尾项


## step3 imputing zero value
alra假设： 通过发现zero在重建后的矩阵中围绕重建矩阵的0对称分布，取最小值的绝对值就可以涵盖掉所有的0.

1. 用最佳K值重新计算alra，得到重建的矩阵

2. 取占比0.1%的数值的绝对值，小于该绝对值的数值全部转变为0

3. scale 数据
$$
\frac{x_{svd}-\mu _{cor}}{\sigma _{cor}} = \frac{x_{exp}-\mu _{ori}}{\sigma _{ori}}
$$
可变形为
$$
x_{exp}=\frac{x_{svd}*\sigma _{ori}}{\sigma _{cor}}-\frac{\mu _{cor}*\sigma _{ori}}{\sigma _{cor}}+\mu _{ori}
$$



```{r alra, eval=F}
alra <- function( A_norm, k=0,q=10, quantile.prob = 0.001, use.mkl = F, mkl.seed=-1) {
  # Computes the k-rank approximation to A_norm and adjusts it according to the
  # error distribution learned from the negative values.
  #
  # Args:
  #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
  #   k : the rank of the rank-k approximation. Set to 0 for automated choice of k.
  #   q : the number of additional power iterations in randomized SVD
  #   use.mkl : Use the Intel MKL based implementation of SVD. Needs to be
  #             installed from https://github.com/KlugerLab/rpca-mkl.
  #   mkl.seed : Only relevant if use.mkl=T. Set the seed for the random
  #   generator for the Intel MKL implementation of SVD. Any number <0 will
  #   use the current timestamp. If use.mkl=F, set the seed using
  #   set.seed() function as usual.
  #
  # Returns:
  #   A list with three items
  #       1) The rank k approximation of A_norm.
  #       2) The rank k approximation of A_norm, adaptively thresholded
  #       3) The rank k approximation of A_norm, adaptively thresholded and
  #       with the first two moments of the non-zero values matched to the
  #       first two moments of the non-zeros of A_norm. This is the completed
  #       matrix most people will want to work with
  # Example:
  #     result.completed <- adjusted_svd(A_norm,15)
  #     A_norm_rank15 <- result.completed[[1]]     # The low rank approximation for reference purposes...not suggested for matrix completion
  #     A_norm_rank15_cor <- result.completed[[3]] # The actual adjusted, completed matrix
  
  cat(sprintf("Read matrix with %d cells and %d genes\n", nrow(A_norm), ncol(A_norm)))
  if (!(any(grepl("matrix", class(A_norm), ignore.case = TRUE)))) { 
    stop(sprintf("A_norm is of class %s, but it should be of class matrix. Did you forget to run as.matrix()?",class(A_norm)))
  }
  
  if (k ==0 ) {
    k_choice <- choose_k(A_norm)
    k <-  k_choice$k
    cat(sprintf("Chose k=%d\n",k))
  }
  
  cat("Getting nonzeros\n")
  originally_nonzero <- A_norm >0 
  
  cat("Randomized SVD\n")
  if (!use.mkl) {
    fastDecomp_noc <- randomized.svd(A_norm,k,q=q)
  }else {
    fastDecomp_noc <- randomized.svd(A_norm,k,q=q, method='rsvd-mkl', mkl.seed=mkl.seed)
  }
# 重建矩阵# 重建矩阵# 重建矩阵# 重建矩阵# 重建矩阵# 重建矩阵
  A_norm_rank_k <- fastDecomp_noc$u[,1:k]%*%diag(fastDecomp_noc$d[1:k])%*% t(fastDecomp_noc$v[,1:k])    
  
  
  cat(sprintf("Find the %f quantile of each gene\n", quantile.prob))
  #A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,min))
  
#取占比0.1%的数值的绝对值，小于该绝对值的数值全部转变为0
  A_norm_rank_k_mins <- abs(apply(A_norm_rank_k,2,FUN=function(x) quantile(x,quantile.prob)))
  cat("Sweep\n")
  A_norm_rank_k_cor <- replace(A_norm_rank_k, A_norm_rank_k <= A_norm_rank_k_mins[col(A_norm_rank_k)], 0)
  
#scale 数据#scale 数据#scale 数据#scale 数据#scale 数据#scale 数据#scale 数据
  sd_nonzero <- function(x) sd(x[!x == 0])
  sigma_1 <- apply(A_norm_rank_k_cor, 2, sd_nonzero)
  sigma_2 <- apply(A_norm, 2, sd_nonzero)
  mu_1 <- colSums(A_norm_rank_k_cor)/colSums(!!A_norm_rank_k_cor)
  mu_2 <- colSums(A_norm)/colSums(!!A_norm)
  
  toscale <- !is.na(sigma_1) & !is.na(sigma_2) & !(sigma_1 == 0 & sigma_2 == 0) & !(sigma_1 == 0)
  
  cat(sprintf("Scaling all except for %d columns\n", sum(!toscale)))
  
  sigma_1_2 <- sigma_2/sigma_1
  toadd  <- -1*mu_1*sigma_2/sigma_1 + mu_2
  
  A_norm_rank_k_temp <- A_norm_rank_k_cor[,toscale]
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, sigma_1_2[toscale],FUN = "*")
  A_norm_rank_k_temp <- sweep(A_norm_rank_k_temp,2, toadd[toscale],FUN = "+")
  
  A_norm_rank_k_cor_sc <- A_norm_rank_k_cor
  A_norm_rank_k_cor_sc[,toscale] <- A_norm_rank_k_temp
  A_norm_rank_k_cor_sc[A_norm_rank_k_cor==0] = 0
  
  lt0 <- A_norm_rank_k_cor_sc  <0
  A_norm_rank_k_cor_sc[lt0] <- 0 
  cat(sprintf("%.2f%% of the values became negative in the scaling process and were set to zero\n", 100*sum(lt0)/(nrow(A_norm)*ncol(A_norm))))
  
  A_norm_rank_k_cor_sc[originally_nonzero & A_norm_rank_k_cor_sc ==0] <- A_norm[originally_nonzero & A_norm_rank_k_cor_sc ==0]
  
  colnames(A_norm_rank_k_cor) <- colnames(A_norm)
  colnames(A_norm_rank_k_cor_sc) <- colnames(A_norm)
  colnames(A_norm_rank_k) <- colnames(A_norm)
  
  original_nz <- sum(A_norm>0)/(nrow(A_norm)*ncol(A_norm))
  completed_nz <- sum(A_norm_rank_k_cor_sc>0)/(nrow(A_norm)*ncol(A_norm))
  cat(sprintf("The matrix went from %.2f%% nonzero to %.2f%% nonzero\n", 100*original_nz, 100*completed_nz))
  
  list(A_norm_rank_k=A_norm_rank_k,A_norm_rank_k_cor =A_norm_rank_k_cor, A_norm_rank_k_cor_sc=A_norm_rank_k_cor_sc)
}
```



