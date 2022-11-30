1. 计算矩阵中每个基因的均值
2. 将每个基因的均值+一个极小的随机的正态分布数值（rnorm(n = length(data.avg))/1e30）
3. 将上面的均值按照大小分为24个bins
```r
cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = 24, labels = FALSE, right = FALSE)
# labels =F： 不以区间展示，以分区数字展示
# right = FALSE：右面为开区间
```
4. 将每个细胞中每个兴趣基因所在的bin随机取100个基因来平均作为背景基因
5. 将每个细胞中每个兴趣基因的值减去对应的背景基因值：相当于得到每个基因在每个细胞中的另类Z-score
6. 将每个细胞中的兴趣基因值均值化：最终结果
