# 上传到linux处理
## 如果我们count的基因是基因名格式，需要添加参数--counts-data=gene_name，如果行名为ensemble名称的话，可以不添加这个参数，使用默认值即可。
cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt      --counts-data=gene_name

## cellphonedb 自己的绘图
cellphonedb plot dot_plot 
cellphonedb plot heatmap_plot cellphonedb_meta.txt 
