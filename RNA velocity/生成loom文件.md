- https://mp.weixin.qq.com/s/TJsdj7qesZGlvmEkNLAy0A
# 生成 loom文件
- https://velocyto.org/velocyto.py/tutorial/cli.html ####官方

rmsk=/home/chan87/hg38_rmsk.gtf       #UCSC table browser 下载
outdir=/home/chan87/out      # cellranger生成
cellranger=/home/chan87/index/GRCh38/genes/genes.gtf         #cellranger生成

## 先用samtools 在out文件夹中把bam文件转制
```shell
nohup samtools sort -@ 10  -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam &
```
## 正式程序运行
- http://velocyto.org/velocyto.py/tutorial/cli.html#run-smartseq2-run-on-smartseq2-samples
nohup velocyto run10x -@ 10 -m $rmsk $outdir $cellranger &
### 如果显示index error，可能是重复序列的问题，注意重新下载重复序列hg38_rmsk.gtf 
#### 自己制作重复序列，下面网址
- https://www.jianshu.com/p/8e3f5e96d690
