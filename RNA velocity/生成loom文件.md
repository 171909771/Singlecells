- https://mp.weixin.qq.com/s/TJsdj7qesZGlvmEkNLAy0A
# 生成 loom文件

rmsk=/home/chan87/hg38_rmsk.gtf       #UCSC table browser 下载
outdir=/home/chan87/out      # cellranger生成
cellranger=/home/chan87/index/GRCh38/genes/genes.gtf         #cellranger生成

## 先用samtools 在out文件夹中把bam文件转制
nohup samtools sort -@ 10  -t CB -O BAM -o cellsorted_possorted_genome_bam.bam possorted_genome_bam.bam &
## 正式程序运行
- http://velocyto.org/velocyto.py/tutorial/cli.html#run-smartseq2-run-on-smartseq2-samples
nohup velocyto run10x -@ 10 -m $rmsk $outdir $cellranger &
