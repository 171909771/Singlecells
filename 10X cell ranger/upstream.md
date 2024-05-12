### 下载cellranger和index
- https://www.10xgenomics.com/support/software/cell-ranger/downloads
### 运行cellranger count
```
cellranger count --id=ctrl \
                 --transcriptome=/home/chan87/index/10xgenomicsindex/refdata-gex-mm10-2020-A \
                 --fastqs=/home/chan87/pengli2/rawdata \
                 --sample=CTRL_D3 \
                 --localcores=30 \
                 --localmem=256
```