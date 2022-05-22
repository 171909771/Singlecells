- https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484355&idx=1&sn=7860fe0c46073a55d2d3700822c3103b&scene=21#wechat_redirect

# 构建注释,包含了想要的基因
cellranger mkgtf Homo_sapiens.GRCh38.106.gtf Homo_sapiens.GRCh38.106.filtered.gtf \  ###下次只要编码蛋白基因
                --attribute=gene_biotype:protein_coding \
                --attribute=gene_biotype:lincRNA \
                --attribute=gene_biotype:antisense \
                --attribute=gene_biotype:IG_LV_gene \
                --attribute=gene_biotype:IG_V_gene \
                --attribute=gene_biotype:IG_V_pseudogene \
                --attribute=gene_biotype:IG_D_gene \
                --attribute=gene_biotype:IG_J_gene \
                --attribute=gene_biotype:IG_J_pseudogene \
                --attribute=gene_biotype:IG_C_gene \
                --attribute=gene_biotype:IG_C_pseudogene \
                --attribute=gene_biotype:TR_V_gene \
                --attribute=gene_biotype:TR_V_pseudogene \
                --attribute=gene_biotype:TR_D_gene \
                --attribute=gene_biotype:TR_J_gene \
                --attribute=gene_biotype:TR_J_pseudogene \
                --attribute=gene_biotype:TR_C_gene  


# 构建索引
cellranger mkref --genome=GRCh38 \
                --nthreads=12 \
                --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                --genes=Homo_sapiens.GRCh38.106.filtered.gtf &

# 改名字
mv SRR7692286_1.fastq.gz SRR7692286_S1_L001_I1_001.fastq.gz 
mv SRR7692286_2.fastq.gz SRR7692286_S1_L001_R1_001.fastq.gz 
mv SRR7692286_3.fastq.gz SRR7692286_S1_L001_R2_001.fastq.gz

# 表达矩阵
nohup cellranger count --id=out \
                  --transcriptome=/home/chan87/index/GRCh38 \
                  --fastqs=/home/chan87/test1/SRR7692286 \
                  --sample=SRR7692286 \
                  --nosecondary \
                  --localcores=12 \
                  --localmem=28 &

# 如果是4个SRR为一个样本时用下面的代码
- https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.html#gsc.tab=0  #参考网页
## cellranger

### 改名
```shell
ls SRR*|while read id ; do echo ${id%_*};done|uniq>txt1   # 取值每一个SRR值
cat txt1|while read id;do mv ${id}_1.fastq.gz ${id}_S1_L00${a}_R1_001.fastq.gz;mv ${id}_2.fastq.gz ${id}_S1_L00${a}_R2_001.fastq.gz; a=`expr $a + 1` ;done    #改名
```
### 启动cellranger count 生成矩阵文件
```shell
nohup cellranger count --id=out  \
--transcriptome=/home/chan87/index/mm39 \
--fastqs=/home/chan87/bio.test/SRP320164/  \
--sample=SRR14570592,SRR14570593,SRR14570594,SRR14570595 \
--nosecondary  \
--localcores=12  \
--localmem=28 &
```
#### 如果是技术重复样本，用并行，先用分别建立文件夹，用find >出绝对路径
```shell
a=1
cat /home/chan87/tmp/123|while read id; do  ( nohup cellranger count \
--id=out$a  \
--transcriptome=/home/chan87/index/mm39  \
--fastqs=${id}/ \
--sample=${id##*/}  \
--nosecondary  \
--localcores=12  \
--localmem=28  &); a=`expr $a + 1`;done
```
### 启动cellranger aggr 合并矩阵文件
- https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
- https://www.jianshu.com/p/ca726a8979d7    #中文
#### csv文件建立
![image](https://user-images.githubusercontent.com/41554601/168627469-90b2067c-f9d6-43a3-837e-56179ac38ea1.png)

cellranger aggr --id=merge --csv=merge.csv
[merge.csv](https://github.com/171909771/Singlecells/files/8749018/merge.csv)

