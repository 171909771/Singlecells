- https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247484355&idx=1&sn=7860fe0c46073a55d2d3700822c3103b&scene=21#wechat_redirect

# 构建注释,包含了想要的基因
cellranger mkgtf Homo_sapiens.GRCh38.106.gtf Homo_sapiens.GRCh38.106.filtered.gtf \
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
