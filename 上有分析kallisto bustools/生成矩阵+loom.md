# 构建index
nohup kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno -n 8 \
Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
Mus_musculus.GRCm39.106.gtf.gz &

# 生成counts
## 注意index之间用逗号隔开不能再加空格
nohup kb count --h5ad -i index.idx_cdna,index.idx_intron.0,index.idx_intron.1,index.idx_intron.2,index.idx_intron.3,index.idx_intron.4,index.idx_intron.5,index.idx_intron.6 \
-g t2g.txt -x 10xv2 -o out \
-c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno --verbose --filter bustools -t 14  \
fastqs/SI-GA-A1_1/SI-GA-A1_1_S1_L001_R1_001.fastq.gz \
fastqs/SI-GA-A1_1/SI-GA-A1_1_S1_L001_R2_001.fastq.gz \
fastqs/SI-GA-A1_2/SI-GA-A1_2_S2_L001_R1_001.fastq.gz \
fastqs/SI-GA-A1_2/SI-GA-A1_2_S2_L001_R2_001.fastq.gz \
fastqs/SI-GA-A1_3/SI-GA-A1_3_S3_L001_R1_001.fastq.gz \
fastqs/SI-GA-A1_3/SI-GA-A1_3_S3_L001_R2_001.fastq.gz \
fastqs/SI-GA-A1_4/SI-GA-A1_4_S4_L001_R1_001.fastq.gz \
fastqs/SI-GA-A1_4/SI-GA-A1_4_S4_L001_R2_001.fastq.gz &
