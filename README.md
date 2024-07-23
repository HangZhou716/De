## Developmental correspondence of juvenile stages across the locust, harlequin ladybird, and diamondback moth

Hang Zhou1*, Runguo Shu1†, Chaowei Zhang1, Yiqi Xiao1, Dong Jing1, Jiejing Tang1, Zixiong Cao2, Xi Chen1, Yang Mei1, Fei Li1*
1Key Laboratory of Biology of Crop Pathogens and Insects of Zhejiang Province, Institute of Insect Sciences, Zhejiang University, Hangzhou 310058, Zhejiang Province, China. 
2Object Research Systems (ORS) Inc. 
†These authors contributed equally: Hang Zhou, Runguo Shu. 
*Correspondence author: lifei18@zju.edu.cn, zhouhang716@zju.edu.cn.

Code version:1.0

The code in "Code S1.md" is primarily used for standardization of expression levels and clustering analysis. The expression TPM matrix of different insects is available in Data S4.

### Trim_galore
```shell
trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 25 -e 0.1 --paired $dir/cmp/01raw_data/$fq1 $dir/cmp/01raw_data/$fq2 --gzip -o cleandata

```
### Hisat2
#### Create index
```shell
extract_splice_sites.py genome.gtf > ss.txt
extract_exons.py genome.gtf > exon.txt
hisat2-build --ss ss.txt --exon exon.txt Haxy_genome.fa Haxy_index

```
#### Map to genome
```shell
for i in cleandata/*_R1.fq.gz   
filename=${i#*/}
do
hisat2 --new-summary -p 40 -x Haxy_index -1 ${i} -2 ${i%_1*}_R2.fq.gz 2>log/${filename%_R1*}.log | samtools view -b -S - | samtools sort -o hisat2_out/${filename}.sorted.bam -

```
#### Caculate the read counts
```shell
featureCounts -T 10 -t exon -g gene_id -a Haxy_genome.gtf -o Haxy_counts.txt `ls hisat2_out/`
```
### DESeq2
DEG analysis
```r
library('DESeq2')
setwd('C:/Users/xmm/Desktop/fsdownload')
countData <- read.table(file='Haxy_counts.txt',header = T, row.names = 1) 
control <- 7 # count value started from colume 7
treatment <- 10

L2_vs_L1_countData <- all_data[, c(control:(control+2), treatment:(treatment+2))]

condition <- factor(c('L1','L1','L1','L2','L2','L2')) #set group name
countData <- round(as.matrix(countData)) 
coldata <- data.frame(row.names = colnames(countData),condition)

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition)
dds$condition <- relevel(dds$condition, ref = 'D9CK') #设置对照组
dds <- dds[ rowSums(counts(dds)) > 1, ] #数据筛选，行的count数之和小于1的删除
dds <- DESeq(dds) #计算FoldChange
res <- results(dds)
table(res$pvalue <0.05) #统计pvalue<0.05的结果
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by='row.names',sort=FALSE)
write.csv(resdata,file='D9-L_vs_D9-CK.csv')
