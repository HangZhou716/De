# Developmental correspondence of juvenile stages across insects

Hang Zhou1*, Runguo Shu1†, Chaowei Zhang1, Yiqi Xiao1, Dong Jing1, Jiejing Tang1, Zixiong Cao2, Xi Chen1, Yang Mei1, Fei Li1*
1Key Laboratory of Biology of Crop Pathogens and Insects of Zhejiang Province, Institute of Insect Sciences, Zhejiang University, Hangzhou 310058, Zhejiang Province, China. 
2Object Research Systems (ORS) Inc. 
†These authors contributed equally: Hang Zhou, Runguo Shu. 
*Correspondence author: lifei18@zju.edu.cn, zhouhang716@zju.edu.cn.

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
#### DEG analysis
```r
library('DESeq2')
setwd('C:/Users/xmm/Desktop/fsdownload')
countData <- read.table(file='Haxy_counts.txt',header = T, row.names = 1) 
control <- 7 # count value started from colume 7
treatment <- 10

L2_vs_L1_countData <- all_data[, c(control:(control+2), treatment:(treatment+2))]
countData <- L2_vs_L1_countData
condition <- factor(c('L1','L1','L1','L2','L2','L2')) #set group name
countData <- round(as.matrix(countData)) 
coldata <- data.frame(row.names = colnames(countData),condition)

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition)
dds$condition <- relevel(dds$condition, ref = 'N5d3')
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
table(res$pvalue <0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by='row.names',sort=FALSE)
write.csv(resdata,file='Haxy_L2_vs_L1_DEGs.csv', row.names = FALSE)
```
#### variance-stabilised counts
```r
library(DESeq2)
library(openxlsx)

all_data <- read.table("Haxy_counts.txt",header = T, row.names = 1)
coldata <- data.frame(condition=factor(rep(c("L1","L2","L3","L4","PP1","PP2","PP3","PP4","P1","P2","P3"), each=3)))

dds <- DESeqDataSetFromMatrix(all_data, colData = coldata, design= ~ condition)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

vst_matrix <- assay(vsd)
vst_df <- as.data.frame(vst_matrix)
vst_df$Gene <- rownames(vst_df)
vst_df <- vst_df[, c(ncol(vst_df), 1:(ncol(vst_df) - 1))]

wb <- createWorkbook()
addWorksheet(wb, "VST Data")
writeData(wb, sheet = "VST Data", vst_df)

saveWorkbook(wb, "Haxy_VST_normalized_data.xlsx", overwrite = TRUE)

```
#### time-specific index
```python
import pandas as pd
import numpy as np
import openpyxl

def calculate_tau(row):
    x = row.values
    x_hat = x / np.max(x)
    tau = (np.sum(1 - x_hat)) / (len(x) - 1)
    return tau

# Read the expression matrix from an xlsx file
# Assuming your file is named 'expression_matrix.xlsx'
df = pd.read_excel('Haxy_VST_normalized_data.xlsx', index_col=0, engine='openpyxl')

# Calculate Tau index for each gene
tau_values = df.apply(calculate_tau, axis=1)

# Create a new DataFrame containing gene names and corresponding Tau indices
result = pd.DataFrame({'Gene': df.index, 'Tau': tau_values})

# Sort by Tau index in descending order
result = result.sort_values('Tau', ascending=False)

# Save the results to an xlsx file
result.to_excel('tau_index_results.xlsx', index=False, engine='openpyxl')

print("Tau index calculation completed. Results have been saved to 'tau_index_results.xlsx'.")

```

## The code in "Code S1.md" is primarily used for standardization of expression levels and clustering analysis. The vs-counts of different insects is available in Data S4.
