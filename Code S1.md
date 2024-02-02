### Sample PCA
```r
library(Mfuzz)
library(readxl)
library(ggpubr)
library(ggthemes)
library(gmodels)
library(ggplot2)

rm(list=ls())
# read table
data <- read_excel("Data S4.xlsx")

# Split data into four groups: Dm, Ha, Px, Lm
groups <- list(
  Dm = data[, grepl("Dm-", colnames(data))],
  Ha = data[, grepl("Ha-", colnames(data))],
  Px = data[, grepl("Px-", colnames(data))],
  Lm = data[, grepl("Lm-", colnames(data))]
)

# Z-score standardization and filtering
standardized_groups <- lapply(groups, function(group) {
  eset <- new("ExpressionSet", exprs = as.matrix(group))
  eset <- filter.NA(eset, thres = 0.25)
  eset <- filter.std(eset, min.std = 0)
  return(standardise(eset))
})

# Identify genes common across all groups
common_genes <- Reduce(intersect, lapply(standardized_groups, rownames))

# Retain only these common genes
filtered_groups <- lapply(standardized_groups, function(eset) {
  eset[common_genes, ]
})

# Merge the data
merged_data <- do.call(cbind, lapply(filtered_groups, exprs))

index <- read.table("../../sample_PCA/pca_index.txt", header = T, row.names = 1)

# run PCA
pca.info <- fast.prcomp(merged_data)
pca.data <- data.frame(sample = rownames(pca.info$rotation),
                       Type = c(rep("D.melanogaster", 10),
                                rep("H.axyridis", 13),
                                rep("P.xylostella", 11),
                                rep("L.migratoria", 7)
                               ), pca.info$rotation)

pca.data$sample <- as.numeric(index["GROUP",])

set.seed(123) # random seed

variance_explained <- summary(pca.info)$importance[2, ] * 100

x_label <- paste("PC1 (", round(variance_explained[1], 2), "% variance explained)")
y_label <- paste("PC2 (", round(variance_explained[2], 2), "% variance explained)")

kmeans_result <- kmeans(pca.data[, c("PC1", "PC2")], centers = 3) # k=3
pca.data$cluster <- as.factor(kmeans_result$cluster) # add the clustering result

ggplot(pca.data, aes(x = PC1, y = PC2, color = cluster)) +
geom_point(size = 1) +
geom_text(aes(label = rownames(pca.data)), vjust = -1, size = 3) +  # add labels
scale_color_manual(values = c("#cb0000","#3398cb","#319700","#cb3398")) +
labs(x = x_label, y = y_label) +
theme_base() +
theme(legend.position = "right")


```
### UMAP clustering
```r
library(Mfuzz)
library(umap)
library(ggplot2)
library(readxl)

rm(list=ls())
# read table
data <- read_excel("Data S4.xlsx")

# Split data into four groups: Dm, Ha, Px, Lm
groups <- list(
  Dm = data[, grepl("Dm-", colnames(data))],
  Ha = data[, grepl("Ha-", colnames(data))],
  Px = data[, grepl("Px-", colnames(data))],
  Lm = data[, grepl("Lm-", colnames(data))]
)

# Mfuzz standardization and filtering
standardized_groups <- lapply(groups, function(group) {
  eset <- new("ExpressionSet", exprs = as.matrix(group))
  eset <- filter.NA(eset, thres = 0.25)
  eset <- filter.std(eset, min.std = 0)
  return(standardise(eset))
})

# Identify genes common across all groups
common_genes <- Reduce(intersect, lapply(standardized_groups, rownames))

# Retain only these common genes
filtered_groups <- lapply(standardized_groups, function(eset) {
  eset[common_genes, ]
})

# Merge the data
merged_data <- do.call(cbind, lapply(filtered_groups, exprs))

set.seed(2023)

# UMAP clustering analysis
umap_result <- umap(t(merged_data),n_components=2, n_neighbors = 8)

kmeans_result <- kmeans(umap_result$layout, centers=3)

# Create a dataframe containing UMAP coordinates and clustering results
umap_df <- as.data.frame(umap_result$layout)
umap_df$cluster <- as.factor(kmeans_result$cluster)
umap_df$sample <- rownames(umap_df)

# Plot using ggplot2 and add data labels
ggplot(umap_df, aes(x=V1, y=V2, label=sample, color=cluster)) +
  geom_point() +
  geom_text(aes(label=sample), check_overlap = FALSE, nudge_x = 0.1, nudge_y = 0.1, size=3) +
  scale_color_manual(values=c("#BE1A21", "#348FBE", "#359139")) + 
  theme_bw() +
  theme(panel.grid = element_blank())+
  labs(title="UMAP projection with K-means Clustering", x="UMAP 1", y="UMAP 2") +
  theme(legend.position= "none") +
  coord_fixed(ratio = 0.5)

```
### Hierarchical clustering
```r
library(Mfuzz)
library(vegan)
library(readxl)
library(dendextend)

rm(list=ls())
# read table
data <- read_excel("Data S4.xlsx")

# Split data into four groups: Dm, Ha, Px, Lm
groups <- list(
  Dm = data[, grepl("Dm-", colnames(data))],
  Ha = data[, grepl("Ha-", colnames(data))],
  Px = data[, grepl("Px-", colnames(data))],
  Lm = data[, grepl("Lm-", colnames(data))]
)

# Mfuzz z-score standardization and filtering
standardized_groups <- lapply(groups, function(group) {
  eset <- new("ExpressionSet", exprs = as.matrix(group))
  eset <- filter.NA(eset, thres = 0.25)
  eset <- filter.std(eset, min.std = 0)
  return(standardise(eset))
})

# Identify genes common across all groups
common_genes <- Reduce(intersect, lapply(standardized_groups, rownames))

# Retain only these common genes
filtered_groups <- lapply(standardized_groups, function(eset) {
  eset[common_genes, ]
})

# Merge the data
merged_data <- do.call(cbind, lapply(filtered_groups, exprs))

dist_matrix <- dist(t(merged_data))

hc_result <- hclust(dist_matrix, method = "complete")

plot(hc_result)

rect.hclust(hc_result, k =3)

```
