# Installations and setup:

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

library("DESeq2", "dplyr", "ggplot2")

setwd("/Volumes/Apple/rnaseq_AC_proper")
getwd()
# ======================================================================================================================
#======================================================================================================================

# Loading the raw counts file:
raw_counts_file <- read.delim("counts/counts.txt", comment.char="#")

# Cleaning the column names:
# Remove 'alignment/' prefix and '.sorted.bam' suffix:
colnames(raw_counts_file)[-c(1:6)] <- gsub("alignment\\.|\\.sorted\\.bam", "", colnames(raw_counts_file)[-c(1:6)])

# Exclude first 6 columns (annotation info from featureCounts)
count_data <- raw_counts_file[, -(1:6)]
gene_ids <- raw_counts_file$Geneid

# Merge technical replicates (adjust names as per your real column names)
merged_counts <- data.frame(
  Geneid = gene_ids,
  
  WT_1 = rowSums(count_data[, c("WT_1A", "WT_1B", "WT_1C")]),
  WT_2 = rowSums(count_data[, c("WT_2A", "WT_2B", "WT_2C")]),
  WT_3 = rowSums(count_data[, c("WT_3A", "WT_3B", "WT_3C")]),
  WT_4 = rowSums(count_data[, c("WT_4A", "WT_4B", "WT_4C")]),
  
  rrp6_1 = rowSums(count_data[, c("rrp6_1A", "rrp6_1B", "rrp6_1C")]),
  rrp6_2 = rowSums(count_data[, c("rrp6_2A", "rrp6_2B", "rrp6_2C")]),
  rrp6_3 = rowSums(count_data[, c("rrp6_3A", "rrp6_3B", "rrp6_3C")]),
  rrp6_4 = rowSums(count_data[, c("rrp6_4A", "rrp6_4B", "rrp6_4C")])
)

# ========================================================================================================================

# Prepping for DESeq2:
# Step 1. Deseq2 requires geneIDs as rownames. Not as a column:

count_matrix <- merged_counts[,-1]
rownames(count_matrix) <- merged_counts$Geneid

# Step 2. Creating the metadata table:
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = c("WT", "WT", "WT", "WT", "rrp6", "rrp6", "rrp6", "rrp6")
)

# Step 3. Setting factors and their levels:
sample_info$condition <- factor(sample_info$condition)
sample_info$condition <- factor(sample_info$condition, levels = c("WT", "rrp6"))

# Step 4. Creating the DESeq2 dataset:
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)

# Step 5. Running DESeq2 analysis:
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]

# Step 6. Saving the results:
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")

# Step 7. QC plots:
# MA Plot:
plotMA(res, ylim = c(-8, 8))

# PCA Plot:
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition", returnData = FALSE)

# Outlier spotted in both conditions. 1 of 4 replicates.

##================================================================================================================================

# Identifying outlier:
vsd_detect <- vsd
colData(vsd_detect)$name <- colnames(dds)
colData(vsd_detect)
plotPCA(vsd_detect, intgroup = "name")

# A more refined plot to identify the outlier:
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, show.legend = FALSE) +
  labs(title = "PCA of RNA-seq samples",
       x = paste0("PC1: ", round(100 * (attr(pcaData, "percentVar")[1])), "% variance"),
       y = paste0("PC2: ", round(100 * (attr(pcaData, "percentVar")[2])), "% variance")
  ) +
  theme_classic()

# Removing the outliers:
filt_count_matrix <- count_matrix[, !colnames(count_matrix) %in% c("WT_1", "rrp6_1")]

filt_sample_info <- sample_info
filt_sample_info$name <- rownames(sample_info)
filt_sample_info <- filt_sample_info |> dplyr::filter(!name %in% c("WT_1", "rrp6_1")) |> dplyr::select(!name)

# Re-running deseq2:
dds <- DESeqDataSetFromMatrix(countData = filt_count_matrix,
                              colData = filt_sample_info,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]

# Saving the results:
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")

# QC plots:
# MA Plot:
plotMA(res, ylim = c(-8, 8))

# PCA Plot:
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition", returnData = FALSE)

# Volcano Plot prepr:
res_df <- as.data.frame(resOrdered)
res_df$gene <- rownames(res_df)
res_df <- na.omit(res_df)  # Remove rows with NA padj or log2FoldChange
library("ggplot2")

# Add a column to classify points as up/down/non-significant
res_df$threshold <- "Not Sig"
res_df$threshold[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$threshold[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
res_df$threshold <- factor(res_df$threshold, levels = c("Up", "Not Sig", "Down"))

# Basic volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 FC",
       y = "-Log10 p-value")+
  theme(legend.title = element_blank())

# Labelling top genes:
# install.packages("ggrepel")
library(ggrepel)

top_genes <- res_df[res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 2, ]

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 10) +
  theme_classic() +
  labs(title = "Volcano Plot with Top Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme(legend.title = element_blank())

## ===============================================================================================================================
## ===============================================================================================================================

### Filtering and working with snRNAs:

# Filtering out only snR genes from counts matrix:
snR_counts <- filt_count_matrix[grep("snR", rownames(filt_count_matrix)),]

# Filtering out snR genes from deseq2 results table:
res_snR <- res_df[grep("snR", rownames(res_df)),]

write.csv(res_snR, file = "snR_folds.csv")

# Marking top 10 snR genes by FC:
top_10_snR <- res_snR |> 
  filter(padj < 0.05) |> 
  arrange(desc(log2FoldChange)) |>
  slice_head(n = 10)

# Plotting snR only:
custom_colors <- c(
  "Up" = "red",
  "Down" = "blue",
  "Not Sig" = "grey"
)

ggplot(res_snR, aes(x = log2FoldChange, y = -log10(pvalue), colour = threshold))+
  geom_point()+
  scale_color_manual(values = custom_colors,
                     breaks = intersect(c("Up", "Not Sig", "Down"), unique(res_df$threshold)),
                     name = "Expression change")+
  labs(title = "Volcano Plot snR genes",
       x = "Log2FC",
       y = "-Log10 p-value") +
  geom_text_repel(data = top_10_snR,
                  aes(label = gene),
                  show.legend = FALSE,
                  max.overlaps = Inf,
                  box.padding = 0.8,
                  point.padding = 0.5)

#=====================================================================================================================================

# Creating a heatmap of snR genes:

# Creating the genes list
snR_list <- grep("snR", rownames(snR_counts), value = TRUE)

# Extracting the snR counts from variance stabilized data (vsd):
snR_vsd <- assay(vsd)[snR_list, ]

# Z-score normalisation:
snr_scaled <- na.omit(t(scale(t(snR_vsd))))

# Prepping annotation column:

annotation_col <- as.data.frame(colData(vsd)[, "condition"])

library(pheatmap)

pheatmap(
  snr_scaled,
  # annotation_col = colnames(snr_scaled),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 7,
  scale = "row",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  angle_col = 45,
  legend_breaks = c(1.5, 0, -1.5)
)



