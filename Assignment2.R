# load packages needed
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(DESeq2)
library(ggplot2)
library(openxlsx)
library(knitr)
library(kableExtra)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(pathview)
library(pheatmap)
library(RColorBrewer)
library(glmnet)
library(survival)
library(survminer)
library(caret)

# 1. set my own directory
path = "/Users/afra/Desktop/Bio - R Assignment/Assignment 2"
setwd(path)

# 2. untar folder
folder_name = "brca_tcga_pan_can_atlas_2018.tar.gz"
folder = paste(path, folder_name, sep = "/")
untar(folder)
# go to new path
new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
setwd(new_dir)

# 3. read the RNA-seq file, using read delim
data_RNAseq <- read.delim("data_mrna_seq_v2_rsem.txt")

# 4. read the Patient Data file and skip 4 rows of column descriptions
data_patient <- read.delim("data_clinical_patient.txt", skip = 4, header = TRUE)

# 5. read the Copy Number Aberrations Data
data_cna <- read.delim("data_cna.txt")

# 6. match the RNASeq patient ids with the CNA ids and the Patient Data ids
for (i in 3:dim(data_RNAseq)[2]) { 
  # get the column name of the i-th column
  pat_barcode <- colnames(data_RNAseq)[i]
  # extract the first 12 characters of the column name
  pat_barcode = substr(pat_barcode, 1, 12)
  # replace "." with "-" to match the TCGA barcode
  pat_barcode = gsub("\\.", "-", pat_barcode)
  # update the column name with the cleaned barcode
  colnames(data_RNAseq)[i] <- pat_barcode  
}
for (i in 3:dim(data_cna)[2]){
  pat_barcode <- colnames(data_cna)[i] 
  pat_barcode <- substr(pat_barcode, 1, 12)
  pat_barcode <- gsub("\\.", "-",pat_barcode)
  colnames(data_cna)[i] <- pat_barcode
}
# match RNA-seq patient ids with the CNA ids
rna_cna_id <- which(colnames(data_RNAseq) %in% colnames(data_cna))
# filter patient data that are present in both RNA-seq and CNA
data_RNAseq <- data_RNAseq[, rna_cna_id]
# remove duplicated data
data_RNAseq <- data_RNAseq[!duplicated(data_RNAseq[, 1]), ]
# filter rows with non-empty Hugo_Symbol data
data_RNAseq <- data_RNAseq[data_RNAseq$Hugo_Symbol != "", ]

# 7. create metadata using the CNA level of ERBB2+
erbb2_cna <- data_cna |>
  filter(Hugo_Symbol == "ERBB2") |> # extract CNA data of ERBB2
  pivot_longer(cols = colnames(data_cna)[3]:colnames(data_cna)[dim(data_cna)[2]]
               , names_to = "PATIENT_ID", values_to = "ERBB2_CNA_LEVEL")
erbb2_cna$Amplified_Level <- ifelse(erbb2_cna$ERBB2_CNA_LEVEL > 0, "Amplified", "Not Amplified")  # label
# cols 1 and 2 are gene names
assay <- tibble(data_RNAseq[,-2])
colnames(assay)[1] <- "Gene" 
assay <- assay |>
  column_to_rownames(var = 'Gene')       
# build metadata
metadata <- matrix(0, dim(assay)[2], 2) 
pat_ids <- erbb2_cna$PATIENT_ID
col_cna <- which(colnames(erbb2_cna) == "ERBB2_CNA_LEVEL")

for (i in 1:dim(assay)[2]){
  idx = which(colnames(assay)[i] == pat_ids)
  metadata[i,1] = colnames(assay)[i]
  metadata[i,2] = 1*(as.numeric(erbb2_cna[idx, 4]) > 0)
}
metadata[is.na(metadata)] =0
colnames(metadata) <- c("PATIENT_ID","Amplified_Level")
metadata <- metadata |>
  data.frame() |>  # convert to data frame
  column_to_rownames(var = 'PATIENT_ID')
all(colnames(data_RNAseq[, -c(1, 2)]) %in% rownames(metadata))
metadata$Amplified_Level <- as.factor(ifelse(metadata$Amplified_Level > 0, "Amplified", "Not Amplified"))
table(metadata) # check number of ERBB2 CND levels -> (Amplified 328)

# 8. normalize data using DESeq2
assay[is.na(assay)] = 0  # impute with zeros the NA
assay[assay<0] = 0 
# create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(assay), 
  colData = metadata,       
  design = ~ Amplified_Level              
)
dds <- DESeq(dds)  # normalizes
resultsNames(dds)  # lists the coefficient

# 9. Obtain Differentially Expressed Genes.
dds$Amplified_Level <- factor(dds$Amplified_Level, levels = c("Amplified", "Not Amplified"))  # set factor level
res <- results(dds, contrast = c("Amplified_Level", "Amplified", "Not Amplified"), alpha = 0.05)
summary(res)  # DEG result summary

# (9.1) log fold change shrinkage for MA Plot
res_shrunken <- lfcShrink(dds, res = res, contrast = c("Amplified_Level", "Amplified", "Not Amplified"), type = "normal")
# unshrunken MA Plot
plotMA(res, ylim = c(-3, 3), xlab = "Mean of Normalized Counts", ylab = "Log Fold Change"
       , colNonSig = "grey", colSig = "steelblue2") 
abline(h=c(-1,1), col="grey", lty=2)
title(main = "(a)  MA Plot: Unshrunken Results", font.main = 1, adj = 0)

# shrunken MA Plot
plotMA(res_shrunken, ylim = c(-3, 3), xlab = "Mean of Normalized Counts", ylab = "Log Fold Change"
       , colNonSig = "grey", colSig = "tomato3") 
abline(h=c(-1,1), col="grey", lty=2)
title(main = "(b)  MA Plot: Shrunken Results", font.main = 1, adj = 0)

# (9.2) top 10 Differentially Expressed Genes Ranked by Fold Change
# set a significance threshold
padj_cutoff <- 0.05
lfc_cutoff <- 0  # changed from 0.58 -> 0.2 -> 0 
# top 10 most differentially expressed
significant <- which(res_shrunken$padj < padj_cutoff & abs(res_shrunken$log2FoldChange) >= lfc_cutoff)
significant_genes <- res_shrunken[significant, ] |>
  data.frame() |>
  rownames_to_column(var = "Gene") |>
  as_tibble() |>
  arrange(desc(abs(log2FoldChange)))
top10_genes <- significant_genes[1:10, ]
# table of top 10
top10_genes <- top10_genes |>
  mutate(
    `P-value (e-format)` = formatC(pvalue, format = "e", digits = 3),
    `Adj P-value (e-format)` = formatC(padj, format = "e", digits = 3)
  )
top10_table <- top10_genes |>
  dplyr::select(Gene, baseMean, log2FoldChange, lfcSE, `P-value (e-format)`, `Adj P-value (e-format)`) |>
  kbl(
    caption = "(a)  Top 10 Differentially Expressed Genes Ranked by Fold Change",
    col.names = c("Gene", "Base Mean", "Log2 Fold Change", "LFC SE", "P-value", "Adj P-value"),
    align = "c"
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))
top10_table

# extract normalised counts
normalized_counts <- counts(dds, normalized = TRUE) |>
  as.data.frame() |>
  rownames_to_column(var = "Gene")

top10_genes_data <- normalized_counts |>
  filter(Gene %in% top10_genes$Gene) |>
  pivot_longer(-Gene, names_to = "Sample", values_to = "Normalized_Count") |>
  mutate(Amplified_Level = metadata$Amplified_Level[match(Sample, rownames(metadata))])

# boxplot of top 10
ggplot(top10_genes_data) +
  geom_boxplot(aes(x = Gene, y = Normalized_Count, color = Amplified_Level)) + 
  scale_y_log10() +
  scale_fill_manual(values = c("steelblue2", "tomato3")) +
  labs(
    title = "(b)  Boxplot: Top 10 Differentially Expressed Gene",
    subtitle = "ERBB2 Amplified vs Not Amplified",
    x = "Genes",
    y = "log10 Normalized Counts",
    fill = "Amplified Level"
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 10, face = "italic"),
    legend.position = "top"
  )

# 10. perform a Pathway Enrichment Analysis
DE_over <- significant_genes %>% filter(log2FoldChange > 0)
DE_under <- significant_genes %>% filter(log2FoldChange < 0)

# convert Entrez ID
entrez_over <- data_RNAseq %>%
  filter(Hugo_Symbol %in% DE_over$Gene) %>%
  pull(Entrez_Gene_Id)

entrez_under <- data_RNAseq %>%
  filter(Hugo_Symbol %in% DE_under$Gene) %>%
  pull(Entrez_Gene_Id)

# perform Kegg pathway enrichment
kegg_results_over <- enrichKEGG(
  gene          = entrez_over,
  organism      = "hsa", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

kegg_results_under <- enrichKEGG(
  gene          = entrez_under,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
# check results
head(kegg_results_over)
head(kegg_results_under)
# plot results
dotplot(kegg_results_over, showCategory = 10) + ggtitle("(a)  Kegg Pathway Enrichment Over Expressed")
dotplot(kegg_results_under, showCategory= 10) + ggtitle("(b)  Kegg Pathway Enrichment Under Expressed") 

# 11. get the variance stabilised transformed expression values
vsd = vst(dds, blind = FALSE)

# 12. with the vst values obtain a PCA plot
pcaData <- plotPCA(vsd, intgroup = "Amplified_Level", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Amplified_Level)) +
  geom_point(size = 2, alpha=0.7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(
    title = "Principal Component Analysis (PCA) Plot",
    subtitle = "ERBB2 Amplified vs Not Amplified") +
  theme_light() +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 10, face = "italic"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "top"
  )

# 13. with the vst values obtain a heatmap
sig_genes <- significant_genes$Gene
vst_data <- assay(vsd) 
vst_sig_genes <- vst_data[rownames(vst_data) %in% sig_genes, ] 

#distance_matrix <- dist(vst_sig_genes, method = "euclidean")
#clustering_result <- hclust(distance_matrix, method = "ward.D2")
#plot(clustering_result, main = "Hierarchical Clustering of Significant Genes", xlab = "", sub = "")

#heatmap_colors <- colorRampPalette(c("steelblue", "white", "tomato3"))(50)
heatmap_colors <- brewer.pal(11,"RdYlBu")
breaks <- seq(-1, 1, length.out = length(heatmap_colors) + 1)
# heatmap of top 10 genes
top10_sig_genes <- significant_genes$Gene[1:10]
vst_top10_sig_genes <- vst_sig_genes[rownames(vst_sig_genes) %in% top10_sig_genes, ]
top10_correlation <- cor(t(vst_top10_sig_genes), method = "pearson")

pheatmap(
  top10_correlation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = breaks,
  scale = "row", 
  color = heatmap_colors,
  main = "(a)  Heatmap: Top 10 Differentially Expressed Genes"
)

# heatmap of significant genes
top100_sig_genes <- significant_genes$Gene[1:100]
vst_top100_sig_genes <- vst_sig_genes[rownames(vst_sig_genes) %in% top100_sig_genes, ]
top100_correlation <- cor(t(vst_top100_sig_genes), method = "pearson")
pheatmap(
  top100_correlation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  breaks = breaks,
  scale = "row", 
  color = heatmap_colors,
  main = "(b)  Heatmap: Significant Differentially Expressed Genes"
)

# 14. with the vst values of the DE genes generate an overall survival model
common_samples <- intersect(data_patient$PATIENT_ID, colnames(vst_sig_genes))
data_patient <- data_patient |> filter(PATIENT_ID %in% common_samples)
vst_sig_genes <- vst_sig_genes[, common_samples, drop = FALSE]
# prepare survival data
clinical_data <- data_patient |>
  dplyr::select(PATIENT_ID, OS_MONTHS, OS_STATUS) |>
  mutate(
    OS_STATUS = as.numeric(substr(OS_STATUS, 1, 1))  # convert OS_STATUS to numeric
  )
clinical_data <- clinical_data |> filter(OS_MONTHS > 0)
common_samples <- clinical_data$PATIENT_ID
vst_sig_genes <- vst_sig_genes[, common_samples, drop = FALSE]
# align vst_sig_genes with clinical data
# vst_sig_genes <- vst_sig_genes[, match(clinical_data$PATIENT_ID, colnames(vst_sig_genes))]
# create survival object
surv_obj <- Surv(time = clinical_data$OS_MONTHS, event = clinical_data$OS_STATUS)
# transpose and scale gene expression data for glmnet
vst_survival_data <- t(vst_sig_genes)
vst_survival_data_scaled <- scale(vst_survival_data)

# fit Lasso Regularized Cox Regression
lasso_model <- glmnet(
  x = vst_survival_data_scaled,
  y = surv_obj,
  family = "cox",
  alpha = 1  # Lasso
)
# -erform cross-validation to find the optimal lambda
cv_lasso <- cv.glmnet(
  x = vst_survival_data_scaled,
  y = surv_obj,
  family = "cox",
  alpha = 1  # Lasso Regularization
)
# plot cross-validation results
#plot(cv_lasso)
# extract the best lambda value (minimizes cross-validated error)
best_lambda <- cv_lasso$lambda.min
cat("Best lambda:", best_lambda, "\n")  # lambda 0.04616271

# extract coefficients at the best lambda value
lasso_coefs <- coef(cv_lasso, s = best_lambda)
lasso_coefs_matrix <- as.matrix(lasso_coefs)
nonzero_indices <- which(lasso_coefs_matrix != 0)
important_genes <- rownames(lasso_coefs_matrix)[nonzero_indices]
coefs <- lasso_coefs_matrix[nonzero_indices]

#cat("Selected genes by Lasso:\n")
#print(important_genes)

vst_important_genes <- vst_sig_genes[important_genes, common_samples, drop = FALSE]
vst_important_genes_t <- t(vst_important_genes)
# calculate risk scores
risk_scores <- vst_important_genes_t %*% coefs 
# combine risk scores with clinical data
combined_data <- clinical_data |> 
  mutate(
    risk_score = as.numeric(risk_scores)
  )
combined_data <- combined_data |> 
  mutate(
    risk_group = ifelse(risk_score > median(risk_score), "High Risk", "Low Risk")
  )

# build survival cure
surv_fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ risk_group, data = combined_data)
# Plot survival curves
ggsurvplot(
  surv_fit,
  data = combined_data,
  pval = TRUE, 
  title = "(a)  Survival Curve Based on Lasso Cox Regression",
  legend.title = "Risk Group",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = c("steelblue2", "tomato3")
)

combined_data <- combined_data |> 
  mutate(Amplified_Level = metadata$Amplified_Level[match(PATIENT_ID, rownames(metadata))])
# build survival object
surv_obj_KM <- Surv(time = combined_data$OS_MONTHS, event = combined_data$OS_STATUS)
# fit Kaplan-Meier model
surv_KM_fit <- survfit(surv_obj_KM ~ Amplified_Level, data = combined_data)
# plot Kaplan-Meier survival curves
ggsurvplot(
  surv_KM_fit,
  data = combined_data,
  conf.int = TRUE,
  pval = TRUE,  
  title = "(b)  Kaplan-Meier Survival Curve by ERBB2 Amplified Level",
  legend.title = "Level",
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  palette = c("steelblue2", "tomato3"), 
  risk.table = TRUE, 
  risk.table.height = 0.3, 
  ggtheme = theme_light() 
)
