# Final QC

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stats)
})

project_dir     <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
normalized_dir  <- file.path(project_dir, "results", "normalized_matrix")
qc_dir          <- file.path(project_dir, "results", "QC")

cat("Locating and reading sample QC metadata\n")

step02_candidates <- list.files(
  qc_dir, 
  pattern = "^sample_qc_.*1106\\.csv$", 
  full.names = TRUE
)

if (length(step02_candidates) == 0) {
  # Fallback to standard explicit filenames
  step02_candidates <- file.path(qc_dir, c("sample_qc_metrics_1106.csv", "sample_qc_summary_1106.csv"))
  step02_candidates <- step02_candidates[file.exists(step02_candidates)]
}

if (length(step02_candidates) == 0) {
  stop("CRITICAL ERROR: Step 02 QC summary CSV file not found in: ", qc_dir)
}

step02_file <- step02_candidates[1]
cat("  Found Step 02 QC file:", basename(step02_file), "\n")

step02_qc <- fread(step02_file)

# Assert column existence with explicit mapping to library_size
if (!"library_size" %in% colnames(step02_qc)) {
  if ("total_counts" %in% colnames(step02_qc)) {
    setnames(step02_qc, "total_counts", "library_size")
  } else {
    stop("CRITICAL ERROR: 'library_size' (or 'total_counts') column missing from Step 02 QC file.")
  }
}

stopifnot(all(c("sample_id", "library_size", "detected_genes", "zero_fraction", "low_depth_flag") %in% colnames(step02_qc)))

cat("  Successfully loaded Step 02 QC metrics for", nrow(step02_qc), "samples.\n")


#load normalized log2-cpm data
cat("\nLoading log2-CPM expression matrix\n")

log2_cpm <- readRDS(file.path(normalized_dir, "log2_cpm_1106.rds"))

n_genes   <- nrow(log2_cpm)
n_samples <- ncol(log2_cpm)

#Ensure sample ordering alignment between matrix and QC table
stopifnot(all(colnames(log2_cpm) %in% step02_qc$sample_id))
step02_qc <- step02_qc[match(colnames(log2_cpm), step02_qc$sample_id), ]
stopifnot(identical(colnames(log2_cpm), step02_qc$sample_id))

cat("  Loaded log2-CPM matrix:", format(n_genes, big.mark = ","), "genes x", 
    format(n_samples, big.mark = ","), "samples\n")

# 4. HIGHLY VARIABLE GENE (HVG) SELECTION
cat("\nExtracting Top 1,000 Highly Variable Genes (HVGs) via MAD\n")

gene_mad      <- apply(log2_cpm, 1, mad)
top_hvg_count <- 1000
hvg_genes     <- names(sort(gene_mad, decreasing = TRUE))[1:top_hvg_count]
hvg_matrix    <- log2_cpm[hvg_genes, ]

cat("  Selected top", top_hvg_count, "HVGs (MAD range:", 
    round(min(gene_mad[hvg_genes]), 3), "to", round(max(gene_mad[hvg_genes]), 3), ")\n")

#Principal component analysis (pca)

cat("\nPerforming PCA on Top 1,000 HVGs\n")

pca_res <- prcomp(t(hvg_matrix), center = TRUE, scale. = TRUE)

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pc1_pve <- round(var_explained[1] * 100, 2)
pc2_pve <- round(var_explained[2] * 100, 2)
pc3_pve <- round(var_explained[3] * 100, 2)

cat("  Variance Explained -> PC1:", pc1_pve, "% | PC2:", pc2_pve, "% | PC3:", pc3_pve, "%\n")


#Sample correlation & hierarchical clustering
cat("\nComputing inter-sample Pearson correlations and hierarchical clustering...\n")

#Calculate Pearson correlation matrix
sample_cor <- cor(hvg_matrix, method = "pearson")

#Remove self-correlation diagonal (r = 1.0) before calculating mean
sample_cor_no_diag <- sample_cor
diag(sample_cor_no_diag) <- NA
mean_sample_cor <- rowMeans(sample_cor_no_diag, na.rm = TRUE)

#Data-driven IQR threshold for candidate identification
q1_cor  <- quantile(mean_sample_cor, 0.25)
q3_cor  <- quantile(mean_sample_cor, 0.75)
iqr_cor <- IQR(mean_sample_cor)

low_cor_threshold <- q1_cor - (1.5 * iqr_cor)

low_cor_candidates <- names(mean_sample_cor[mean_sample_cor < low_cor_threshold])

cat("  Mean Pearson Correlation Stats:\n")
cat("  -> Min      :", round(min(mean_sample_cor), 4), "\n")
cat("  -> Median   :", round(median(mean_sample_cor), 4), "\n")
cat("  -> Max      :", round(max(mean_sample_cor), 4), "\n")
cat("  -> 1.5xIQR Lower Bound Threshold:", round(low_cor_threshold, 4), "\n")
cat("  -> Low Correlation Candidates Identified:", length(low_cor_candidates), "\n")

#Hierarchical Clustering (Average Linkage on Correlation Distance)
sample_dist <- as.dist(1 - sample_cor)
hc_res      <- hclust(sample_dist, method = "average")

#Integrated qc metrics dataframe & candidate identification

cat("\n Building integrated sample QC metadata dataframe...\n")

qc_df <- data.frame(
  sample_id                 = colnames(log2_cpm),
  library_size              = step02_qc$library_size,
  detected_genes            = step02_qc$detected_genes,
  zero_fraction             = step02_qc$zero_fraction,
  step02_flag               = step02_qc$low_depth_flag,
  qc_status_step02          = ifelse(step02_qc$low_depth_flag, "Step 02 Flagged (Low Depth)", "Not Flagged in Step 02"),
  PC1                       = pca_res$x[, 1],
  PC2                       = pca_res$x[, 2],
  PC3                       = pca_res$x[, 3],
  mean_inter_cor            = mean_sample_cor,
  low_correlation_candidate = mean_sample_cor < low_cor_threshold
)

#Identify overlap between Step 02 flags and low correlation candidates
overlap_candidates <- qc_df[qc_df$step02_flag & qc_df$low_correlation_candidate, "sample_id"]
cat("  Overlap (Step 02 Flagged AND Low Correlation Candidate):", length(overlap_candidates), "\n")

#Ddiagnostic visualizations
cat("\nGenerating diagnostic QC visualizations\n")

#Plot A: PCA Scatter Plot colored by Step 02 Status
p_pca <- ggplot(qc_df, aes(x = PC1, y = PC2, color = qc_status_step02)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("Not Flagged in Step 02" = "#2B5C8F", "Step 02 Flagged (Low Depth)" = "#E41A1C")) +
  labs(
    title    = "PCA Score Plot (Top 1,000 HVGs)",
    subtitle = paste0("TCGA-BRCA Cohort (N = ", n_samples, ")"),
    x        = paste0("PC1 (", pc1_pve, "% Variance)"),
    y        = paste0("PC2 (", pc2_pve, "% Variance)"),
    color    = "Step 02 QC Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(qc_dir, "08_pca_pc1_pc2.png"), plot = p_pca, width = 7, height = 5.5, dpi = 300)

# Plot B: Scree Plot (Top 10 PCs)
scree_df <- data.frame(PC = paste0("PC", 1:10), PVE = var_explained[1:10] * 100)
scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)

p_scree <- ggplot(scree_df, aes(x = PC, y = PVE)) +
  geom_col(fill = "#2B5C8F", width = 0.6) +
  geom_text(aes(label = paste0(round(PVE, 1), "%")), vjust = -0.5, size = 3.2) +
  labs(
    title = "PCA Scree Plot (Top 10 Principal Components)",
    x     = "Principal Component",
    y     = "Proportion of Variance Explained (%)"
  ) +
  ylim(0, max(scree_df$PVE) * 1.15) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(qc_dir, "09_pca_scree_plot.png"), plot = p_scree, width = 7, height = 4.5, dpi = 300)

#Plot C: Mean Pearson Correlation Distribution
p_cor <- ggplot(qc_df, aes(x = mean_inter_cor)) +
  geom_histogram(bins = 40, fill = "#2B5C8F", color = "white", alpha = 0.8) +
  geom_vline(xintercept = low_cor_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    title    = "Distribution of Mean Inter-Sample Pearson Correlations",
    subtitle = paste0("Red dashed line = IQR threshold (r = ", round(low_cor_threshold, 4), ")"),
    x        = "Mean Inter-Sample Pearson Correlation (r)",
    y        = "Sample Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(qc_dir, "10_sample_correlation_distribution.png"), plot = p_cor, width = 7, height = 4.5, dpi = 300)

#Plot D: Hierarchical Clustering Dendrogram (Clean Structural View)
png(file.path(qc_dir, "11_sample_hierarchical_clustering_dendrogram.png"), width = 2400, height = 1200, res = 300)
par(mar = c(3, 4, 3, 1))
plot(
  hc_res, 
  labels = FALSE, 
  main   = paste0("Sample Hierarchical Clustering Dendrogram (N = ", n_samples, ")"),
  sub    = "Average Linkage on Correlation Distance (1 - r)",
  xlab   = "", 
  ylab   = "Height (1 - r)"
)
dev.off()

cat("\nExporting integrated QC metric tables and audit summaries\n")

# Save integrated metrics table
write.csv(qc_df, file = file.path(qc_dir, "sample_level_step05_qc_metrics.csv"), row.names = FALSE)

# Save R data object
pca_output <- list(
  pca_object             = pca_res,
  qc_dataframe           = qc_df,
  var_explained          = var_explained,
  hvg_genes              = hvg_genes,
  sample_cor             = sample_cor,
  hierarchical_clustering = hc_res,
  low_cor_threshold      = low_cor_threshold
)

saveRDS(pca_output, file = file.path(qc_dir, "pca_qc_results.rds"))

summary_file <- file.path(qc_dir, "processing_summary_05.txt")
sink(summary_file)

cat("Method:                     Median Absolute Deviation (MAD)\n")
cat("Selected Features:          Top 1,000 HVGs\n")

cat("\nPCA variance explained (pve) (PVE)\n")
cat("----------------------------\n")
cat("PC1:                        ", pc1_pve, "%\n", sep = "")
cat("PC2:                        ", pc2_pve, "%\n", sep = "")
cat("PC3:                        ", pc3_pve, "%\n", sep = "")
cat("Cumulative (PC1-PC3):       ", round(pc1_pve + pc2_pve + pc3_pve, 2), "%\n", sep = "")

cat("\nInter-sample pearson correlation\n")
cat("--------------------------------\n")
cat("Min Mean Correlation:       ", round(min(mean_sample_cor), 4), "\n")
cat("Median Mean Correlation:    ", round(median(mean_sample_cor), 4), "\n")
cat("Max Mean Correlation:       ", round(max(mean_sample_cor), 4), "\n")
cat("1.5xIQR Lower Bound Threshold: ", round(low_cor_threshold, 4), "\n")
cat("Low Correlation Candidates: ", length(low_cor_candidates), "\n")

cat("\nMulti-evidence overlap audit\n")
cat("----------------------------\n")
cat("Step 02 Low Depth Flagged:  ", sum(qc_df$step02_flag), "\n")
cat("Low Correlation Candidate:  ", sum(qc_df$low_correlation_candidate), "\n")
cat("Dual-Flagged Overlap Count: ", length(overlap_candidates), "\n")

cat("\nOutput files\n")
cat("------------\n")
cat("results/QC/08_pca_pc1_pc2.png\n")
cat("results/QC/09_pca_scree_plot.png\n")
cat("results/QC/10_sample_correlation_distribution.png\n")
cat("results/QC/11_sample_hierarchical_clustering_dendrogram.png\n")
cat("results/QC/sample_level_step05_qc_metrics.csv\n")
cat("results/QC/pca_qc_results.rds\n")
cat("results/QC/processing_summary_05.txt\n")

sink()

cat("QC plots and integrated metrics saved to:\n")
cat(qc_dir, "\n")
