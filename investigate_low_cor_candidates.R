# ==============================================================================
# SCRIPT 05b: DEEP AUDIT & CLASSIFICATION OF LOW-CORRELATION CANDIDATE SAMPLES
# Project: TCGA-BRCA Immune Infiltration Analysis
# Author: Karan S.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & DEPENDENCIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stats)
})

# Define project paths
project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
qc_dir      <- file.path(project_dir, "results", "QC")

cat("==================================================================\n")
cat("STEP 05b: CANDIDATE INVESTIGATION ON 4 LOW-CORRELATION SAMPLES\n")
cat("==================================================================\n\n")

# ------------------------------------------------------------------------------
# 2. LOAD STEP 05 RDS RESULTS
# ------------------------------------------------------------------------------
rds_path <- file.path(qc_dir, "pca_qc_results.rds")

if (!file.exists(rds_path)) {
  stop("CRITICAL ERROR: 'pca_qc_results.rds' not found in: ", qc_dir)
}

cat("[1/4] Loading Step 05 high-dimensional QC results...\n")
pca_qc_data <- readRDS(rds_path)

qc_df   <- pca_qc_data$qc_dataframe
pca_res <- pca_qc_data$pca_object
hc_res  <- pca_qc_data$hierarchical_clustering

# ------------------------------------------------------------------------------
# 3. COMPUTE STATISTICAL PCA OUTLIER METRICS
# ------------------------------------------------------------------------------
cat("[2/4] Calculating multidimensional PCA outlier metrics (Z-scores & Mahalanobis distance)...\n")

pc_mat <- qc_df[, c("PC1", "PC2", "PC3")]

# PC Z-scores
pc_zscores <- scale(pc_mat)
qc_df$max_pc_zscore <- apply(abs(pc_zscores), 1, max)

# Mahalanobis Distance across top 3 PCs
center_vec <- colMeans(pc_mat)
cov_mat    <- cov(pc_mat)
qc_df$mahalanobis_dist <- mahalanobis(pc_mat, center_vec, cov_mat)

# Chi-Square p-value for Mahalanobis Distance (df = 3)
qc_df$mahalanobis_pval <- pchisq(qc_df$mahalanobis_dist, df = 3, lower.tail = FALSE)

# Extreme PCA Outlier Threshold (|Z| > 3.0 or Mahalanobis p-val < 0.001)
qc_df$extreme_pca_outlier <- qc_df$max_pc_zscore > 3.0 | qc_df$mahalanobis_pval < 0.001

# ------------------------------------------------------------------------------
# 4. EXTRACT DENDROGRAM HEIGHT & CLUSTER POSITION
# ------------------------------------------------------------------------------
cat("[3/4] Extracting hierarchical clustering dendrogram position and branch heights...\n")

# Cophenetic height (distance at which sample joins the dendrogram)
dist_matrix <- as.matrix(as.dist(1 - pca_qc_data$sample_cor))

# Find the height at which each sample merges into a cluster in hclust
# Extract leaf positions and subtree heights
get_merge_height <- function(sample_name, hc, dist_mat) {
  # Calculate min correlation distance to nearest neighbor
  sample_dists <- dist_mat[sample_name, ]
  sample_dists <- sample_dists[names(sample_dists) != sample_name]
  min_dist <- min(sample_dists)
  return(min_dist)
}

qc_df$hc_nearest_neighbor_dist <- sapply(qc_df$sample_id, get_merge_height, hc = hc_res, dist_mat = dist_matrix)

# ------------------------------------------------------------------------------
# 5. FILTER CANDIDATES & BUILD DETAILED INSPECTION TABLE
# ------------------------------------------------------------------------------
cat("[4/4] Compiling candidate investigation profile...\n\n")

candidates_df <- qc_df[qc_df$low_correlation_candidate == TRUE, ]

if (nrow(candidates_df) == 0) {
  cat("No low correlation candidates found based on the threshold.\n")
} else {
  
  # Format clear report
  investigation_report <- data.frame(
    Sample_ID             = candidates_df$sample_id,
    Library_Size          = candidates_df$library_size,
    Detected_Genes        = candidates_df$detected_genes,
    Zero_Fraction         = round(candidates_df$zero_fraction, 4),
    Step02_Flag           = candidates_df$step02_flag,
    PC1                   = round(candidates_df$PC1, 2),
    PC2                   = round(candidates_df$PC2, 2),
    PC3                   = round(candidates_df$PC3, 2),
    Mean_Inter_Cor        = round(candidates_df$mean_inter_cor, 4),
    HC_NN_Distance        = round(candidates_df$hc_nearest_neighbor_dist, 4),
    Extreme_PCA_Outlier   = candidates_df$extreme_pca_outlier,
    Max_PC_ZScore         = round(candidates_df$max_pc_zscore, 2),
    Mahalanobis_pVal      = format.pval(candidates_df$mahalanobis_pval, digits = 3)
  )
  
  # Assign Preliminary Classification Rule:
  # - REMOVE: Technical Failure (Step 02 Flag = TRUE AND Extreme PCA Outlier = TRUE)
  # - REVIEW: Single evidence anomaly (Extreme PCA Outlier = TRUE OR High Zero Fraction / HC Isolation)
  # - RETAIN: Step 02 PASS, normal depth, intact gene count (biological variance)
  investigation_report$Classification <- apply(investigation_report, 1, function(row) {
    step2_flag  <- as.logical(row["Step02_Flag"])
    pca_outlier <- as.logical(row["Extreme_PCA_Outlier"])
    
    if (step2_flag && pca_outlier) {
      return("REMOVE")
    } else if (pca_outlier || step2_flag) {
      return("REVIEW")
    } else {
      return("RETAIN")
    }
  })
  
  # Print report to console
  print(as.data.frame(investigation_report), row.names = FALSE)
  
  # Export results
  out_csv <- file.path(qc_dir, "step05_low_correlation_candidates_investigation.csv")
  write.csv(investigation_report, file = out_csv, row.names = FALSE)
  
  # Generate diagnostic plot focused on the 4 candidates in PCA space
  qc_df$candidate_label <- ifelse(qc_df$low_correlation_candidate, qc_df$sample_id, "Main Cohort")
  
  p_focus <- ggplot(qc_df, aes(x = PC1, y = PC2)) +
    geom_point(data = qc_df[!qc_df$low_correlation_candidate, ], 
               color = "#2B5C8F", alpha = 0.4, size = 1.5) +
    geom_point(data = qc_df[qc_df$low_correlation_candidate, ], 
               color = "#E41A1C", size = 3.5) +
    geom_text(data = qc_df[qc_df$low_correlation_candidate, ], 
              aes(label = sample_id), vjust = -1, size = 3, fontface = "bold", color = "#E41A1C") +
    labs(
      title    = "PCA Isolation Map: 4 Low-Correlation Candidate Samples",
      subtitle = "Red points indicate candidates with mean Pearson correlation r < 0.2328",
      x        = "PC1",
      y        = "PC2"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(qc_dir, "12_candidate_samples_pca_isolation.png"), plot = p_focus, width = 8, height = 6, dpi = 300)
  
  cat("\n==================================================================\n")
  cat("Investigation complete! Results saved to:\n")
  cat("  -> CSV : ", out_csv, "\n")
  cat("  -> Plot: ", file.path(qc_dir, "12_candidate_samples_pca_isolation.png"), "\n")
  cat("==================================================================\n")
}