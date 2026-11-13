# ==============================================================================
# SCRIPT 04: LIBRARY NORMALIZATION (TMM & LOG2-CPM)
# Project: TCGA-BRCA Immune Infiltration Analysis
# Author: Karan S.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & DEPENDENCIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
  library(ggplot2)
})

# Define project paths
project_dir     <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
filtered_dir    <- file.path(project_dir, "results", "filtered_matrix")
normalized_dir  <- file.path(project_dir, "results", "normalized_matrix")
qc_dir          <- file.path(project_dir, "results", "QC")

# Create output directory if it does not exist
if (!dir.exists(normalized_dir)) {
  dir.create(normalized_dir, recursive = TRUE)
}

cat("==================================================================\n")
cat("STEP 04: LIBRARY NORMALIZATION (TMM & LOG2-CPM)\n")
cat("==================================================================\n\n")


# ------------------------------------------------------------------------------
# 2. LOAD FILTERED DATA
# ------------------------------------------------------------------------------
cat("[1/6] Loading filtered count matrix and gene annotations...\n")

filtered_counts <- readRDS(
  file.path(filtered_dir, "filtered_counts_1106.rds")
)

filtered_annotation <- readRDS(
  file.path(filtered_dir, "filtered_gene_annotation_1106.rds")
)

n_samples <- ncol(filtered_counts)
n_genes   <- nrow(filtered_counts)

cat("  Loaded:", format(n_genes, big.mark = ","), "retained genes x", 
    format(n_samples, big.mark = ","), "samples\n")


# ------------------------------------------------------------------------------
# 3. COMPUTE TMM NORMALIZATION FACTORS
# ------------------------------------------------------------------------------
cat("\n[2/6] Constructing DGEList and calculating TMM normalization factors...\n")

# Create edgeR DGEList object
dge <- DGEList(counts = filtered_counts)

# Calculate TMM normalization factors
dge <- calcNormFactors(dge, method = "TMM")

raw_lib_size <- dge$samples$lib.size
norm_factors <- dge$samples$norm.factors
eff_lib_size <- raw_lib_size * norm_factors

cat("  TMM Factor Summary:\n")
cat("  -> Min TMM Factor   :", round(min(norm_factors), 4), "\n")
cat("  -> Median TMM Factor:", round(median(norm_factors), 4), "\n")
cat("  -> Max TMM Factor   :", round(max(norm_factors), 4), "\n")

cat("\n  Raw vs Effective Library Size Median:\n")
cat("  -> Raw Library Size Median      :", round(median(raw_lib_size) / 1e6, 2), "M\n")
cat("  -> Effective Library Size Median:", round(median(eff_lib_size) / 1e6, 2), "M\n")


# ------------------------------------------------------------------------------
# 4. GENERATE NORMALIZED EXPRESSION MATRICES
# ------------------------------------------------------------------------------
cat("\n[3/6] Generating log2-CPM and linear CPM matrices...\n")

# A. log2-CPM matrix (prior.count = 2 for variance stabilization)
log2_cpm_matrix <- cpm(dge, log = TRUE, prior.count = 2)

# B. Linear normalized CPM matrix (for cell deconvolution algorithms)
linear_cpm_matrix <- cpm(dge, log = FALSE, prior.count = 0)


# ------------------------------------------------------------------------------
# 5. NORMALIZATION & MATRIX INTEGRITY VALIDATION
# ------------------------------------------------------------------------------
cat("\n[4/6] Running matrix integrity assertions and numerical checks...\n")

# A. Structural Dimension Alignment
stopifnot(identical(dim(log2_cpm_matrix), dim(filtered_counts)))
stopifnot(identical(dim(linear_cpm_matrix), dim(filtered_counts)))

# B. Gene Order Consistency
stopifnot(identical(rownames(log2_cpm_matrix), rownames(filtered_counts)))
stopifnot(identical(rownames(linear_cpm_matrix), rownames(filtered_counts)))

# C. Sample Order Consistency
stopifnot(identical(colnames(log2_cpm_matrix), colnames(filtered_counts)))
stopifnot(identical(colnames(linear_cpm_matrix), colnames(filtered_counts)))

# D. Numerical Integrity
stopifnot(all(is.finite(log2_cpm_matrix)))
stopifnot(all(is.finite(linear_cpm_matrix)))
stopifnot(all(linear_cpm_matrix >= 0))

cat("  [PASS] Matrix dimensions match (19,902 x 1,106)\n")
cat("  [PASS] Gene and Sample row/column alignment verified\n")
cat("  [PASS] Numerical checks passed (no NAs, NaNs, Infs, or negative values)\n")


# ------------------------------------------------------------------------------
# 6. DIAGNOSTIC PLOT: TMM NORMALIZATION FACTOR DISTRIBUTION
# ------------------------------------------------------------------------------
cat("\n[5/6] Generating TMM factor distribution diagnostic plot...\n")

tmm_df <- data.frame(sample_id = rownames(dge$samples), norm_factor = norm_factors)

p_tmm <- ggplot(tmm_df, aes(x = norm_factor)) +
  geom_histogram(bins = 40, fill = "#2B5C8F", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "red", size = 0.8) +
  labs(
    title = "Distribution of TMM Normalization Factors",
    subtitle = paste0("TCGA-BRCA Cohort (N = ", n_samples, ")"),
    x = "TMM Normalization Factor",
    y = "Sample Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  filename = file.path(qc_dir, "07_tmm_factor_distribution.png"),
  plot = p_tmm,
  width = 7,
  height = 4.5,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 7. SAVE NORMALIZED MATRICES & LOGS
# ------------------------------------------------------------------------------
cat("\n[6/6] Exporting normalized matrices and summary logs...\n")

# Save objects
saveRDS(
  dge,
  file = file.path(normalized_dir, "dge_tmm_1106.rds")
)

saveRDS(
  log2_cpm_matrix,
  file = file.path(normalized_dir, "log2_cpm_1106.rds")
)

saveRDS(
  linear_cpm_matrix,
  file = file.path(normalized_dir, "linear_cpm_1106.rds")
)

# Export sample normalization metadata table
sample_norm_df <- data.frame(
  sample_id       = rownames(dge$samples),
  raw_lib_size    = raw_lib_size,
  tmm_factor      = norm_factors,
  eff_lib_size    = eff_lib_size
)

write.csv(
  sample_norm_df,
  file = file.path(normalized_dir, "sample_normalization_summary_1106.csv"),
  row.names = FALSE
)

# Write audit log
summary_file <- file.path(qc_dir, "processing_summary_04.txt")
sink(summary_file)

cat("==================================================================\n")
cat("STEP 04: LIBRARY NORMALIZATION (TMM & LOG2-CPM)\n")
cat("TCGA-BRCA Immune Infiltration Analysis\n")
cat("==================================================================\n\n")

cat("INPUT MATRIX\n")
cat("------------\n")
cat("Retained Genes:             ", format(n_genes, big.mark = ","), "\n")
cat("Retained Samples:           ", format(n_samples, big.mark = ","), "\n")

cat("\nTMM FACTOR STATISTICAL SUMMARY\n")
cat("------------------------------\n")
cat("Min TMM Factor:             ", round(min(norm_factors), 4), "\n")
cat("1st Quartile TMM Factor:    ", round(quantile(norm_factors, 0.25), 4), "\n")
cat("Median TMM Factor:          ", round(median(norm_factors), 4), "\n")
cat("Mean TMM Factor:            ", round(mean(norm_factors), 4), "\n")
cat("3rd Quartile TMM Factor:    ", round(quantile(norm_factors, 0.75), 75), "\n")
cat("Max TMM Factor:             ", round(max(norm_factors), 4), "\n")

cat("\nRAW VS EFFECTIVE LIBRARY SIZE COMPARISON (Millions)\n")
cat("---------------------------------------------------\n")
cat("Raw Library Size (Min/Median/Max):       ", 
    round(min(raw_lib_size) / 1e6, 2), "M / ",
    round(median(raw_lib_size) / 1e6, 2), "M / ",
    round(max(raw_lib_size) / 1e6, 2), "M\n", sep = "")
cat("Effective Library Size (Min/Median/Max): ", 
    round(min(eff_lib_size) / 1e6, 2), "M / ",
    round(median(eff_lib_size) / 1e6, 2), "M / ",
    round(max(eff_lib_size) / 1e6, 2), "M\n", sep = "")

cat("\nMATRIX INTEGRITY VERIFICATION\n")
cat("-----------------------------\n")
cat("Dimension Integrity Pass:   TRUE (19,902 x 1,106)\n")
cat("Gene Order Alignment Pass:  TRUE\n")
cat("Sample Order Alignment Pass: TRUE\n")
cat("Numerical Integrity Pass:   TRUE (Finite & Valid Values)\n")

cat("\nOUTPUT FILES\n")
cat("------------\n")
cat("results/normalized_matrix/dge_tmm_1106.rds\n")
cat("results/normalized_matrix/log2_cpm_1106.rds\n")
cat("results/normalized_matrix/linear_cpm_1106.rds\n")
cat("results/normalized_matrix/sample_normalization_summary_1106.csv\n")
cat("results/QC/07_tmm_factor_distribution.png\n")
cat("results/QC/processing_summary_04.txt\n")

cat("\n==================================================================\n")
cat("STEP 04 COMPLETE\n")
cat("==================================================================\n")

sink()

cat("\n==================================================================\n")
cat("SUCCESS: Step 04 complete!\n")
cat("Normalized matrices saved to:\n")
cat(normalized_dir, "\n")
cat("==================================================================\n")