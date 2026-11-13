#Normalization

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
  library(ggplot2)
})

project_dir     <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
filtered_dir    <- file.path(project_dir, "results", "filtered_matrix")
normalized_dir  <- file.path(project_dir, "results", "normalized_matrix")
qc_dir          <- file.path(project_dir, "results", "QC")

if (!dir.exists(normalized_dir)) {
  dir.create(normalized_dir, recursive = TRUE)
}

cat("Loading filtered count matrix and gene annotations\n")

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

#Compute tmm normalization factors
cat("\n Constructing DGEList and calculating TMM normalization factors\n")

#Create edgeR DGEList object
dge <- DGEList(counts = filtered_counts)

#Calculate TMM normalization factors
dge <- calcNormFactors(dge, method = "TMM")

raw_lib_size <- dge$samples$lib.size
norm_factors <- dge$samples$norm.factors
eff_lib_size <- raw_lib_size * norm_factors

cat("  TMM factor summary:\n")
cat("  -> Min TMM factor   :", round(min(norm_factors), 4), "\n")
cat("  -> Median TMM factor:", round(median(norm_factors), 4), "\n")
cat("  -> Max TMM factor   :", round(max(norm_factors), 4), "\n")

cat("\n  Raw vs Effective Library Size Median:\n")
cat("  -> Raw Library Size Median      :", round(median(raw_lib_size) / 1e6, 2), "M\n")
cat("  -> Effective Library Size Median:", round(median(eff_lib_size) / 1e6, 2), "M\n")

#Generate normalized expression matrices
cat("\n Generating log2-CPM and linear CPM matrices...\n")

#log2-CPM matrix (prior.count = 2 for variance stabilization)
log2_cpm_matrix <- cpm(dge, log = TRUE, prior.count = 2)

#Linear normalized CPM matrix (for cell deconvolution algorithms)
linear_cpm_matrix <- cpm(dge, log = FALSE, prior.count = 0)

#Normalization & matrix integrity validation
cat("\n Running matrix integrity assertions and numerical checks...\n")

#Structural Dimension Alignment
stopifnot(identical(dim(log2_cpm_matrix), dim(filtered_counts)))
stopifnot(identical(dim(linear_cpm_matrix), dim(filtered_counts)))

#Gene Order Consistency
stopifnot(identical(rownames(log2_cpm_matrix), rownames(filtered_counts)))
stopifnot(identical(rownames(linear_cpm_matrix), rownames(filtered_counts)))

#Sample Order Consistency
stopifnot(identical(colnames(log2_cpm_matrix), colnames(filtered_counts)))
stopifnot(identical(colnames(linear_cpm_matrix), colnames(filtered_counts)))

#Numerical Integrity
stopifnot(all(is.finite(log2_cpm_matrix)))
stopifnot(all(is.finite(linear_cpm_matrix)))
stopifnot(all(linear_cpm_matrix >= 0))

cat("  [PASS] Matrix dimensions match (19,902 x 1,106)\n")
cat("  [PASS] Gene and Sample row/column alignment verified\n")
cat("  [PASS] Numerical checks passed (no NAs, NaNs, Infs, or negative values)\n")


#Diagnostic plot: tmm normalization factor distribution
cat("\n Generating TMM factor distribution diagnostic plot...\n")

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

#Save normalized matrices & logs
cat("\n Exporting normalized matrices and summary logs...\n")

#Save objects
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

#Export sample normalization metadata table
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

summary_file <- file.path(qc_dir, "processing_summary_04.txt")
sink(summary_file)

cat("Retained Genes:             ", format(n_genes, big.mark = ","), "\n")
cat("Retained Samples:           ", format(n_samples, big.mark = ","), "\n")

cat("\nTMM factor statistical Summary\n")
cat("------------------------------\n")
cat("Min TMM Factor:             ", round(min(norm_factors), 4), "\n")
cat("1st Quartile TMM Factor:    ", round(quantile(norm_factors, 0.25), 4), "\n")
cat("Median TMM Factor:          ", round(median(norm_factors), 4), "\n")
cat("Mean TMM Factor:            ", round(mean(norm_factors), 4), "\n")
cat("3rd Quartile TMM Factor:    ", round(quantile(norm_factors, 0.75), 75), "\n")
cat("Max TMM Factor:             ", round(max(norm_factors), 4), "\n")

cat("\nraw vs effective library size comparison (millions)\n")
cat("---------------------------------------------------\n")
cat("Raw Library Size (Min/Median/Max):       ", 
    round(min(raw_lib_size) / 1e6, 2), "M / ",
    round(median(raw_lib_size) / 1e6, 2), "M / ",
    round(max(raw_lib_size) / 1e6, 2), "M\n", sep = "")
cat("Effective Library Size (Min/Median/Max): ", 
    round(min(eff_lib_size) / 1e6, 2), "M / ",
    round(median(eff_lib_size) / 1e6, 2), "M / ",
    round(max(eff_lib_size) / 1e6, 2), "M\n", sep = "")

cat("\nmatrix integrity verification\n")
cat("-----------------------------\n")
cat("Dimension Integrity Pass:   TRUE (19,902 x 1,106)\n")
cat("Gene Order Alignment Pass:  TRUE\n")
cat("Sample Order Alignment Pass: TRUE\n")
cat("Numerical Integrity Pass:   TRUE (Finite & Valid Values)\n")

cat("\nOutput files\n")
cat("------------\n")
cat("results/normalized_matrix/dge_tmm_1106.rds\n")
cat("results/normalized_matrix/log2_cpm_1106.rds\n")
cat("results/normalized_matrix/linear_cpm_1106.rds\n")
cat("results/normalized_matrix/sample_normalization_summary_1106.csv\n")
cat("results/QC/07_tmm_factor_distribution.png\n")
cat("results/QC/processing_summary_04.txt\n")

sink()

cat("Normalized matrices saved to:\n")
cat(normalized_dir, "\n")
