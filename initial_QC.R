#Initial QC
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

project_dir   <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
expr_dir      <- file.path(project_dir, "results", "expression_matrix")
metadata_path <- file.path(project_dir, "metadata", "primary_tumor_1106_metadata.csv")
qc_dir        <- file.path(project_dir, "results", "QC")

#Create QC directory
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

cat("Sample level intial qc\n")

cat("Loading expression matrix, annotation and metadata...\n")

raw_counts <- readRDS(
  file.path(expr_dir, "raw_counts_1106.rds")
)

gene_annotation <- readRDS(
  file.path(expr_dir, "gene_annotation_1106.rds")
)

metadata <- read.csv(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

n_samples <- ncol(raw_counts)
n_genes   <- nrow(raw_counts)

cat("  Expression matrix:",
    format(n_genes, big.mark = ","),
    "genes x",
    format(n_samples, big.mark = ","),
    "samples\n")

#STRUCTURAL & INTEGRITY VALIDATION
cat("\n Validating matrix and metadata integrity...\n")

#Expected dimensions
stopifnot(n_genes == 60616)
stopifnot(n_samples == 1106)

#Metadata dimensions
stopifnot(nrow(metadata) == 1106)

#sample IDs must be unique
stopifnot(!anyDuplicated(colnames(raw_counts)))
stopifnot(!anyDuplicated(metadata$sample_id))

stopifnot(
  identical(
    colnames(raw_counts),
    metadata$sample_id
  )
)

# Gene IDs in expression matrix must be unique
stopifnot(!anyDuplicated(rownames(raw_counts)))

#Check gene annotation alignment
if (!identical(rownames(raw_counts), rownames(gene_annotation))) {
  cat("  [!] Annotation row names do not match raw count row names identically.\n")
  cat("      Attempting alignment verification...\n")
  
  #Check if all raw_count genes exist in gene_annotation
  if (all(rownames(raw_counts) %in% rownames(gene_annotation))) {
    cat("      -> All count matrix genes present in annotation. Reordering annotation...\n")
    gene_annotation <- gene_annotation[rownames(raw_counts), , drop = FALSE]
  } else if ("gene_id" %in% colnames(gene_annotation) && all(rownames(raw_counts) %in% gene_annotation$gene_id)) {
    cat("      -> Matching via 'gene_id' column in annotation...\n")
    rownames(gene_annotation) <- gene_annotation$gene_id
    gene_annotation <- gene_annotation[rownames(raw_counts), , drop = FALSE]
  } else {
    stop("CRITICAL ERROR: Gene IDs in raw_counts do not set-match gene_annotation!")
  }
}

#Final assertion after alignment
stopifnot(identical(rownames(raw_counts), rownames(gene_annotation)))

#No missing values
stopifnot(!anyNA(raw_counts))

#No negative counts
stopifnot(!any(raw_counts < 0))

#Counts should be numeric
stopifnot(is.numeric(raw_counts))

cat("sample IDs: PASS\n")
cat("gene IDs: PASS\n")
cat("annotation alignment: PASS (Aligned)\n")
cat("no missing values: PASS\n")
cat("no negative counts: PASS\n")
cat("matrix dimensions: PASS\n")

#sample level qc
cat("\n Computing sample-level sequencing metrics...\n")

#Total counts
library_sizes <- colSums(raw_counts)

#Number of detected genes (count > 0)
detected_genes <- colSums(raw_counts > 0)

#Number of zero-count genes
zero_genes <- colSums(raw_counts == 0)

#Fraction of genes with zero counts
zero_fractions <- zero_genes / n_genes

#Fraction of detected genes
detected_fraction <- detected_genes / n_genes

#Log10 library size (for robust outlier assessment)
log10_library_size <- log10(library_sizes)

# 5. Create qc table

qc_metrics <- data.frame(
  sample_id             = colnames(raw_counts),
  file_id               = metadata$file_id,
  case_id               = metadata$case_id,
  library_size          = as.numeric(library_sizes),
  detected_genes        = as.numeric(detected_genes),
  zero_genes            = as.numeric(zero_genes),
  zero_fraction         = as.numeric(zero_fractions),
  detected_fraction     = as.numeric(detected_fraction),
  log10_library_size    = as.numeric(log10_library_size),
  stringsAsFactors      = FALSE
)


cat("\n initial QC \n")

cat("\nNumber of samples evaluated:", format(n_samples, big.mark = ","), "\n")
cat("Number of genes evaluated:  ", format(n_genes, big.mark = ","), "\n")

cat("\nLibrary Size / Total Counts:\n")
cat("  Minimum :", format(min(library_sizes), big.mark = ","), "\n")
cat("  Median  :", format(median(library_sizes), big.mark = ","), "\n")
cat("  Mean    :", format(round(mean(library_sizes)), big.mark = ","), "\n")
cat("  Maximum :", format(max(library_sizes), big.mark = ","), "\n")

cat("\nDetected Genes:\n")
cat("  Minimum :", format(min(detected_genes), big.mark = ","), "\n")
cat("  Median  :", format(median(detected_genes), big.mark = ","), "\n")
cat("  Mean    :", format(round(mean(detected_genes)), big.mark = ","), "\n")
cat("  Maximum :", format(max(detected_genes), big.mark = ","), "\n")

cat("\nZero-Count Fraction:\n")
cat("  Minimum :", round(min(zero_fractions), 4), "\n")
cat("  Median  :", round(median(zero_fractions), 4), "\n")
cat("  Mean    :", round(mean(zero_fractions), 4), "\n")
cat("  Maximum :", round(max(zero_fractions), 4), "\n")

cat("\n Evaluating potential technical outliers...\n")

#Flag 1: Hard cutoff below 10M counts (Flag only)
low_depth_flag <- library_sizes < 10000000
cat("  Samples with library size < 10,000,000:", sum(low_depth_flag), "\n")

#Flag 2: Robust MAD cutoff on log10 scale
lib_log_median <- median(log10_library_size)
lib_log_mad    <- mad(log10_library_size, center = lib_log_median)
mad_threshold  <- lib_log_median - (3 * lib_log_mad)

mad_low_depth_flag <- log10_library_size < mad_threshold
cat("  Samples flagged by 3-MAD log10 library-size criterion:", sum(mad_low_depth_flag), "\n")

#Assign flags to dataframe
qc_metrics$low_depth_flag     <- low_depth_flag
qc_metrics$mad_low_depth_flag <- mad_low_depth_flag
qc_metrics$initial_technical_flag <- qc_metrics$low_depth_flag | qc_metrics$mad_low_depth_flag

cat("  Samples with at least one initial technical flag:", sum(qc_metrics$initial_technical_flag), "\n")

cat("\n Calculating raw count distribution \n")

total_count_values <- as.numeric(raw_counts)

raw_count_summary <- data.frame(
  statistic = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"),
  value     = c(
    min(total_count_values),
    quantile(total_count_values, 0.25),
    median(total_count_values),
    mean(total_count_values),
    quantile(total_count_values, 0.75),
    max(total_count_values)
  )
)

rm(total_count_values)
cat("  Raw count distribution summary calculated.\n")

#Generate QC plot
cat("\n Generating diagnostic QC plots...\n")

# Plot 1: Library size distribution
p1 <- ggplot(qc_metrics, aes(x = library_size / 1e6)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white") +
  geom_vline(xintercept = 10, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "TCGA-BRCA Library Size Distribution",
    subtitle = "1,106 Primary Tumor RNA-seq samples",
    x = "Library Size (Million Counts)",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(qc_dir, "01_library_size_distribution.png"), plot = p1, width = 8, height = 5, dpi = 300)

# Plot 2: Library size vs detected genes
p2 <- ggplot(qc_metrics, aes(x = library_size / 1e6, y = detected_genes)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = median(detected_genes), linetype = "dotted", linewidth = 0.8) +
  labs(
    title = "Library Size vs. Detected Genes",
    subtitle = "Horizontal line represents median detected-gene count",
    x = "Library Size (Million Counts)",
    y = "Detected Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(qc_dir, "02_library_size_vs_detected_genes.png"), plot = p2, width = 8, height = 5, dpi = 300)

# Plot 3: Zero-count fraction distribution
p3 <- ggplot(qc_metrics, aes(x = zero_fraction)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white") +
  labs(
    title = "Zero-Count Fraction Across Samples",
    subtitle = "Fraction of 60,616 raw-count features with zero counts",
    x = "Zero-Count Fraction",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(qc_dir, "03_zero_fraction_distribution.png"), plot = p3, width = 8, height = 5, dpi = 300)

# Plot 4: Detected genes distribution
p4 <- ggplot(qc_metrics, aes(x = detected_genes)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white") +
  geom_vline(xintercept = median(detected_genes), linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Detected Gene Distribution Across Samples",
    subtitle = "Genes with raw count > 0",
    x = "Detected Genes",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(qc_dir, "04_detected_gene_distribution.png"), plot = p4, width = 8, height = 5, dpi = 300)

# Plot 5: Log10 library size distribution
p5 <- ggplot(qc_metrics, aes(x = log10_library_size)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white") +
  geom_vline(xintercept = lib_log_median, linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Log10 Library Size Distribution",
    subtitle = "Used for robust MAD-based technical outlier assessment",
    x = "log10(Total Counts)",
    y = "Number of Samples"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(qc_dir, "05_log10_library_size_distribution.png"), plot = p5, width = 8, height = 5, dpi = 300)

# Save QC outputs
cat("\nExporting QC metrics and summary files...\n")

write.csv(qc_metrics, file = file.path(qc_dir, "sample_qc_metrics_1106.csv"), row.names = FALSE)
saveRDS(qc_metrics, file = file.path(qc_dir, "sample_qc_metrics_1106.rds"))
write.csv(raw_count_summary, file = file.path(qc_dir, "raw_count_distribution_summary.csv"), row.names = FALSE)


summary_file <- file.path(qc_dir, "processing_summary_02.txt")
sink(summary_file)


cat("Samples evaluated: ", format(n_samples, big.mark = ","), "\n", sep = "")
cat("Genes evaluated:   ", format(n_genes, big.mark = ","), "\n", sep = "")

cat("\n Structural validation\n")
cat("Matrix dimensions: PASS\n")
cat("Sample ID alignment: PASS\n")
cat("Gene ID uniqueness: PASS\n")
cat("Gene annotation alignment: PASS\n")
cat("No NA values: PASS\n")
cat("No negative counts: PASS\n")

cat("\nLibrary size\n")
cat("------------\n")
cat("Minimum: ", format(min(library_sizes), big.mark = ","), "\n", sep = "")
cat("Median:  ", format(median(library_sizes), big.mark = ","), "\n", sep = "")
cat("Mean:    ", format(round(mean(library_sizes)), big.mark = ","), "\n", sep = "")
cat("Maximum: ", format(max(library_sizes), big.mark = ","), "\n", sep = "")

cat("\nDetected genes \n")
cat("--------------\n")
cat("Minimum: ", format(min(detected_genes), big.mark = ","), "\n", sep = "")
cat("Median:  ", format(median(detected_genes), big.mark = ","), "\n", sep = "")
cat("Mean:    ", format(round(mean(detected_genes)), big.mark = ","), "\n", sep = "")
cat("Maximum: ", format(max(detected_genes), big.mark = ","), "\n", sep = "")

cat("\nzero count fraction\n")
cat("-------------------\n")
cat("Minimum: ", round(min(zero_fractions), 4), "\n", sep = "")
cat("Median:  ", round(median(zero_fractions), 4), "\n", sep = "")
cat("Mean:    ", round(mean(zero_fractions), 4), "\n", sep = "")
cat("Maximum: ", round(max(zero_fractions), 4), "\n", sep = "")

cat("\nTechnical flags\n")
cat("---------------\n")
cat("Library size < 10,000,000: ", sum(low_depth_flag), "\n", sep = "")
cat("3-MAD low-depth flag:      ", sum(mad_low_depth_flag), "\n", sep = "")
cat("Any initial technical flag: ", sum(qc_metrics$initial_technical_flag), "\n", sep = "")

cat("\nsample removal\n")
cat("--------------\n")
cat("Samples removed in Step 02: 0\n")
cat("All samples retained for downstream normalized-data QC.\n")

cat("\noutput files\n")
cat("sample_qc_metrics_1106.csv\n")
cat("sample_qc_metrics_1106.rds\n")
cat("raw_count_distribution_summary.csv\n")
cat("01_library_size_distribution.png\n")
cat("02_library_size_vs_detected_genes.png\n")
cat("03_zero_fraction_distribution.png\n")
cat("04_detected_gene_distribution.png\n")
cat("05_log10_library_size_distribution.png\n")
cat("processing_summary_02.txt\n")

sink()

cat("QC metrics, diagnostic plots, and processing log saved to:\n")
cat(qc_dir, "\n")
