# ==============================================================================
# Script: 03b_extract_gdc_tpm.R
# Task 15.11: Extract Official GDC TPM from Raw STAR-Counts for EPIC
# Project: TCGA-BRCA Immune Infiltration Analysis
# Cohort: N = 1,106 TCGA-BRCA Primary Tumor Samples
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------------------
# 1. SETUP PATHS & LOAD VALIDATED COHORT METADATA
# ------------------------------------------------------------------------------
project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
norm_dir    <- file.path(project_dir, "results", "normalized_matrix")
deconv_dir  <- file.path(project_dir, "data", "deconvolution_inputs")

dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deconv_dir, recursive = TRUE, showWarnings = FALSE)

cat("[1/5] Loading validated 1,106 sample metadata CSV...\n")
metadata_file <- file.path(project_dir, "metadata", "primary_tumor_1106_metadata.csv")

stopifnot(file.exists(metadata_file))

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Cohort & Metadata Integrity Assertions
stopifnot(nrow(sample_metadata) == 1106)
stopifnot(length(unique(sample_metadata$sample_id)) == 1106)
stopifnot(length(unique(sample_metadata$file_id)) == 1106)

target_paths <- sample_metadata$file_path
if (!all(file.exists(target_paths))) {
  target_paths <- file.path(project_dir, sample_metadata$file_path)
}

missing_files <- target_paths[!file.exists(target_paths)]
if (length(missing_files) > 0) {
  cat("ERROR: The following STAR count files were not found on disk:\n")
  print(head(missing_files, 10))
  stop("Terminating: File path mismatch between metadata CSV and local disk.")
}

cat("   - Successfully verified all 1,106 STAR-Counts files on disk.\n")

# ------------------------------------------------------------------------------
# 2. EXTRACT tpm_unstranded WITH GENE ORDER ASSERTION
# ------------------------------------------------------------------------------
cat("[2/5] Reading tpm_unstranded from 1,106 files and checking gene order...\n")

tpm_list <- vector("list", length = 1106)
names(tpm_list) <- sample_metadata$sample_id

for (i in seq_along(target_paths)) {
  # Column 1 = gene_id, Column 7 = tpm_unstranded
  dt <- fread(target_paths[i], skip = 1, select = c(1, 7))
  colnames(dt) <- c("gene_id", "tpm_unstranded")
  
  # Remove GDC summary metadata rows
  dt <- dt[!grepl("^N_", gene_id)]
  
  # Strip Ensembl version suffix
  dt[, gene_id := sub("\\.\\d+$", "", gene_id)]
  
  # QC Check 1: Verify gene order for every file against File #1
  if (i == 1) {
    gene_ids <- dt$gene_id
    # QC Check 3: Ensure version stripping produced no duplicate Ensembl IDs
    if (anyDuplicated(gene_ids)) {
      stop("Duplicate Ensembl gene IDs detected after version stripping.")
    }
  } else {
    if (!identical(dt$gene_id, gene_ids)) {
      stop("Gene ID or gene order mismatch detected in file: ", target_paths[i])
    }
  }
  
  s_id <- sample_metadata$sample_id[i]
  tpm_list[[s_id]] <- dt$tpm_unstranded
}

# Construct genome-wide TPM matrix
gdc_tpm_raw <- do.call(cbind, tpm_list)
rownames(gdc_tpm_raw) <- gene_ids

stopifnot(ncol(gdc_tpm_raw) == 1106)
stopifnot(all(colnames(gdc_tpm_raw) == sample_metadata$sample_id))

# QC Check 2: Verify finite non-negative values & TPM sum per sample (~1,000,000)
stopifnot(all(is.finite(gdc_tpm_raw)))
stopifnot(all(gdc_tpm_raw >= 0))

tpm_column_sums <- colSums(gdc_tpm_raw)

cat(sprintf("   - Raw Genome-Wide GDC TPM Matrix: %d genes x %d samples\n", 
            nrow(gdc_tpm_raw), ncol(gdc_tpm_raw)))
cat(sprintf("   - TPM column-sum range: %.2f - %.2f\n", 
            min(tpm_column_sums), max(tpm_column_sums)))

# Save provenance matrix
saveRDS(gdc_tpm_raw, file.path(norm_dir, "gdc_tpm_1106.rds"))

# ------------------------------------------------------------------------------
# 3. MAP ENSEMBL TO HGNC SYMBOLS
# ------------------------------------------------------------------------------
cat("[3/5] Mapping Ensembl IDs to HGNC symbols...\n")

ensembl_ids <- rownames(gdc_tpm_raw)

mapped_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = ensembl_ids,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

mapping_df <- data.frame(
  Ensembl_ID  = ensembl_ids,
  HGNC_Symbol = mapped_symbols,
  Mean_TPM    = rowMeans(gdc_tpm_raw),
  stringsAsFactors = FALSE
)

mapped_df <- mapping_df %>% filter(!is.na(HGNC_Symbol) & HGNC_Symbol != "")

# ------------------------------------------------------------------------------
# 4. RESOLVE DUPLICATE SYMBOLS & ENHANCED AUDIT METRICS
# ------------------------------------------------------------------------------
cat("[4/5] Auditing duplicates and applying highest Mean_TPM rule...\n")

duplicate_rows_removed <- nrow(mapped_df) - length(unique(mapped_df$HGNC_Symbol))

duplicate_symbols_summary <- mapped_df %>%
  count(HGNC_Symbol) %>%
  filter(n > 1)

n_duplicate_symbols <- nrow(duplicate_symbols_summary)

dedup_df <- mapped_df %>%
  group_by(HGNC_Symbol) %>%
  slice_max(order_by = Mean_TPM, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("--------------------------------------------------\n")
cat(sprintf("  Genome-wide Ensembl Genes     : %d\n", nrow(mapping_df)))
cat(sprintf("  Mapped Ensembl Rows           : %d\n", nrow(mapped_df)))
cat(sprintf("  Unmapped Ensembl Genes        : %d\n", nrow(mapping_df) - nrow(mapped_df)))
cat(sprintf("  Duplicated HGNC Symbols       : %d\n", n_duplicate_symbols))
cat(sprintf("  Duplicate Ensembl Rows Pruned : %d\n", duplicate_rows_removed))
cat(sprintf("  Final Unique HGNC Matrix      : %d genes x %d samples\n", 
            nrow(dedup_df), ncol(gdc_tpm_raw)))
cat("--------------------------------------------------\n")

selected_ensembl <- dedup_df$Ensembl_ID
final_symbols    <- dedup_df$HGNC_Symbol

epic_tpm_matrix <- gdc_tpm_raw[selected_ensembl, ]
rownames(epic_tpm_matrix) <- final_symbols

# ------------------------------------------------------------------------------
# 5. EXPORT EPIC MATRIX & AUDIT LOG
# ------------------------------------------------------------------------------
cat("[5/5] Performing assertions and exporting epic_input_tpm.rds...\n")

stopifnot(all(!is.na(epic_tpm_matrix)))
stopifnot(all(!is.infinite(epic_tpm_matrix)))
stopifnot(all(epic_tpm_matrix >= 0))
stopifnot(ncol(epic_tpm_matrix) == 1106)

saveRDS(epic_tpm_matrix, file.path(deconv_dir, "epic_input_tpm.rds"))

audit_summary <- data.frame(
  Metric = c("GenomeWide_Ensembl_Genes", "Mapped_Ensembl_Rows", "Unmapped_Ensembl_Genes",
             "Distinct_Duplicated_HGNC_Symbols", "Duplicate_Ensembl_Rows_Pruned", 
             "Final_EPIC_HGNC_Genes", "Sample_Count", 
             "TPM_Column_Sum_Min", "TPM_Column_Sum_Max"),
  Value  = c(nrow(mapping_df), nrow(mapped_df), nrow(mapping_df) - nrow(mapped_df),
             n_duplicate_symbols, duplicate_rows_removed, 
             nrow(epic_tpm_matrix), ncol(epic_tpm_matrix),
             round(min(tpm_column_sums), 2), round(max(tpm_column_sums), 2))
)

write.csv(audit_summary, file.path(deconv_dir, "gdc_tpm_extraction_audit.csv"), row.names = FALSE)

cat("\n======================================================================\n")
cat("  [SUCCESS] GDC TPM EXTRACTION & AUDIT COMPLETE\n")
cat("  Saved EPIC TPM matrix to:\n")
cat(sprintf("  %s/epic_input_tpm.rds\n", deconv_dir))
cat("======================================================================\n")