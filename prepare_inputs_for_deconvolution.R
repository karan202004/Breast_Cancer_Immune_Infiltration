# ==============================================================================
# SCRIPT 05 (TASK 15): PREPARE METHOD-SPECIFIC DECONVOLUTION INPUTS
# Project: TCGA-BRCA Immune Infiltration Analysis
# Location: C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration
# Cohort: N = 1,106 samples | 19,902 genes
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------------------
# 1. SETUP PATHS
# ------------------------------------------------------------------------------
project_dir    <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
normalized_dir <- file.path(project_dir, "results", "normalized_matrix")
deconv_dir     <- file.path(project_dir, "data", "deconvolution_inputs")

dir.create(deconv_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. LOAD VALIDATED STEP 04 MATRICES
# ------------------------------------------------------------------------------
cat("[1/6] Loading validated Step 04 normalized matrices...\n")

linear_cpm_matrix <- readRDS(file.path(normalized_dir, "linear_cpm_1106.rds"))
log2_cpm_matrix   <- readRDS(file.path(normalized_dir, "log2_cpm_1106.rds"))

# Confirm matrix integrity
stopifnot(ncol(linear_cpm_matrix) == 1106)
stopifnot(ncol(log2_cpm_matrix) == 1106)
stopifnot(nrow(linear_cpm_matrix) == 19902)

# ------------------------------------------------------------------------------
# 3. MAP ENSEMBL -> HGNC & APPLY DUPLICATE RULE (slice_max MeanExpr)
# ------------------------------------------------------------------------------
cat("[2/6] Mapping Ensembl IDs to HGNC Symbols & auditing duplicates...\n")

ensembl_ids <- rownames(linear_cpm_matrix)

# Map Ensembl to HGNC Symbol
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
  MeanExpr    = rowMeans(linear_cpm_matrix),
  stringsAsFactors = FALSE
)

# Filtering unmapped genes
mapped_df <- mapping_df %>%
  filter(!is.na(HGNC_Symbol) & HGNC_Symbol != "")

# Apply rule: Keep single Ensembl ID per HGNC symbol with highest mean linear expression
dedup_df <- mapped_df %>%
  group_by(HGNC_Symbol) %>%
  slice_max(order_by = MeanExpr, n = 1, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("   - Total Ensembl genes       : %d\n", nrow(mapping_df)))
cat(sprintf("   - Successfully mapped genes : %d\n", nrow(mapped_df)))
cat(sprintf("   - Unmapped Ensembl genes    : %d\n", nrow(mapping_df) - nrow(mapped_df)))
cat(sprintf("   - Duplicated HGNC symbols   : %d\n", nrow(mapped_df) - nrow(dedup_df)))
cat(sprintf("   - Final unique HGNC genes   : %d\n", nrow(dedup_df)))

# Extract matrices matching deduplicated unique HGNC genes
selected_ensembl <- dedup_df$Ensembl_ID
final_symbols    <- dedup_df$HGNC_Symbol

linear_cpm_hgnc <- linear_cpm_matrix[selected_ensembl, ]
rownames(linear_cpm_hgnc) <- final_symbols

log2_cpm_hgnc <- log2_cpm_matrix[selected_ensembl, ]
rownames(log2_cpm_hgnc) <- final_symbols

# ------------------------------------------------------------------------------
# 4. PREPARE METHOD-SPECIFIC INPUTS
# ------------------------------------------------------------------------------
cat("[3/6] Exporting method-specific input files...\n")

# A. CIBERSORTx (Requires Tab-Delimited TXT with non-log linear values)
cibersort_out <- data.frame(GeneSymbol = rownames(linear_cpm_hgnc), linear_cpm_hgnc, check.names = FALSE)
write.table(
  cibersort_out,
  file      = file.path(deconv_dir, "cibersortx_input_linear_cpm.txt"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

# B. quanTIseq (Non-log linear matrix)
saveRDS(linear_cpm_hgnc, file.path(deconv_dir, "quantiseq_input_linear_cpm.rds"))

# C. MCP-counter (Non-log linear matrix)
saveRDS(linear_cpm_hgnc, file.path(deconv_dir, "mcpcounter_input_linear_cpm.rds"))

# D. xCell (Expects linear or log2 CPM; xCell performs its own SSGSEA ranking)
saveRDS(linear_cpm_hgnc, file.path(deconv_dir, "xcell_input_linear_cpm.rds"))
saveRDS(log2_cpm_hgnc,   file.path(deconv_dir, "xcell_input_log2_cpm.rds"))

# ------------------------------------------------------------------------------
# 5. AUDIT LOG & INPUT QC
# ------------------------------------------------------------------------------
cat("[4/6] Verifying input numerical integrity...\n")

stopifnot(all(!is.na(linear_cpm_hgnc)))
stopifnot(all(!is.infinite(linear_cpm_hgnc)))
stopifnot(all(linear_cpm_hgnc >= 0))

cat("   [PASS] All generated matrices passed numerical integrity checks.\n")

# Save audit table
audit_log <- data.frame(
  Metric = c("Original_Ensembl_Genes", "Mapped_Genes", "Unmapped_Genes", 
             "Duplicated_Symbols_Resolved", "Final_Unique_HGNC_Genes", "Sample_Count"),
  Value  = c(nrow(mapping_df), nrow(mapped_df), nrow(mapping_df) - nrow(mapped_df),
             nrow(mapped_df) - nrow(dedup_df), nrow(dedup_df), ncol(linear_cpm_hgnc))
)
write.csv(audit_log, file.path(deconv_dir, "task15_mapping_audit_summary.csv"), row.names = FALSE)

cat("\n==================================================================\n")
cat("TASK 15 PREPARATION COMPLETE\n")
cat("==================================================================\n")