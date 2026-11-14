# ==============================================================================
# Script: 04_run_batch_deconvolution.R
# Task 16: Local Batch Immune Deconvolution (quanTIseq, EPIC, MCP-counter, xCell)
# Project: TCGA-BRCA Immune Infiltration Analysis
# Cohort: N = 1,106 TCGA-BRCA Primary Tumor Samples
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(immunedeconv)
})

# ------------------------------------------------------------------------------
# 1. SETUP PATHS
# ------------------------------------------------------------------------------
project_dir <- file.path("C:/Users/karan/OneDrive/Documents", "Breast_Cancer_Immune_Infiltration")
input_dir   <- file.path(project_dir, "data", "deconvolution_inputs")
meta_file   <- file.path(project_dir, "metadata", "primary_tumor_1106_metadata.csv")
output_dir  <- file.path(project_dir, "results", "deconvolution_results")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. VALIDATE SAMPLE METADATA & INPUT FILE
# ------------------------------------------------------------------------------
cat("[1/4] Verifying metadata and input TPM matrix...\n")

stopifnot(file.exists(meta_file))
metadata <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(metadata) == 1106)
target_samples <- metadata$sample_id

tpm_file <- file.path(input_dir, "epic_input_tpm.rds")
stopifnot(file.exists(tpm_file))

# Load the validated non-log HGNC TPM matrix
tpm_mat <- readRDS(tpm_file)

# Structural & numerical assertions
stopifnot(
  is.matrix(tpm_mat) || is.data.frame(tpm_mat),
  ncol(tpm_mat) == 1106,
  !anyNA(tpm_mat),
  all(is.finite(as.matrix(tpm_mat))),
  all(as.matrix(tpm_mat) >= 0),
  setequal(colnames(tpm_mat), target_samples)
)

# Enforce strict column order matching sample_metadata
tpm_mat <- as.matrix(tpm_mat)
tpm_mat <- tpm_mat[, target_samples, drop = FALSE]
stopifnot(identical(colnames(tpm_mat), target_samples))

cat(sprintf("   - Input Matrix Verified: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

# ------------------------------------------------------------------------------
# 3. UNIFIED INPUT MAPPING ACROSS ALGORITHMS
# ------------------------------------------------------------------------------
# immunedeconv handles internal scale requirements when provided non-log TPM
methods_to_run <- c("quantiseq", "epic", "mcp_counter", "xcell")

# ------------------------------------------------------------------------------
# 4. EXECUTION ENGINE
# ------------------------------------------------------------------------------
run_method <- function(method_name) {
  cat("\n==================================================")
  cat(sprintf("\nRunning Deconvolution: %s\n", method_name))
  cat("==================================================\n")
  
  tryCatch({
    # Execute method via immunedeconv::deconvolute
    result <- switch(
      method_name,
      
      quantiseq = immunedeconv::deconvolute(
        gene_expression = tpm_mat,
        method          = "quantiseq",
        tumor           = TRUE,
        arrays          = FALSE,
        scale_mrna      = TRUE
      ),
      
      epic = immunedeconv::deconvolute(
        gene_expression = tpm_mat,
        method          = "epic",
        tumor           = TRUE,
        scale_mrna      = TRUE
      ),
      
      mcp_counter = immunedeconv::deconvolute(
        gene_expression = tpm_mat,
        method          = "mcp_counter"
      ),
      
      xcell = immunedeconv::deconvolute(
        gene_expression = tpm_mat,
        method          = "xcell"
      )
    )
    
    # Save RDS matrix object
    rds_out <- file.path(output_dir, paste0(method_name, "_results_1106.rds"))
    saveRDS(result, rds_out)
    
    # Save CSV spreadsheet (cell_type in col 1, 1106 sample columns)
    csv_out <- file.path(output_dir, paste0(method_name, "_results_1106.csv"))
    write.csv(as.data.frame(result), csv_out, row.names = FALSE)
    
    n_cells   <- nrow(result)
    n_samples <- ncol(result) - 1 # subtracting 'cell_type' column
    
    cat(sprintf("[SUCCESS] %s execution complete.\n", method_name))
    cat(sprintf("Output shape: %d cell types x %d samples (+1 label column)\n", n_cells, n_samples))
    
    list(
      method  = method_name,
      status  = "SUCCESS",
      message = "",
      result  = result
    )
    
  }, error = function(e) {
    cat(sprintf("[FAILED] %s execution encountered an error.\n", method_name))
    cat("Error Message:", conditionMessage(e), "\n")
    
    list(
      method  = method_name,
      status  = "FAILED",
      message = conditionMessage(e),
      result  = NULL
    )
  })
}

# ------------------------------------------------------------------------------
# 5. EXECUTE BATCH PIPELINE & SAVE AUDIT LOG
# ------------------------------------------------------------------------------
cat("\n[2/4] Starting batch execution across 4 methods...\n")
start_time <- Sys.time()

batch_results <- setNames(
  lapply(methods_to_run, run_method),
  methods_to_run
)

end_time <- Sys.time()
cat(sprintf("\nTotal Execution Time: %.2f minutes\n", as.numeric(difftime(end_time, start_time, units = "mins"))))

# Export Summary Audit
cat("\n[3/4] Exporting execution summary log...\n")
run_summary <- data.frame(
  method  = methods_to_run,
  status  = vapply(batch_results, function(x) x$status, character(1)),
  message = vapply(batch_results, function(x) x$message, character(1)),
  stringsAsFactors = FALSE
)

summary_file <- file.path(output_dir, "deconvolution_run_summary.csv")
write.csv(run_summary, summary_file, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. PRINT EXECUTION REPORT
# ------------------------------------------------------------------------------
cat("\n[4/4] Execution Complete. Summary:\n")
cat("======================================================================\n")
print(run_summary)
cat("======================================================================\n")
cat(sprintf("All results saved to: %s\n", output_dir))