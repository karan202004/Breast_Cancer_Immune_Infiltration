# CIBERSORT (LM22 deconvolution)

project_dir <- file.path(
  "C:/Users/karan/OneDrive/Documents",
  "Breast_Cancer_Immune_Infiltration"
)

ciber_dir <- file.path(project_dir, "data", "CIBERSORT")
output_dir <- file.path(project_dir, "results", "deconvolution_results")

script_file <- file.path(ciber_dir, "CIBERSORT.R")
sig_file <- file.path(ciber_dir, "LM22.txt")
mixture_file <- file.path(ciber_dir, "BRCA_TPM_1106.txt")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(script_file),
  file.exists(sig_file),
  file.exists(mixture_file)
)

# Load the original CIBERSORT function
source(script_file)

stopifnot(exists("CIBERSORT", mode = "function"))

#Run with 100 permutations; disable quantile normalization for RNA-seq TPM
setwd(output_dir)

cibersort_results <- CIBERSORT(
  sig_matrix = sig_file,
  mixture_file = mixture_file,
  perm = 100,
  QN = FALSE
)

#Save R object
saveRDS(
  cibersort_results,
  file = file.path(output_dir, "cibersort_LM22_results_1106.rds")
)

#Export tabular results
write.csv(
  data.frame(
    sample_id = rownames(cibersort_results),
    cibersort_results,
    check.names = FALSE
  ),
  file = file.path(output_dir, "cibersort_LM22_results_1106.csv"),
  row.names = FALSE
)

cat("CIBERSORT finished.\n")
cat("Result dimensions:", dim(cibersort_results), "\n")
cat("Results saved in:", output_dir, "\n")
