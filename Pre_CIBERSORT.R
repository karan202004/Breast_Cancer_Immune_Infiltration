
# ==============================================================================
# CIBERSORT - Step 3: Prepare and validate input files
# ==============================================================================

project_dir <- file.path(
  "C:/Users/karan/OneDrive/Documents",
  "Breast_Cancer_Immune_Infiltration"
)

ciber_dir <- file.path(project_dir, "data", "CIBERSORT")
input_file <- file.path(
  project_dir, "data", "deconvolution_inputs",
  "epic_input_tpm.rds"
)
meta_file <- file.path(
  project_dir, "metadata",
  "primary_tumor_1106_metadata.csv"
)

dir.create(ciber_dir, recursive = TRUE, showWarnings = FALSE)

sig_file <- file.path(ciber_dir, "LM22.txt")
script_file <- file.path(ciber_dir, "CIBERSORT.R")
mixture_file <- file.path(ciber_dir, "BRCA_TPM_1106.txt")

stopifnot(
  file.exists(sig_file),
  file.exists(script_file),
  file.exists(input_file),
  file.exists(meta_file)
)

# Load expression matrix and metadata
tpm_mat <- as.matrix(readRDS(input_file))
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)

stopifnot(
  nrow(metadata) == 1106,
  ncol(tpm_mat) == 1106,
  !is.null(rownames(tpm_mat)),
  !is.null(colnames(tpm_mat)),
  !anyDuplicated(rownames(tpm_mat)),
  !anyDuplicated(colnames(tpm_mat)),
  all(is.finite(tpm_mat)),
  all(tpm_mat >= 0),
  setequal(colnames(tpm_mat), metadata$sample_id)
)

# Match metadata sample order
tpm_mat <- tpm_mat[, metadata$sample_id, drop = FALSE]
stopifnot(identical(colnames(tpm_mat), metadata$sample_id))

# Read LM22 and inspect the gene overlap
lm22 <- read.delim(
  sig_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

lm22 <- as.matrix(lm22)
storage.mode(lm22) <- "numeric"

stopifnot(
  nrow(lm22) > 0,
  ncol(lm22) == 22,
  !anyDuplicated(rownames(lm22)),
  all(is.finite(lm22))
)

shared_genes <- intersect(rownames(lm22), rownames(tpm_mat))

cat("TPM matrix:", nrow(tpm_mat), "genes x",
    ncol(tpm_mat), "samples\n")
cat("LM22 signature:", nrow(lm22), "genes x",
    ncol(lm22), "cell types\n")
cat("Shared genes:", length(shared_genes), "\n")

# Save mixture in CIBERSORT's required tab-delimited format
write.table(
  data.frame(GeneSymbol = rownames(tpm_mat), tpm_mat,
             check.names = FALSE),
  file = mixture_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Mixture file saved to:\n", mixture_file, "\n")
