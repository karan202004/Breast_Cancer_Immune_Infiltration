project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"

cibersort_file <- file.path(
  project_dir,
  "data",
  "deconvolution_inputs",
  "cibersortx_input_linear_cpm.txt"
)

stopifnot(file.exists(cibersort_file))

x <- read.delim(
  cibersort_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Rows:", nrow(x), "\n")
cat("Columns:", ncol(x), "\n")
cat("First column:", colnames(x)[1], "\n")
cat("Samples:", ncol(x) - 1, "\n")
cat("Duplicate genes:", anyDuplicated(x[[1]]), "\n")
cat("NA:", anyNA(x), "\n")
