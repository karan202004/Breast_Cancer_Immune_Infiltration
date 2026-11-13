# ==============================================================================
# SCRIPT 01: BUILD RAW COUNT EXPRESSION MATRIX & GENE ANNOTATION (FORTIFIED)
# Project: TCGA-BRCA Immune Infiltration Analysis
# Author: Karan S.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & DEPENDENCIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

# Define Paths
project_dir   <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
metadata_path <- file.path(project_dir, "metadata", "primary_tumor_1106_metadata.csv")
data_dir      <- file.path(project_dir, "data")
results_dir   <- file.path(project_dir, "results", "expression_matrix")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("[INIT] Created directory:", results_dir, "\n")
}

cat("==================================================================\n")
cat("STEP 01: FORTIFIED RAW EXPRESSION MATRIX CONSTRUCTION (PAR AGGREGATED)\n")
cat("==================================================================\n\n")

# ------------------------------------------------------------------------------
# 2. LOAD & VALIDATE METADATA
# ------------------------------------------------------------------------------
cat("[1/8] Auditing metadata...\n")
if (!file.exists(metadata_path)) {
  stop("FATAL ERROR: Metadata file not found at: ", metadata_path)
}

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
n_samples <- nrow(metadata)

stopifnot(n_samples == 1106)
stopifnot(length(unique(metadata$sample_id)) == 1106)
stopifnot(length(unique(metadata$file_id)) == 1106)

cat("  -> Metadata validated:", n_samples, "unique primary tumor samples.\n")

# ------------------------------------------------------------------------------
# 3. LOCATE & AUDIT DISK FILES
# ------------------------------------------------------------------------------
cat("[2/8] Locating & auditing TSV files on disk...\n")

all_disk_tsvs <- list.files(
  path = data_dir,
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(all_disk_tsvs) != 1106) {
  stop("FATAL ERROR: Expected 1106 TSV files on disk, but found ", length(all_disk_tsvs))
}

disk_file_ids <- basename(dirname(all_disk_tsvs))
if (anyDuplicated(disk_file_ids) > 0) {
  stop("FATAL ERROR: Duplicate file IDs (UUIDs) detected on disk!")
}

disk_lookup <- setNames(all_disk_tsvs, disk_file_ids)
file_paths  <- disk_lookup[metadata$file_id]

if (any(is.na(file_paths)) || any(!file.exists(file_paths))) {
  stop("FATAL ERROR: One or more file_ids in metadata could not be mapped to disk paths!")
}

cat("  -> Exactly 1,106 TSVs successfully mapped to metadata sample IDs.\n")

# ------------------------------------------------------------------------------
# 4. TEST FILE 1 PARSING & BUILD AGGREGATED GENE ANNOTATION
# ------------------------------------------------------------------------------
cat("[3/8] Building aggregated gene annotation template (handling PAR Y genes)...\n")

first_file <- fread(
  file_paths[1], 
  skip = "gene_id", 
  header = TRUE, 
  showProgress = FALSE
)

# Remove GDC summary rows (N_unmapped, etc.)
first_file <- first_file[!grepl("^N_", gene_id)]

# Base Ensembl IDs (stripped of version and _PAR_Y suffix)
base_gene_ids <- gsub("\\..*$", "", first_file$gene_id)
first_file$base_gene_id <- base_gene_ids

# Unique base gene annotation table (preserving first occurrence metadata)
unique_annotation <- first_file[!duplicated(base_gene_id), .(
  base_gene_id,
  gene_name,
  gene_type
)]

unique_gene_ids <- unique_annotation$base_gene_id
n_unique_genes  <- length(unique_gene_ids)

cat("  -> Filtered summary rows and aggregated PAR genes.\n")
cat("  -> Total unique features post-PAR aggregation:", n_unique_genes, "genes.\n")

# ------------------------------------------------------------------------------
# 5. ALLOCATE MATRIX & ITERATE WITH AGGREGATION
# ------------------------------------------------------------------------------
cat("[4/8] Processing 1,106 samples with on-the-fly PAR summation...\n")

raw_matrix <- matrix(
  0L, 
  nrow = n_unique_genes, 
  ncol = n_samples,
  dimnames = list(unique_gene_ids, metadata$sample_id)
)

pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

for (i in seq_len(n_samples)) {
  dt <- fread(
    file_paths[i], 
    skip = "gene_id",
    select = c("gene_id", "unstranded"),
    showProgress = FALSE
  )
  
  dt <- dt[!grepl("^N_", gene_id)]
  dt[, base_id := gsub("\\..*$", "", gene_id)]
  
  # Sum counts by base_gene_id (combining PAR_Y rows)
  dt_agg <- dt[, .(count = sum(unstranded)), by = base_id]
  
  # Align to master gene order
  aligned_counts <- dt_agg$count[match(unique_gene_ids, dt_agg$base_id)]
  
  # Integrity checks
  stopifnot(!anyNA(aligned_counts))
  stopifnot(all(aligned_counts >= 0))
  
  raw_matrix[, i] <- as.integer(aligned_counts)
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n  -> Matrix construction and PAR aggregation complete.\n")

# ------------------------------------------------------------------------------
# 6. STRUCTURAL ASSERTIONS
# ------------------------------------------------------------------------------
cat("[5/8] Running final structural assertions on aggregated matrix...\n")

stopifnot(nrow(raw_matrix) == n_unique_genes)
stopifnot(ncol(raw_matrix) == n_samples)
stopifnot(anyNA(raw_matrix) == FALSE)
stopifnot(all(raw_matrix >= 0))
stopifnot(identical(colnames(raw_matrix), metadata$sample_id))
stopifnot(identical(rownames(raw_matrix), unique_gene_ids))
stopifnot(anyDuplicated(rownames(raw_matrix)) == 0)

cat("  -> Matrix passed all structural and integrity checks.\n")

# ------------------------------------------------------------------------------
# 7. EXPORT DATASETS
# ------------------------------------------------------------------------------
cat("[6/8] Writing output datasets to disk...\n")

gene_annotation <- as.data.frame(unique_annotation)
colnames(gene_annotation)[1] <- "gene_id"

saveRDS(raw_matrix, file = file.path(results_dir, "raw_counts_1106.rds"))
saveRDS(gene_annotation, file = file.path(results_dir, "gene_annotation_1106.rds"))

write.csv(raw_matrix, file = file.path(results_dir, "raw_counts_1106.csv"), quote = FALSE)
write.csv(gene_annotation, file = file.path(results_dir, "gene_annotation_1106.csv"), row.names = FALSE)

# Processing Log
log_path <- file.path(results_dir, "processing_summary_01.txt")
log_conn <- file(log_path, "w")
cat("TCGA-BRCA EXPRESSION MATRIX CONSTRUCTION SUMMARY\n", file = log_conn)
cat("================================================\n", file = log_conn)
cat("Execution Date:        ", as.character(Sys.time()), "\n", file = log_conn)
cat("Files processed:       ", n_samples, "\n", file = log_conn)
cat("Samples retained:      ", n_samples, "\n", file = log_conn)
cat("Total genes retained:  ", n_unique_genes, "\n", file = log_conn)
cat("PAR genes aggregated:  44 base IDs (88 rows -> 44 rows)\n", file = log_conn)
close(log_conn)

cat("[7/8] Processing log written to:", log_path, "\n")
cat("==================================================================\n")
cat("SUCCESS: Step 01 complete! 60,616 clean genes x 1,106 samples assembled.\n")
cat("==================================================================\n")