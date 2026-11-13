library(TCGAbiolinks)

# Your original tumor query results
# (tumor_query_results should already exist)

# Check the available sample types
table(tumor_query_results$sample_type)
# See all sample types
unique(tumor_query_results$sample_type)

# Count each sample type
sample_type_counts <- as.data.frame(
  table(tumor_query_results$sample_type)
)

sample_type_counts

tumor_samples <- tumor_query_results[
  tumor_query_results$sample_type == "Primary Tumor",
]

normal_samples <- tumor_query_results[
  tumor_query_results$sample_type == "Solid Tissue Normal",
]

nrow(tumor_samples)
nrow(normal_samples)
setdiff(
  unique(tumor_query_results$sample_type),
  c("Primary Tumor", "Solid Tissue Normal")
)
table(tumor_query_results$sample_type)


# Find samples that have multiple files
duplicate_samples <- tumor_query_results[
  duplicated(tumor_query_results$sample.submitter_id) |
    duplicated(tumor_query_results$sample.submitter_id, fromLast = TRUE),
]

# Show the duplicate samples
duplicate_samples[, c(
  "sample.submitter_id",
  "id",
  "file_name",
  "file_size",
  "md5sum"
)]
# Count how many files each sample has
duplicate_counts <- table(tumor_query_results$sample.submitter_id)

duplicate_counts[duplicate_counts > 1]

# Compare total file records with unique samples
nrow(tumor_query_results)
length(unique(tumor_query_results$sample.submitter_id))

# Look at the important columns
names(tumor_query_results)

# Show how many records exist for each sample.submitter_id
head(
  sort(
    table(tumor_query_results$sample.submitter_id),
    decreasing = TRUE
  ),
  20
)
# Look at all file-level identifiers
tumor_query_results[, c(
  "id",
  "sample.submitter_id",
  "file_name",
  "file_size",
  "md5sum"
)][1:20, ]

library(TCGAbiolinks)

# ============================================================
# 1. Main directories
# ============================================================

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"

data_dir <- file.path(project_dir, "data")

# ============================================================
# 2. List all downloaded STAR-Counts files
# ============================================================

all_tsv_files <- list.files(
  data_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

length(all_tsv_files)

batch_counts <- sapply(
  1:11,
  function(i) {
    batch_dir <- file.path(data_dir, paste0("Batch", i))
    
    length(list.files(
      batch_dir,
      pattern = "\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    ))
  }
)

batch_counts
all_downloaded_ids <- sub(
  ".*/([^/\\\\]+)/[^/\\\\]+$",
  "\\1",
  gsub("\\\\", "/", all_tsv_files)
)

length(all_downloaded_ids)
expected_ids <- tumor_query_results$id

extra_ids <- setdiff(
  all_downloaded_ids,
  expected_ids
)

length(extra_ids)
extra_ids
missing_ids <- setdiff(
  expected_ids,
  all_downloaded_ids
)

length(missing_ids)
missing_ids
cat("Expected samples:", length(unique(tumor_query_results$sample.submitter_id)), "\n")
cat("Expected files:", nrow(tumor_query_results), "\n")
cat("Downloaded TSV files:", length(all_tsv_files), "\n")
cat("Extra downloaded files:", length(extra_ids), "\n")
cat("Missing expected files:", length(missing_ids), "\n")

extra_info <- tumor_query_results[
  tumor_query_results$id %in% extra_ids,
]

extra_info
extra_files <- all_tsv_files[
  grepl(
    paste(extra_ids, collapse = "|"),
    all_tsv_files
  )
]

extra_files
extra_files

library(jsonlite)

extra_ids
extra_info_list <- lapply(
  extra_ids,
  function(x) {
    url <- paste0(
      "https://api.gdc.cancer.gov/files/",
      x
    )
    
    fromJSON(url)$data
  }
)
extra_info <- do.call(
  rbind,
  lapply(extra_info_list, as.data.frame)
)

extra_info
extra_info <- do.call(
  rbind,
  lapply(extra_info_list, as.data.frame)
)

extra_info
extra_info[, c(
  "file_id",
  "file_name",
  "file_size",
  "cases"
)]
names(extra_info)


library(jsonlite)

test_id <- extra_ids[1]

url <- paste0(
  "https://api.gdc.cancer.gov/files/",
  test_id,
  "?expand=cases.samples"
)

test_result <- fromJSON(url, flatten = TRUE)$data

names(test_result)
str(test_result, max.level = 3)

samples_test <- test_result$cases[[1]]$samples[[1]]

names(samples_test)
str(samples_test, max.level = 2)
library(jsonlite)

extra_sample_info <- lapply(
  extra_ids,
  function(file_id) {
    
    url <- paste0(
      "https://api.gdc.cancer.gov/files/",
      file_id,
      "?expand=cases.samples"
    )
    
    result <- fromJSON(url)$data
    
    # Extract case
    case <- result$cases[[1]]
    
    # Extract sample
    sample <- case$samples[[1]]
    
    data.frame(
      file_id = file_id,
      file_name = result$file_name,
      sample_submitter_id = sample$submitter_id,
      sample_type = sample$sample_type,
      case_submitter_id = case$submitter_id,
      stringsAsFactors = FALSE
    )
  }
)

extra_sample_info <- do.call(
  rbind,
  extra_sample_info
)

extra_sample_info
table(extra_sample_info$sample_type)
extra_sample_info[, c(
  "file_id",
  "file_name",
  "sample_submitter_id",
  "sample_type",
  "case_submitter_id"
)]

extra_file_info <- lapply(
  extra_files,
  function(f) {
    
    # Read first 10 lines
    lines <- readLines(f, n = 10)
    
    data.frame(
      file = f,
      header = paste(lines, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
)

extra_file_info <- do.call(
  rbind,
  extra_file_info
)

extra_file_info
# ============================================================
# CREATE CLEAN DATA DIRECTORY
# ============================================================

clean_dir <- file.path(
  data_dir,
  "Primary_Tumor_1106"
)

dir.create(
  clean_dir,
  showWarnings = FALSE
)
# ============================================================
# CREATE FILE-ID → DOWNLOADED FILE MAPPING
# ============================================================

downloaded_map <- data.frame(
  file_id = all_downloaded_ids,
  file_path = all_tsv_files,
  stringsAsFactors = FALSE
)

head(downloaded_map)
# ============================================================
# KEEP ONLY EXPECTED 1106 FILES
# ============================================================

expected_files <- downloaded_map[
  downloaded_map$file_id %in% tumor_query_results$id,
]

nrow(expected_files)
length(unique(expected_files$file_id))
# ============================================================
# ADD SAMPLE INFORMATION
# ============================================================

clean_metadata <- tumor_query_results[
  match(
    expected_files$file_id,
    tumor_query_results$id
  ),
]

nrow(clean_metadata)
table(clean_metadata$sample_type)
length(unique(clean_metadata$sample.submitter_id))

library(TCGAbiolinks)
library(jsonlite)
# ============================================================
# SAVE CLEAN METADATA FOR THE 1106 PRIMARY TUMOR SAMPLES
# ============================================================

clean_metadata <- data.frame(
  sample_id = clean_metadata$sample.submitter_id,
  file_id = clean_metadata$id,
  file_path = expected_files$file_path,
  sample_type = clean_metadata$sample_type,
  case_id = clean_metadata$cases.submitter_id,
  stringsAsFactors = FALSE
)

# Check
head(clean_metadata)

dim(clean_metadata)

table(clean_metadata$sample_type)

length(unique(clean_metadata$sample_id))
metadata_file <- file.path(
  project_dir,
  "metadata",
  "primary_tumor_1106_metadata.csv"
)

dir.create(
  dirname(metadata_file),
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  clean_metadata,
  metadata_file,
  row.names = FALSE
)

metadata_file

# ============================================================
# TEST ONE STAR-COUNTS FILE
# ============================================================

test_file <- clean_metadata$file_path[1]

test_lines <- readLines(test_file, n = 10)

cat(test_lines, sep = "\n")
test_counts <- read.delim(
  test_file,
  comment.char = "#",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(test_counts)

head(test_counts)
colnames(test_counts)
nrow(test_counts)

# ============================================================
# TEST GENE EXTRACTION FROM ONE STAR-COUNTS FILE
# ============================================================

test_gene_data <- test_counts[
  !grepl("^N_", test_counts$gene_id),
]

# Check dimensions
dim(test_gene_data)

# Check first genes
head(
  test_gene_data[, c(
    "gene_id",
    "gene_name",
    "gene_type",
    "unstranded"
  )]
)

# Check number of unique genes
length(unique(test_gene_data$gene_id))

# Check missing gene names
sum(
  is.na(test_gene_data$gene_name) |
    test_gene_data$gene_name == ""
)

# Check duplicate gene IDs
sum(duplicated(test_gene_data$gene_id))

# Check gene types
table(test_gene_data$gene_type)


# ============================================================
# SAFE CLEANUP & STRICT AUDIT OF GDC DATA DIRECTORY
# ============================================================

library(dplyr)

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
metadata_path <- file.path(project_dir, "metadata", "primary_tumor_1106_metadata.csv")
data_dir <- file.path(project_dir, "data")

# 1. Load authoritative clean manifest (1,106 samples)
clean_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
valid_file_ids <- clean_metadata$file_id

cat("============================================================\n")
cat("PRE-CLEANUP AUDIT\n")
cat("============================================================\n")
cat("Target manifest sample count (authoritative):", length(valid_file_ids), "\n")

# 2. Identify all downloaded TSV files on disk
all_tsv_files <- list.files(
  path = data_dir, 
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", 
  recursive = TRUE, 
  full.names = TRUE
)

# Extract parent folder UUIDs (file_ids) from file paths
disk_file_ids <- basename(dirname(all_tsv_files))

cat("Total TSV files on disk before cleanup:      ", length(all_tsv_files), "\n")

# Check for duplicate file IDs on disk
duplicate_disk_ids <- unique(disk_file_ids[duplicated(disk_file_ids)])
cat("Duplicate file IDs on disk:                 ", length(duplicate_disk_ids), "\n")

if (length(duplicate_disk_ids) > 0) {
  cat("  -> Duplicate UUIDs detected:\n")
  print(duplicate_disk_ids)
}

# 3. Derive extra file IDs by set difference
extra_file_ids <- setdiff(disk_file_ids, valid_file_ids)
cat("Extra file IDs identified for removal:       ", length(extra_file_ids), "\n")

if (length(extra_file_ids) > 0) {
  cat("\nTargeted UUIDs for deletion:\n")
  print(extra_file_ids)
  
  cat("\nDeleting extra file directories...\n")
  for (extra_id in extra_file_ids) {
    target_tsv <- all_tsv_files[basename(dirname(all_tsv_files)) == extra_id]
    target_dir <- dirname(target_tsv)
    
    if (dir.exists(target_dir)) {
      unlink(target_dir, recursive = TRUE, force = TRUE)
      cat("  - Successfully removed:", target_dir, "\n")
    }
  }
} else {
  cat("\nNo extra files found to delete.\n")
}

# ============================================================
# POST-CLEANUP STRICT RE-AUDIT
# ============================================================

remaining_tsv_files <- list.files(
  path = data_dir, 
  pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", 
  recursive = TRUE, 
  full.names = TRUE
)

remaining_file_ids <- basename(dirname(remaining_tsv_files))

cat("\n============================================================\n")
cat("POST-CLEANUP VERIFICATION\n")
cat("============================================================\n")
cat("Target manifest files (expected): ", length(valid_file_ids), "\n")
cat("Downloaded TSV files remaining:   ", length(remaining_tsv_files), "\n")
cat("Unique file IDs on disk:         ", length(unique(remaining_file_ids)), "\n")
cat("Missing expected files:           ", length(setdiff(valid_file_ids, remaining_file_ids)), "\n")
cat("Extra unmapped files:             ", length(setdiff(remaining_file_ids, valid_file_ids)), "\n")
cat("============================================================\n")

if (
  length(remaining_tsv_files) == 1106 &&
  length(unique(remaining_file_ids)) == 1106 &&
  length(setdiff(valid_file_ids, remaining_file_ids)) == 0 &&
  length(setdiff(remaining_file_ids, valid_file_ids)) == 0
) {
  cat("\nSUCCESS: Data directory is 100% clean and matches clean_metadata exactly!\n")
} else {
  cat("\nWARNING: File count mismatch remaining. Re-examine audit log.\n")
}