#Inspect duplicate base ENSEMBL IDs

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
results_dir <- file.path(project_dir, "results", "expression_matrix")

annotation_path <- file.path(results_dir, "gene_annotation_1106.rds")
matrix_path     <- file.path(results_dir, "raw_counts_1106.rds")

if (!file.exists(annotation_path) || !file.exists(matrix_path)) {
  stop("FATAL ERROR: Save files not found! Run Script 01 first.")
}

cat("loading saved RDS objects...\n")
gene_annotation <- readRDS(annotation_path)
raw_matrix      <- readRDS(matrix_path)

cat("Extracting duplicate base Ensembl IDs...\n")
dup_base_ids <- unique(
  gene_annotation$base_gene_id[
    duplicated(gene_annotation$base_gene_id)
  ]
)

cat("Number of unique duplicated base Ensembl IDs:", length(dup_base_ids), "\n")
cat("Total matrix rows involved in duplication:  ", sum(gene_annotation$base_gene_id %in% dup_base_ids), "\n\n")

# Filter duplicate records and sort by base_gene_id
duplicate_annotation <- gene_annotation[
  gene_annotation$base_gene_id %in% dup_base_ids, 
]
duplicate_annotation <- duplicate_annotation[
  order(duplicate_annotation$base_gene_id), 
]

cat("Detailed duplicate annotation table\n")
print(head(duplicate_annotation, 30)) # Prints first 30 rows

if (nrow(duplicate_annotation) > 30) {
  cat("\n... [top 30 rows out of", nrow(duplicate_annotation), "total duplicate entries]\n")
}

cat("\nChecking PAR / Y suffix presence...\n")
par_count <- sum(grepl("_PAR_Y", duplicate_annotation$gene_id_version))
cat("  -> Rows containing '_PAR_Y' tag:", par_count, "\n\n")

cat("Unique row names in matrix: ", length(unique(rownames(raw_matrix))), "\n")
cat("Matrix dimensions:          ", nrow(raw_matrix), "genes x", ncol(raw_matrix), "samples\n")

