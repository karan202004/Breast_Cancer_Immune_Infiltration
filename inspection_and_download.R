query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
results <- getResults(query)
nrow(results)
table(results$sample_type)
head(results)
colnames(results)
head(results[, c("case_id", "sample_type")])
length(unique(results$cases.submitter_id))
length(unique(results$sample.submitter_id))
table(results$sample_type)
table(table(results$cases.submitter_id))
results[, c("cases.submitter_id", "sample.submitter_id", "sample_type")]
head(results[, c(
  "cases.submitter_id",
  "sample.submitter_id",
  "sample_type"
)], 20)
tumor_results <- results[
  results$sample_type == "Primary Tumor",
]
length(unique(tumor_results$cases.submitter_id))
table(table(tumor_results$cases.submitter_id))
nrow(tumor_results)
multi_tumor_patients <- names(
  table(tumor_results$cases.submitter_id)[
    table(tumor_results$cases.submitter_id) > 1
  ]
)

multi_tumor_patients
tumor_results[
  tumor_results$cases.submitter_id %in% multi_tumor_patients,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "sample_type",
    "file_id",
    "file_size"
  )
]
results[
  results$cases.submitter_id == "TCGA-A7-A13E",
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "file_id",
    "file_name",
    "analysis_id",
    "analysis_workflow_version",
    "created_datetime",
    "updated_datetime"
  )
]
results[
  results$cases.submitter_id == "TCGA-A7-A13E",
  c(
    "sample.submitter_id",
    "file_id",
    "analysis_id",
    "file_name"
  )
]
GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  file.id = c(
    "0a511373-8fd3-433a-a5e1-877530d6a239",
    "4b96ee43-874c-4fa6-b5ce-08caf37ef2e1"
  )
)
# 1. Deduplicate by sample ID
tumor_unique <- tumor_results[!duplicated(tumor_results$sample.submitter_id), ]

# 2. Check total unique sample count
nrow(tumor_unique)

# 3. Check total unique patient count
length(unique(tumor_unique$cases.submitter_id))

# 4. Check frequency distribution of samples per patient
table(table(tumor_unique$cases.submitter_id))

tumor_unique[
  tumor_unique$cases.submitter_id %in% multi_tumor_patients,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "sample_type"
  )
]


library(TCGAbiolinks)
packageVersion("TCGAbiolinks")

tumor_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
tumor_query_results <- getResults(tumor_query)
nrow(tumor_query_results)
length(unique(tumor_query_results$sample.submitter_id))
sum(duplicated(tumor_query_results$sample.submitter_id))

final_tumor_results <- tumor_query_results[
  !duplicated(tumor_query_results$sample.submitter_id),
]
nrow(final_tumor_results)
sum(final_tumor_results$file_size) / 1024^3
sum(tumor_query_results$file_size) / 1024^3
range(final_tumor_results$file_size) / 1024^2

args(GDCquery)
tumor_query_results <- getResults(tumor_query)

tumor_query_results <- tumor_query_results[
  !duplicated(tumor_query_results$sample.submitter_id),
]
batch1_barcodes <- tumor_query_results$sample.submitter_id[1:100]
length(batch1_barcodes)
batch1_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch1_barcodes
)
batch1_results <- getResults(batch1_query)

nrow(batch1_results)
length(unique(batch1_results$sample.submitter_id))

GDCdownload(batch1_query)


write.table(
  batch1_results$file_id,
  file = "batch1_file_ids.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
write.csv(
  batch1_results,
  "batch1_metadata.csv",
  row.names = FALSE
)




write.table(
  batch1_results[, c("id", "file_size", "md5sum", "file_name")],
  "batch1_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

data_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data"

files <- list.files(
  data_dir,
  recursive = TRUE,
  full.names = TRUE
)

length(files)
tsv_files <- files[grepl("\\.tsv$", files, ignore.case = TRUE)]

length(tsv_files)
head(basename(tsv_files), 10)
readLines(tsv_files[1], n = 10)

library(TCGAbiolinks)
packageVersion("TCGAbiolinks")
batch2_barcodes <- tumor_query_results$sample.submitter_id[101:200]

length(batch2_barcodes)
batch2_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch2_barcodes
)
batch2_results <- getResults(batch2_query)


data_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data"

batch1_dir <- file.path(data_dir, "Batch1")

dir.create(batch1_dir, showWarnings = FALSE)
list.files(data_dir)

#batch1 folder creation 
data_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data"
batch1_dir <- file.path(data_dir, "Batch1")

# Read Batch 1 file IDs
batch1_ids <- readLines(
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch1_file_ids.txt"
)

# Remove empty lines
batch1_ids <- batch1_ids[nzchar(batch1_ids)]

# Move Batch 1 folders
for (id in batch1_ids) {
  source <- file.path(data_dir, id)
  destination <- file.path(batch1_dir, id)
  
  if (dir.exists(source)) {
    file.rename(source, destination)
  }
}
length(list.dirs(batch1_dir, recursive = FALSE, full.names = FALSE))
length(list.files(
  batch1_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))

data_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data"

batch2_dir <- file.path(data_dir, "Batch2")

dir.create(batch2_dir, showWarnings = FALSE)

write.table(
  batch2_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch2_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch2_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch2"

length(list.dirs(
  batch2_dir,
  recursive = FALSE,
  full.names = FALSE
))
length(list.files(
  batch2_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))


batch3_barcodes <- tumor_query_results$sample.submitter_id[201:300]

length(batch3_barcodes)
batch3_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch3_barcodes
)
batch3_results <- getResults(batch3_query)

nrow(batch3_results)
length(unique(batch3_results$sample.submitter_id))



data_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data"

batch3_dir <- file.path(data_dir, "Batch3")

dir.create(batch3_dir, showWarnings = FALSE)
write.table(
  batch3_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch3_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch3_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch3"

length(list.dirs(
  batch3_dir,
  recursive = FALSE,
  full.names = FALSE
))
length(list.files(
  batch3_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))

batch4_barcodes <- tumor_query_results$sample.submitter_id[301:400]

length(batch4_barcodes)
batch4_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch4_barcodes
)
batch4_results <- getResults(batch4_query)

nrow(batch4_results)
length(unique(batch4_results$sample.submitter_id))

batch4_dir <- file.path(data_dir, "Batch4")
dir.create(batch4_dir, showWarnings = FALSE)
write.table(
  batch4_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch4_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch4_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch4"

length(list.dirs(
  batch4_dir,
  recursive = FALSE,
  full.names = FALSE
))
length(list.files(
  batch4_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))


batch5_barcodes <- tumor_query_results$sample.submitter_id[401:500]

length(batch5_barcodes)
batch5_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch5_barcodes
)
batch5_results <- getResults(batch5_query)

nrow(batch5_results)
length(unique(batch5_results$sample.submitter_id))

batch5_dir <- file.path(data_dir, "Batch5")
dir.create(batch5_dir, showWarnings = FALSE)
write.table(
  batch5_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch5_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch5_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch5"

length(list.dirs(
  batch5_dir,
  recursive = FALSE,
  full.names = FALSE
))
length(list.files(
  batch5_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))

batch6_barcodes <- tumor_query_results$sample.submitter_id[501:600]

length(batch6_barcodes)
batch6_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch6_barcodes
)
batch6_results <- getResults(batch6_query)

nrow(batch6_results)
length(unique(batch6_results$sample.submitter_id))

batch6_dir <- file.path(data_dir, "Batch6")
dir.create(batch6_dir, showWarnings = FALSE)
write.table(
  batch6_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch6_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch6_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch6"

length(list.dirs(batch6_dir, recursive = FALSE, full.names = FALSE))
length(list.files(
  batch6_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))



batch7_barcodes <- tumor_query_results$sample.submitter_id[601:700]

length(batch7_barcodes)
batch7_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch7_barcodes
)

batch7_results <- getResults(batch7_query)

nrow(batch7_results)
length(unique(batch7_results$sample.submitter_id))


batch7_dir <- file.path(data_dir, "Batch7")
dir.create(batch7_dir, showWarnings = FALSE)
write.table(
  batch7_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch7_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch7_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch7"

length(list.dirs(batch7_dir, recursive = FALSE, full.names = FALSE))
length(list.files(
  batch7_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))

batch8_barcodes <- tumor_query_results$sample.submitter_id[701:800]

length(batch8_barcodes)
batch8_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch8_barcodes
)

batch8_results <- getResults(batch8_query)

nrow(batch8_results)
length(unique(batch8_results$sample.submitter_id))

batch8_dir <- file.path(data_dir, "Batch8")
dir.create(batch8_dir, showWarnings = FALSE)
write.table(
  batch8_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch8_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch8_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch8"

length(list.dirs(batch8_dir, recursive = FALSE, full.names = FALSE))
length(list.files(
  batch8_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))



batch9_barcodes <- tumor_query_results$sample.submitter_id[801:900]

length(batch9_barcodes)
batch9_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch9_barcodes
)

batch9_results <- getResults(batch9_query)

nrow(batch9_results)
length(unique(batch9_results$sample.submitter_id))

batch9_dir <- file.path(data_dir, "Batch9")
dir.create(batch9_dir, showWarnings = FALSE)
write.table(
  batch9_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch9_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


batch9_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch9"

length(list.dirs(batch9_dir, recursive = FALSE, full.names = FALSE))
length(list.files(
  batch9_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))


batch9_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch9"

length(list.dirs(batch9_dir, recursive = FALSE, full.names = FALSE))
length(list.files(
  batch9_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))


library(TCGAbiolinks)
packageVersion("TCGAbiolinks")
batch10_barcodes <- tumor_query_results$sample.submitter_id[901:1000]

length(batch10_barcodes)
batch10_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch10_barcodes
)

batch10_results <- getResults(batch10_query)

nrow(batch10_results)
length(unique(batch10_results$sample.submitter_id))


batch10_dir <- file.path(data_dir, "Batch10")
dir.create(batch10_dir, showWarnings = FALSE)
write.table(
  batch10_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch10_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

batch10_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch10"

length(list.dirs(batch10_dir, recursive = FALSE, full.names = FALSE))
length(list.files(
  batch10_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))

batch11_barcodes <- tumor_query_results$sample.submitter_id[1001:1106]

length(batch11_barcodes)
batch11_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = batch11_barcodes
)

batch11_results <- getResults(batch11_query)

nrow(batch11_results)
length(unique(batch11_results$sample.submitter_id))

batch11_dir <- file.path(data_dir, "Batch11")
dir.create(batch11_dir, showWarnings = FALSE)
write.table(
  batch11_results[, c("id", "file_size", "md5sum", "file_name")],
  "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/batch11_manifest.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)



batch11_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration/data/Batch11"

length(list.dirs(
  batch11_dir,
  recursive = FALSE,
  full.names = FALSE
))
length(list.files(
  batch11_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
))
