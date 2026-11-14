project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
raw_dir <- file.path(project_dir, "metadata", "raw_clinical")
master  <- readRDS(file.path(project_dir, "metadata", "clinical_master_table.rds"))

## 1. "Do not use" samples: which platform, which vial, and why
qa <- read.delim(file.path(raw_dir, "PanCanAtlas_merged_sample_quality_annotations.tsv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
col_of <- function(p) grep(p, names(qa), ignore.case = TRUE, value = TRUE)[1]
c_type <- col_of("cancer"); c_plat <- col_of("platform"); c_id <- col_of("aliquot"); c_dnu <- col_of("do_not_use")
cat("Annotation columns:\n"); print(names(qa))

qa <- qa[qa[[c_type]] == "BRCA", ]
qa$sample15 <- substr(qa[[c_id]], 1, 15)
qa$our_sample <- master$sample_id[match(qa$sample15, substr(master$sample_id, 1, 15))]
dnu <- qa[!is.na(qa$our_sample) & toupper(qa[[c_dnu]]) == "TRUE", ]
dnu$same_vial <- substr(dnu[[c_id]], 1, 16) == dnu$our_sample
dnu$in_analysis_set <- master$in_analysis_set[match(dnu$our_sample, master$sample_id)]

cat("\nDo_not_use flags by platform and whether it is the same vial we used:\n")
print(table(platform = dnu[[c_plat]], same_vial = dnu$same_vial))

mrna <- dnu[grepl("RNASeq", dnu[[c_plat]]) & !grepl("miRNA", dnu[[c_plat]]), ]
cat("\nmRNA-seq flags (", nrow(mrna), "rows ). All annotation columns:\n")
print(mrna, row.names = FALSE)

## 2. HER2: how status was decided, and ERBB2 expression
cat("\nHER2 IHC vs FISH:\n")
print(table(IHC = master$her2_ihc, FISH = master$her2_fish, useNA = "ifany"))

tpm <- read.delim(file.path(project_dir, "data", "CIBERSORT", "BRCA_TPM_1106.txt"),
                  check.names = FALSE, stringsAsFactors = FALSE)
erbb2_row <- tpm[tpm[[1]] == "ERBB2", -1]
names(erbb2_row) <- gsub(".", "-", names(erbb2_row), fixed = TRUE)
master$ERBB2 <- log2(as.numeric(erbb2_row[master$sample_id]) + 1)

cat("\nERBB2 log2(TPM+1) by HER2 status:\n")
print(do.call(rbind, tapply(master$ERBB2, master$her2_status, quantile, probs = c(.1, .25, .5, .75, .9))))

neg95 <- quantile(master$ERBB2[master$her2_status %in% "Negative"], 0.95, na.rm = TRUE)
pos <- master[master$her2_status %in% "Positive", ]
cat("\nHER2-positive tumours with ERBB2 above the 95th percentile of HER2-negative:",
    sum(pos$ERBB2 > neg95), "of", nrow(pos), "\n")
cat("\nBy PAM50 (HER2-positive only): n, and n with high ERBB2\n")
print(cbind(n = table(pos$pam50), high_ERBB2 = tapply(pos$ERBB2 > neg95, pos$pam50, sum)))
cat("\nHow Luminal A HER2-positive calls were made (IHC x FISH):\n")
la <- pos[pos$pam50 %in% "Luminal A", ]
print(table(IHC = la$her2_ihc, FISH = la$her2_fish, useNA = "ifany"))

## 3. Lymph node column name in the clinical supplement
bt <- readRDS(file.path(raw_dir, "gdc_biotab_clinical_raw.rds"))
cat("\nLymph-node columns in the clinical supplement:\n")
print(grep("lymph|node", names(bt[["clinical_patient_brca"]]), value = TRUE, ignore.case = TRUE))