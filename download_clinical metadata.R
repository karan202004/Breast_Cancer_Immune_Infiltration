#Download clinical + PAM50 data for 1,106 samples

suppressPackageStartupMessages(library(TCGAbiolinks))

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
meta_dir    <- file.path(project_dir, "metadata")
raw_dir     <- file.path(meta_dir, "raw_clinical")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)


samples <- read.csv(file.path(meta_dir, "primary_tumor_1106_metadata.csv"),
                    stringsAsFactors = FALSE)
if (nrow(samples) != 1106 || anyDuplicated(samples$sample_id))
  stop("Sample metadata is not 1,106 unique samples.")

samples$patient <- substr(samples$sample_id, 1, 12)
if (!all(samples$patient == samples$case_id))
  warning("case_id does not always equal the first 12 characters of sample_id; using sample_id.")

n_pat <- length(unique(samples$patient))
multi <- unique(samples$patient[duplicated(samples$patient)])
cat("Samples:", nrow(samples), "| unique patients:", n_pat,
    "| patients with >1 primary tumour sample:", length(multi), "\n\n")

sources <- list()
pick <- function(df, field, exact, pattern = NULL) {
  hit <- intersect(exact, names(df))
  if (!length(hit) && !is.null(pattern))
    hit <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (!length(hit)) {
    sources[[field]] <<- "NOT FOUND"
    return(rep(NA, nrow(df)))
  }
  sources[[field]] <<- hit[1]
  df[[hit[1]]]
}

#Clinical data from GDC
cat("Downloading GDC clinical data...\n")
clin_raw <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
saveRDS(clin_raw, file.path(raw_dir, "gdc_clinical_raw.rds"))
cat("GDC clinical:", nrow(clin_raw), "rows x", ncol(clin_raw), "columns\n")

clin <- data.frame(
  patient            = pick(clin_raw, "patient", "submitter_id"),
  gender             = pick(clin_raw, "gender", c("gender", "demographic.gender")),
  race               = pick(clin_raw, "race", c("race", "demographic.race")),
  age_at_index       = pick(clin_raw, "age_at_index", c("age_at_index", "demographic.age_at_index")),
  age_at_diag_days   = pick(clin_raw, "age_at_diagnosis", c("age_at_diagnosis", "diagnoses.age_at_diagnosis")),
  vital_status       = pick(clin_raw, "vital_status", c("vital_status", "demographic.vital_status")),
  days_to_death      = pick(clin_raw, "days_to_death", c("days_to_death", "demographic.days_to_death")),
  days_to_last_fu    = pick(clin_raw, "days_to_last_follow_up", c("days_to_last_follow_up", "diagnoses.days_to_last_follow_up")),
  stage_raw          = pick(clin_raw, "ajcc_pathologic_stage", c("ajcc_pathologic_stage", "diagnoses.ajcc_pathologic_stage"), "pathologic_stage$"),
  pT                 = pick(clin_raw, "ajcc_pathologic_t", c("ajcc_pathologic_t", "diagnoses.ajcc_pathologic_t"), "pathologic_t$"),
  pN                 = pick(clin_raw, "ajcc_pathologic_n", c("ajcc_pathologic_n", "diagnoses.ajcc_pathologic_n"), "pathologic_n$"),
  pM                 = pick(clin_raw, "ajcc_pathologic_m", c("ajcc_pathologic_m", "diagnoses.ajcc_pathologic_m"), "pathologic_m$"),
  histology          = pick(clin_raw, "primary_diagnosis", c("primary_diagnosis", "diagnoses.primary_diagnosis")),
  stringsAsFactors = FALSE
)

if (all(is.na(clin$patient))) stop("Could not find patient IDs (submitter_id) in GDC clinical data.")
if (anyDuplicated(clin$patient)) {
  cat("Note:", sum(duplicated(clin$patient)), "duplicated patient rows in GDC clinical; keeping first.\n")
  clin <- clin[!duplicated(clin$patient), ]
}

#Derived variables
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
clin$age <- num(clin$age_at_index)
fill_age <- is.na(clin$age) & !is.na(num(clin$age_at_diag_days))
clin$age[fill_age] <- floor(num(clin$age_at_diag_days[fill_age]) / 365.25)

st <- toupper(trimws(as.character(clin$stage_raw)))
clin$stage <- ifelse(grepl("^STAGE (IV|III|II|I)", st),
                     sub("^STAGE (IV|III|II|I).*$", "\\1", st), NA)
clin$stage <- factor(clin$stage, levels = c("I", "II", "III", "IV"))

dead <- tolower(clin$vital_status) == "dead"
clin$os_event <- ifelse(is.na(clin$vital_status), NA, as.integer(dead))
clin$os_days  <- ifelse(dead, num(clin$days_to_death), num(clin$days_to_last_fu))

##PAM50 subtype
cat("\nDownloading PAM50 subtypes...\n")
pam <- NULL
sub_raw <- tryCatch(TCGAquery_subtype(tumor = "brca"), error = function(e) NULL)
if (!is.null(sub_raw)) {
  saveRDS(sub_raw, file.path(raw_dir, "tcga_brca_subtype_raw.rds"))
  pcol <- grep("PAM50", names(sub_raw), value = TRUE)[1]
  if (!is.na(pcol) && "patient" %in% names(sub_raw)) {
    pam <- data.frame(patient = as.character(sub_raw$patient),
                      pam50_raw = as.character(sub_raw[[pcol]]), stringsAsFactors = FALSE)
    sources$pam50 <- paste0("TCGAquery_subtype: ", pcol)
  }
}
if (is.null(pam)) {   # fallback: PanCancer Atlas mRNA subtypes
  pc <- tryCatch(PanCancerAtlas_subtypes(), error = function(e) NULL)
  if (!is.null(pc)) {
    saveRDS(pc, file.path(raw_dir, "pancanatlas_subtypes_raw.rds"))
    pc <- pc[pc$cancer.type == "BRCA", ]
    pam <- data.frame(patient = substr(pc$pan.samplesID, 1, 12),
                      pam50_raw = as.character(pc$Subtype_mRNA), stringsAsFactors = FALSE)
    sources$pam50 <- "PanCancerAtlas_subtypes: Subtype_mRNA"
  }
}
if (is.null(pam)) stop("Could not download PAM50 subtypes from either source.")

pam <- pam[!is.na(pam$patient) & !duplicated(pam$patient), ]
lab <- c(LumA = "Luminal A", LumB = "Luminal B", Her2 = "HER2-enriched",
         Basal = "Basal-like", Normal = "Normal-like")
pam$pam50 <- ifelse(pam$pam50_raw %in% names(lab), lab[pam$pam50_raw], NA)
odd <- setdiff(unique(pam$pam50_raw[is.na(pam$pam50)]), c(NA, "", "NA"))
if (length(odd)) cat("PAM50 labels set to NA (unrecognised):", paste(odd, collapse = ", "), "\n")

##Merge onto the 1,106 samples (metadata order)
out <- merge(samples[, c("sample_id", "patient")], clin, by = "patient", all.x = TRUE, sort = FALSE)
out <- merge(out, pam[, c("patient", "pam50")], by = "patient", all.x = TRUE, sort = FALSE)
out <- out[match(samples$sample_id, out$sample_id), ]
stopifnot(identical(out$sample_id, samples$sample_id))

out$pam50 <- factor(out$pam50, levels = unname(lab))
out$multi_sample_patient <- out$patient %in% multi

keep <- c("sample_id", "patient", "multi_sample_patient", "pam50", "age", "gender", "race",
          "stage", "stage_raw", "pT", "pN", "pM", "histology",
          "vital_status", "os_event", "os_days")
out <- out[, keep]
rownames(out) <- NULL

write.csv(out, file.path(meta_dir, "clinical_phase4_1106.csv"), row.names = FALSE)
saveRDS(out, file.path(meta_dir, "clinical_phase4_1106.rds"))

cat("Source columns used:\n")
print(data.frame(field = names(sources), column = unlist(sources)), row.names = FALSE)

cat("\nSamples with clinical record:", sum(!is.na(out$vital_status) | !is.na(out$age)), "of", nrow(out), "\n")
cat("\nPAM50:\n");        print(table(out$pam50, useNA = "always"))
cat("\nStage:\n");        print(table(out$stage, useNA = "always"))
cat("\nGender:\n");       print(table(out$gender, useNA = "always"))
cat("\nHistology (top 5):\n"); print(head(sort(table(out$histology, useNA = "always"), decreasing = TRUE), 5))
cat("\nAge:\n");          print(summary(out$age))
cat("\nVital status:\n"); print(table(out$vital_status, useNA = "always"))
cat("Samples with usable OS time:", sum(!is.na(out$os_days) & out$os_days >= 0), "\n")
cat("\nSaved:", file.path(meta_dir, "clinical_phase4_1106.csv"), "\n")
