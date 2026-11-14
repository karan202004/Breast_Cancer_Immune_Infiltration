# ============================================================
# CLINICAL AND MOLECULAR DATA + ANALYSIS SET
# TCGA-BRCA | 1,106 primary tumour RNA-seq samples
#
# Run order:  1. clinical_data_and_analysis_set.R   (this script)
#             2. phase4_immune_landscape.R
#             3. clinical_master_table_qc.R
#
# Replaces phase4_get_clinical_metadata.R, the sex_at_birth fix and the
# clinical part of clinical_data_matching_qc.R.
#
# Data sources:
#   - GDC clinical (TCGAbiolinks::GDCquery_clinic)
#   - PAM50: Berger et al. Cancer Cell 2018;33:690-705 (TCGAbiolinks::TCGAquery_subtype)
#   - GDC BCR Biotab clinical supplement: ER, PR, HER2, menopause
#   - TCGA-CDR, Liu et al. Cell 2018;173:400-416: curated survival endpoints
#   - PanCanAtlas (Thorsson et al. Immunity 2018;48:812-830): leukocyte fraction
#     (DNA methylation), ABSOLUTE purity, merged sample quality annotations
#
# Downloads are skipped when the file already exists in metadata/raw_clinical.
# ============================================================

suppressPackageStartupMessages(library(TCGAbiolinks))
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
options(timeout = 900)

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
meta_dir <- file.path(project_dir, "metadata")
raw_dir  <- file.path(meta_dir, "raw_clinical")
out_dir  <- file.path(project_dir, "results", "clinical_integration")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Helpers -------------------------------------------------------
checks <- list()
log_check <- function(link, check, result, status, detail = "") {
  checks[[length(checks) + 1]] <<- data.frame(
    link = link, check = check, result = as.character(result),
    status = status, detail = detail, stringsAsFactors = FALSE)
}
pass_fail <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"

# TCGA uses bracketed codes for missing values
clean_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "--", "'--", "[Not Available]", "[Not Applicable]",
             "[Not Evaluated]", "[Unknown]", "[Discrepancy]", "[Not Reported]",
             "not reported", "Not Reported", "unknown", "Unknown", "#N/A")] <- NA
  x
}
num <- function(x) suppressWarnings(as.numeric(clean_na(x)))

col_found <- list()
pick <- function(df, field, exact, pattern = NULL) {
  hit <- intersect(exact, names(df))
  if (!length(hit) && !is.null(pattern))
    hit <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  col_found[[field]] <<- if (length(hit)) hit[1] else "NOT FOUND"
  if (!length(hit)) return(rep(NA, nrow(df)))
  df[[hit[1]]]
}

get_file <- function(url, dest) {
  if (!file.exists(dest)) {
    ok <- tryCatch({ download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) { message("Download failed: ", basename(dest),
                                                 " (", conditionMessage(e), ")"); FALSE })
    if (!ok) { if (file.exists(dest)) file.remove(dest); return(NULL) }
  }
  dest
}

## 1. SAMPLES ====================================================
samples <- read.csv(file.path(meta_dir, "primary_tumor_1106_metadata.csv"), stringsAsFactors = FALSE)
ids <- as.character(samples$sample_id)
samples$patient  <- substr(ids, 1, 12)
samples$sample15 <- substr(ids, 1, 15)          # TCGA-XX-XXXX-01, used to match other platforms

# TCGA barcode: TCGA-<site>-<participant>-<sample type 01 = primary tumour><vial>
bc_ok <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01[A-Z]$", ids)
log_check("Sample", "Barcode format is a primary tumour (TCGA-XX-XXXX-01X)",
          sprintf("%d / %d", sum(bc_ok), length(ids)), pass_fail(all(bc_ok)),
          paste(head(ids[!bc_ok], 3), collapse = ", "))
log_check("Sample", "Sample IDs unique", sprintf("%d duplicates", sum(duplicated(ids))),
          pass_fail(!anyDuplicated(ids)))
log_check("Sample", "case_id equals patient part of sample barcode",
          sprintf("%d / %d", sum(samples$case_id == samples$patient), nrow(samples)),
          pass_fail(all(samples$case_id == samples$patient)))
fid_ok <- mapply(function(f, p) grepl(f, p, fixed = TRUE), samples$file_id, samples$file_path)
log_check("Sample", "file_id appears in its file path",
          sprintf("%d / %d", sum(fid_ok), nrow(samples)), pass_fail(all(fid_ok)))
log_check("Sample", "File IDs unique", sprintf("%d duplicates", sum(duplicated(samples$file_id))),
          pass_fail(!anyDuplicated(samples$file_id)))
multi <- unique(samples$patient[duplicated(samples$patient)])
log_check("Sample", "Patients with more than one primary tumour sample", length(multi), "INFO",
          "Only the first sample per patient enters the analysis set")
n_pat <- length(unique(samples$patient))

## 2. GDC CLINICAL ===============================================
gdc_file <- file.path(raw_dir, "gdc_clinical_raw.rds")
if (!file.exists(gdc_file)) saveRDS(GDCquery_clinic(project = "TCGA-BRCA", type = "clinical"), gdc_file)
gdc <- readRDS(gdc_file)
gdc <- gdc[!duplicated(gdc$submitter_id), ]
log_check("Clinical <-> Sample", "Patients with a GDC clinical record",
          sprintf("%d / %d", sum(unique(samples$patient) %in% gdc$submitter_id), n_pat),
          if (all(unique(samples$patient) %in% gdc$submitter_id)) "PASS" else "FLAG",
          paste(setdiff(unique(samples$patient), gdc$submitter_id), collapse = ", "))
gi <- match(samples$patient, gdc$submitter_id)
g  <- function(field, exact, pattern = NULL) pick(gdc, field, exact, pattern)[gi]

sex <- tolower(clean_na(g("sex", c("sex_at_birth", "gender", "demographic.gender"), "sex_at_birth|gender")))
race <- clean_na(g("race", c("race", "demographic.race")))
age <- num(g("age_at_index", c("age_at_index", "demographic.age_at_index")))
age_days <- num(g("age_at_diagnosis", c("age_at_diagnosis", "diagnoses.age_at_diagnosis")))
fill <- is.na(age) & !is.na(age_days)
age[fill] <- floor(age_days[fill] / 365.25)

vital <- clean_na(g("vital_status", c("vital_status", "demographic.vital_status")))
days_death <- num(g("days_to_death", c("days_to_death", "demographic.days_to_death")))
days_fu    <- num(g("days_to_last_follow_up", c("days_to_last_follow_up", "diagnoses.days_to_last_follow_up")))
dead <- tolower(vital) %in% "dead"
os_event_gdc <- ifelse(is.na(vital), NA, as.integer(dead))
os_days_gdc  <- ifelse(dead, days_death, days_fu)

stage_raw <- clean_na(g("ajcc_pathologic_stage", c("ajcc_pathologic_stage"), "pathologic_stage$"))
st <- toupper(stage_raw)
stage <- factor(ifelse(grepl("^STAGE (IV|III|II|I)", st), sub("^STAGE (IV|III|II|I).*$", "\\1", st), NA),
                levels = c("I", "II", "III", "IV"))
ajcc_edition <- clean_na(g("ajcc_staging_system_edition", c("ajcc_staging_system_edition"), "staging_system_edition"))

# TNM: keep the main category only (T1c -> T1, N1a -> N1); TX/NX/MX -> missing
tnm_main <- function(x, letter) {
  x <- toupper(clean_na(x)); x <- sub("^[CP]", "", x)
  ifelse(grepl(paste0("^", letter, "[0-4]"), x), substr(x, 1, 2), NA)
}
pT <- tnm_main(g("ajcc_pathologic_t", c("ajcc_pathologic_t"), "pathologic_t$"), "T")
pN <- tnm_main(g("ajcc_pathologic_n", c("ajcc_pathologic_n"), "pathologic_n$"), "N")
pM <- tnm_main(g("ajcc_pathologic_m", c("ajcc_pathologic_m"), "pathologic_m$"), "M")

histology <- clean_na(g("primary_diagnosis", c("primary_diagnosis")))
histology_group <- ifelse(is.na(histology), NA,
                          ifelse(grepl("lobular", histology, TRUE) & !grepl("duct", histology, TRUE), "Lobular",
                                 ifelse(grepl("duct", histology, TRUE) & !grepl("lobular", histology, TRUE), "Ductal", "Other/mixed")))

# Histologic grade: report whether any value exists (do not assume)
grade_cols <- grep("grade", names(gdc), value = TRUE, ignore.case = TRUE)
grade_n <- if (length(grade_cols)) max(sapply(grade_cols, function(cn) sum(!is.na(clean_na(gdc[[cn]])[gi])))) else 0
log_check("Clinical", "Histologic grade available", sprintf("%d samples with a value", grade_n),
          if (grade_n > 0) "INFO" else "FLAG",
          if (length(grade_cols)) paste("Columns checked:", paste(grade_cols, collapse = ", ")) else "No grade column")

## 3. PAM50 ======================================================
pam <- NULL
sub_file <- file.path(raw_dir, "tcga_brca_subtype_raw.rds")
if (!file.exists(sub_file)) {
  sub_raw <- tryCatch(TCGAquery_subtype(tumor = "brca"), error = function(e) NULL)
  if (!is.null(sub_raw)) saveRDS(sub_raw, sub_file)
}
if (file.exists(sub_file)) {
  sub_raw <- readRDS(sub_file)
  pcol <- grep("PAM50", names(sub_raw), value = TRUE)[1]
  if (!is.na(pcol) && "patient" %in% names(sub_raw)) {
    pam <- data.frame(patient = as.character(sub_raw$patient),
                      pam50_raw = as.character(sub_raw[[pcol]]), stringsAsFactors = FALSE)
    col_found$pam50 <- paste0("TCGAquery_subtype: ", pcol)
  }
}
if (is.null(pam)) stop("PAM50 subtypes not available (TCGAquery_subtype failed).")
pam <- pam[!is.na(pam$patient) & !duplicated(pam$patient), ]
pam_lab <- c(LumA = "Luminal A", LumB = "Luminal B", Her2 = "HER2-enriched",
             Basal = "Basal-like", Normal = "Normal-like")
pam50 <- factor(unname(pam_lab[pam$pam50_raw[match(samples$patient, pam$patient)]]),
                levels = unname(pam_lab))
log_check("Clinical <-> Sample", "Samples with a PAM50 subtype",
          sprintf("%d / %d", sum(!is.na(pam50)), length(ids)), "INFO")

## 4. BIOTAB: RECEPTOR STATUS AND MENOPAUSE ======================
biotab_rds <- file.path(raw_dir, "gdc_biotab_clinical_raw.rds")
if (!file.exists(biotab_rds)) {
  q <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical",
                data.type = "Clinical Supplement", data.format = "BCR Biotab")
  GDCdownload(q, directory = file.path(raw_dir, "GDCdata"))
  saveRDS(GDCprepare(q, directory = file.path(raw_dir, "GDCdata")), biotab_rds)
}
bt <- readRDS(biotab_rds)
pt <- as.data.frame(bt[["clinical_patient_brca"]], stringsAsFactors = FALSE)
pt <- pt[grepl("^TCGA-", pt$bcr_patient_barcode), ]           # drop the two description rows
pt <- pt[!duplicated(pt$bcr_patient_barcode), ]
log_check("Clinical <-> Sample", "Patients with a Biotab record (receptor status)",
          sprintf("%d / %d", sum(unique(samples$patient) %in% pt$bcr_patient_barcode), n_pat),
          if (all(unique(samples$patient) %in% pt$bcr_patient_barcode)) "PASS" else "FLAG")
bi <- match(samples$patient, pt$bcr_patient_barcode)

er_raw    <- clean_na(pick(pt, "ER status (IHC)", "er_status_by_ihc"))[bi]
pr_raw    <- clean_na(pick(pt, "PR status (IHC)", "pr_status_by_ihc"))[bi]
her2_ihc  <- clean_na(pick(pt, "HER2 status (IHC)", "her2_status_by_ihc"))[bi]
her2_fish <- clean_na(pick(pt, "HER2 status (FISH)", "her2_fish_status"))[bi]
menopause <- clean_na(pick(pt, "Menopause status", "menopause_status", "menopause"))[bi]
# No positive-lymph-node count exists in the Biotab file (only nodes examined); nodal status uses pN.

# ER/PR: recorded IHC status. Positive = >=1% staining under ASCO/CAP guidance
# (Allison et al. J Clin Oncol 2020;38:1346-1366).
binary_status <- function(x) ifelse(x %in% c("Positive", "Negative"), x, NA)
er_status <- binary_status(er_raw)
pr_status <- binary_status(pr_raw)

# HER2: ASCO/CAP logic (Wolff et al. J Clin Oncol 2018;36:2105-2122):
# IHC positive or ISH positive -> positive; IHC or ISH negative -> negative;
# IHC equivocal with no ISH result -> unknown.
her2_pos <- her2_ihc %in% "Positive" | her2_fish %in% "Positive"
her2_neg <- !her2_pos & (her2_ihc %in% "Negative" | her2_fish %in% "Negative")
her2_status <- ifelse(her2_pos, "Positive", ifelse(her2_neg, "Negative", NA))

hr_status <- ifelse(er_status %in% "Positive" | pr_status %in% "Positive", "Positive",
                    ifelse(er_status %in% "Negative" & pr_status %in% "Negative", "Negative", NA))
# Joint HR/HER2 groups (Howlader et al. JNCI 2014;106:dju055)
make_receptor_subtype <- function(hr, her2) {
  x <- ifelse(is.na(hr) | is.na(her2), NA,
              ifelse(hr == "Positive" & her2 == "Negative", "HR+/HER2-",
                     ifelse(hr == "Positive" & her2 == "Positive", "HR+/HER2+",
                            ifelse(hr == "Negative" & her2 == "Positive", "HR-/HER2+", "Triple-negative"))))
  factor(x, levels = c("HR+/HER2-", "HR+/HER2+", "HR-/HER2+", "Triple-negative"))
}
receptor_subtype <- make_receptor_subtype(hr_status, her2_status)

menopause_status <- ifelse(is.na(menopause), NA,
                           ifelse(grepl("^Pre|^Peri", menopause), "Pre/peri-menopausal",
                                  ifelse(grepl("^Post", menopause), "Post-menopausal", NA)))

## 5. TCGA-CDR SURVIVAL (Liu et al. 2018) ========================
cdr_file <- get_file("https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
                     file.path(raw_dir, "TCGA_CDR_SupplementalTableS1.xlsx"))
endpoint_cols <- c("OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")
cdr_m <- as.data.frame(matrix(NA_real_, length(ids), length(endpoint_cols),
                              dimnames = list(NULL, endpoint_cols)))
cdr_age <- cdr_stage <- rep(NA, length(ids))
if (!is.null(cdr_file)) {
  cdr <- suppressMessages(as.data.frame(readxl::read_excel(cdr_file, sheet = "TCGA-CDR", guess_max = 20000)))
  cdr <- cdr[cdr$type == "BRCA", ]
  cdr <- cdr[!duplicated(cdr$bcr_patient_barcode), ]
  ci <- match(samples$patient, cdr$bcr_patient_barcode)
  log_check("Clinical <-> Sample", "Patients in TCGA-CDR (curated survival)",
            sprintf("%d / %d", sum(unique(samples$patient) %in% cdr$bcr_patient_barcode), n_pat),
            if (all(unique(samples$patient) %in% cdr$bcr_patient_barcode)) "PASS" else "FLAG")
  for (e in intersect(endpoint_cols, names(cdr))) cdr_m[[e]] <- num(cdr[[e]])[ci]
  if ("age_at_initial_pathologic_diagnosis" %in% names(cdr))
    cdr_age <- num(cdr$age_at_initial_pathologic_diagnosis)[ci]
  if ("ajcc_pathologic_tumor_stage" %in% names(cdr)) {
    s <- toupper(clean_na(cdr$ajcc_pathologic_tumor_stage))[ci]
    cdr_stage <- ifelse(grepl("^STAGE (IV|III|II|I)", s), sub("^STAGE (IV|III|II|I).*$", "\\1", s), NA)
  }
} else {
  log_check("Clinical <-> Sample", "TCGA-CDR download", "failed", "FLAG",
            "Download TCGA-CDR-SupplementalTableS1.xlsx manually into metadata/raw_clinical and rerun")
}

# Cross-source consistency
age_diff <- abs(age - cdr_age)
log_check("Clinical sources", "Age: GDC vs TCGA-CDR (differ by >1 year)",
          sum(age_diff > 1, na.rm = TRUE), if (sum(age_diff > 1, na.rm = TRUE) == 0) "PASS" else "FLAG")
os_mis <- !is.na(os_event_gdc) & !is.na(cdr_m$OS) & os_event_gdc != cdr_m$OS
log_check("Clinical sources", "Vital status: GDC vs TCGA-CDR OS event",
          sum(os_mis), if (sum(os_mis) == 0) "PASS" else "FLAG",
          "TCGA-CDR is the curated source; use it for survival")
st_mis <- !is.na(stage) & !is.na(cdr_stage) & as.character(stage) != cdr_stage
log_check("Clinical sources", "Stage: GDC vs TCGA-CDR", sum(st_mis), if (sum(st_mis) == 0) "PASS" else "FLAG")
ed_tab <- table(ajcc_edition, useNA = "ifany")
log_check("Clinical", "AJCC staging editions used",
          paste(names(ed_tab), ed_tab, sep = ": ", collapse = "; "), "INFO",
          "TCGA cases span several AJCC editions; stage groups are used as recorded")

## 6. PANCANATLAS MOLECULAR DATA =================================
match15 <- function(df_ids, values) {
  key <- substr(df_ids, 1, 15)
  agg <- tapply(values, key, function(v) mean(v, na.rm = TRUE))
  n   <- tapply(values, key, length)
  list(value = unname(agg[samples$sample15]), n = unname(n[samples$sample15]))
}
find_id_col <- function(d) names(d)[sapply(d, function(z) mean(grepl("^TCGA-", z)) > 0.9)][1]

lf_val <- pur_val <- rep(NA_real_, length(ids))
lf_file <- get_file("https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8",
                    file.path(raw_dir, "PanCanAtlas_leukocyte_fraction_methylation.tsv"))
if (!is.null(lf_file)) {
  lf <- read.delim(lf_file, header = FALSE, stringsAsFactors = FALSE)
  idc <- find_id_col(lf)
  vc  <- names(lf)[sapply(lf, function(z) mean(!is.na(suppressWarnings(as.numeric(z)))) > 0.9)][1]
  m <- match15(lf[[idc]], suppressWarnings(as.numeric(lf[[vc]])))
  lf_val <- m$value
  log_check("Sample <-> PanCanAtlas", "Samples with methylation leukocyte fraction",
            sprintf("%d / %d", sum(!is.na(lf_val)), length(ids)), "INFO",
            sprintf("%d samples had >1 aliquot (averaged)", sum(m$n > 1, na.rm = TRUE)))
}

abs_file <- get_file("https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
                     file.path(raw_dir, "PanCanAtlas_ABSOLUTE_purity_ploidy.txt"))
if (!is.null(abs_file)) {
  ab <- read.delim(abs_file, stringsAsFactors = FALSE, check.names = FALSE)
  idc <- find_id_col(ab)
  pc  <- grep("^purity$", names(ab), ignore.case = TRUE, value = TRUE)[1]
  m <- match15(ab[[idc]], suppressWarnings(as.numeric(ab[[pc]])))
  pur_val <- m$value
  log_check("Sample <-> PanCanAtlas", "Samples with ABSOLUTE purity",
            sprintf("%d / %d", sum(!is.na(pur_val)), length(ids)), "INFO")
}

## 7. QUALITY-ANNOTATION EXCLUSIONS ==============================
# Rule (decided after reviewing every flag):
#  - patient-level TCGA notifications (prior systemic treatment, neoadjuvant therapy,
#    not invasive) apply to the patient on every platform -> exclude all their samples;
#  - Do_not_use on the mRNA-seq aliquot of the exact sample we use -> exclude that sample;
#  - Do_not_use on other platforms only (SNP array, methylation, miRNA, RPPA) -> not relevant.
excl_reason <- rep(NA_character_, length(ids))
qa_file <- get_file("https://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf",
                    file.path(raw_dir, "PanCanAtlas_merged_sample_quality_annotations.tsv"))
if (!is.null(qa_file)) {
  qa <- read.delim(qa_file, stringsAsFactors = FALSE, check.names = FALSE)
  col_of <- function(p) grep(p, names(qa), ignore.case = TRUE, value = TRUE)[1]
  c_type <- col_of("cancer"); c_plat <- col_of("platform"); c_id <- col_of("aliquot")
  c_dnu <- col_of("do_not_use"); c_pa <- col_of("patient_annotation")
  if (anyNA(c(c_type, c_plat, c_id, c_dnu, c_pa))) {
    log_check("Sample <-> PanCanAtlas", "Quality annotation columns", "not recognised", "FAIL",
              paste("Columns found:", paste(names(qa), collapse = ", ")))
  } else {
    qa <- qa[qa[[c_type]] == "BRCA", ]
    is_true <- function(v) toupper(as.character(v)) %in% c("TRUE", "T", "1", "YES")
    short_reason <- function(txt) ifelse(grepl("prior treatment", txt, ignore.case = TRUE),
                                         "Prior systemic treatment for another malignancy",
                                         ifelse(grepl("neoadjuvant", txt, ignore.case = TRUE), "Neoadjuvant systemic therapy",
                                                ifelse(grepl("DCIS", txt), "DCIS, not invasive",
                                                       ifelse(nzchar(trimws(txt)), txt, "Do_not_use on mRNA-seq aliquot"))))
    dnu <- qa[is_true(qa[[c_dnu]]), ]
    pat <- dnu[nzchar(trimws(dnu[[c_pa]])), ]
    pat_reason <- tapply(short_reason(pat[[c_pa]]), substr(pat[[c_id]], 1, 12), `[`, 1)
    mrna <- dnu[grepl("RNASeq", dnu[[c_plat]]) & !grepl("miRNA", dnu[[c_plat]]), ]
    mrna_reason <- tapply(short_reason(mrna[[c_pa]]), substr(mrna[[c_id]], 1, 16), `[`, 1)
    excl_reason <- ifelse(ids %in% names(mrna_reason), unname(mrna_reason[ids]),
                          ifelse(samples$patient %in% names(pat_reason), unname(pat_reason[samples$patient]), NA))
    log_check("Sample <-> PanCanAtlas", "Excluded on TCGA/PanCanAtlas quality annotations",
              sprintf("%d samples, %d patients", sum(!is.na(excl_reason)),
                      length(unique(samples$patient[!is.na(excl_reason)]))),
              "INFO", paste(names(table(excl_reason)), table(excl_reason), sep = ": ", collapse = "; "))
  }
} else {
  log_check("Sample <-> PanCanAtlas", "Quality annotation download", "failed", "FAIL",
            "Exclusions cannot be applied; download the file manually and rerun")
}
write.csv(data.frame(sample_id = ids[!is.na(excl_reason)], patient = samples$patient[!is.na(excl_reason)],
                     reason = excl_reason[!is.na(excl_reason)]),
          file.path(meta_dir, "pancanatlas_excluded_samples.csv"), row.names = FALSE)

## 8. ANALYSIS SET ===============================================
# Female, not excluded on quality annotations, first remaining primary tumour sample per patient
cand <- sex %in% "female" & is.na(excl_reason)
ord  <- order(samples$patient, ids)                 # 01A before 01B
keep <- ord[cand[ord]]
first <- logical(length(ids)); first[keep[!duplicated(samples$patient[keep])]] <- TRUE
in_analysis_set <- cand & first

exclusion_reason <- ifelse(is.na(sex), "Sex unknown / no clinical record",
                           ifelse(sex == "male", "Male",
                                  ifelse(!is.na(excl_reason), paste("Quality annotation:", excl_reason),
                                         ifelse(!first, "Additional sample from same patient", NA))))
stopifnot(identical(in_analysis_set, is.na(exclusion_reason)))

set_summary <- data.frame(
  step = c("Total samples", names(table(exclusion_reason)), "Analysis set"),
  n = c(length(ids), as.integer(table(exclusion_reason)), sum(in_analysis_set)))
write.csv(set_summary, file.path(out_dir, "analysis_set_summary.csv"), row.names = FALSE)

## 9. SAVE =======================================================
clinical <- data.frame(
  sample_id = ids, patient = samples$patient, file_id = samples$file_id,
  in_analysis_set = in_analysis_set, exclusion_reason = exclusion_reason,
  multi_sample_patient = samples$patient %in% multi,
  sex = sex, age = age, race = race, menopause_status = menopause_status,
  pam50 = pam50,
  er_status = er_status, pr_status = pr_status, her2_status = her2_status,
  her2_ihc = her2_ihc, her2_fish = her2_fish, hr_status = hr_status,
  receptor_subtype = receptor_subtype,
  stage = stage, stage_raw = stage_raw, ajcc_edition = ajcc_edition,
  pT = pT, pN = pN, pM = pM,
  histology = histology, histology_group = histology_group,
  os_event_gdc = os_event_gdc, os_days_gdc = os_days_gdc,
  cdr_m,
  leukocyte_fraction_methylation = round(lf_val, 4), absolute_purity = round(pur_val, 4),
  check.names = FALSE, stringsAsFactors = FALSE)

write.csv(clinical, file.path(meta_dir, "clinical_variables.csv"), row.names = FALSE, na = "")
saveRDS(clinical, file.path(meta_dir, "clinical_variables.rds"))

checks_df <- do.call(rbind, checks)
write.csv(checks_df, file.path(out_dir, "clinical_data_checks.csv"), row.names = FALSE)
write.csv(data.frame(field = names(col_found), column = unlist(col_found)),
          file.path(out_dir, "clinical_source_columns_used.csv"), row.names = FALSE)

## 10. REPORT ====================================================
op <- options(width = 200)
cat("\n================ CLINICAL DATA AND ANALYSIS SET ================\n")
print(checks_df[, c("link", "check", "result", "status")], row.names = FALSE)
cat("\nSource columns used:\n"); print(unlist(col_found))
cat("\nAnalysis set:\n"); print(set_summary, row.names = FALSE)
cat("\nSaved: metadata/clinical_variables.csv and metadata/pancanatlas_excluded_samples.csv\n")
cat("Next: run phase4_immune_landscape.R\n")
options(op)