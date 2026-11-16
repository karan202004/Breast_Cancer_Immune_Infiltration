# Clinical master table and matching QC

for (p in c("ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(library(ggplot2))

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
meta_dir <- file.path(project_dir, "metadata")
out_dir  <- file.path(project_dir, "results", "clinical_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

checks <- list()
log_check <- function(link, check, result, status, detail = "") {
  checks[[length(checks) + 1]] <<- data.frame(
    link = link, check = check, result = as.character(result),
    status = status, detail = detail, stringsAsFactors = FALSE)
}
pass_fail <- function(ok) if (isTRUE(ok)) "PASS" else "FAIL"

clinical <- readRDS(file.path(meta_dir, "clinical_variables.rds"))
ids <- as.character(clinical$sample_id)
samples <- read.csv(file.path(meta_dir, "primary_tumor_1106_metadata.csv"), stringsAsFactors = FALSE)
log_check("Clinical <-> Sample", "Clinical table has the 1,106 samples in metadata order",
          identical(ids, as.character(samples$sample_id)),
          pass_fail(identical(ids, as.character(samples$sample_id))))

#Expression
tpm <- read.delim(file.path(project_dir, "data", "CIBERSORT", "BRCA_TPM_1106.txt"),
                  check.names = FALSE, stringsAsFactors = FALSE)
genes <- tpm[[1]]; tpm <- as.matrix(tpm[, -1]); storage.mode(tpm) <- "numeric"
rownames(tpm) <- genes
colnames(tpm) <- gsub(".", "-", colnames(tpm), fixed = TRUE)

expr_same <- setequal(colnames(tpm), ids) && ncol(tpm) == length(ids)
log_check("Sample <-> Expression", "Expression matrix has exactly the 1,106 sample IDs",
          sprintf("%d of %d IDs found", sum(ids %in% colnames(tpm)), length(ids)), pass_fail(expr_same))
log_check("Sample <-> Expression", "No duplicated expression columns",
          sum(duplicated(colnames(tpm))), pass_fail(!anyDuplicated(colnames(tpm))))
if (!expr_same) stop("Expression sample IDs do not match the metadata. See the check log.")
tpm <- tpm[, ids]

alt_file <- file.path(project_dir, "results", "normalized_matrix", "gdc_tpm_1106.rds")
if (file.exists(alt_file)) {
  alt_ids <- tryCatch(colnames(readRDS(alt_file)), error = function(e) NULL)
  if (!is.null(alt_ids)) {
    if (sum(alt_ids %in% samples$file_id) > sum(alt_ids %in% ids))
      alt_ids <- samples$sample_id[match(alt_ids, samples$file_id)]
    alt_ids <- gsub(".", "-", alt_ids, fixed = TRUE)
    log_check("Sample <-> Expression", "gdc_tpm_1106.rds has the same 1,106 samples",
              sprintf("%d of %d IDs found", sum(ids %in% alt_ids), length(ids)),
              pass_fail(setequal(alt_ids, ids)))
  }
}

#Immune results
cons <- readRDS(file.path(project_dir, "results", "method_comparison_consensus",
                          "06_consensus_immune_features.rds"))
land <- readRDS(file.path(project_dir, "results", "immune_landscape", "sample_level_results.rds"))
log_check("Expression <-> Immune", "Consensus immune features: same samples, same order",
          identical(as.character(cons$sample_id), ids), pass_fail(identical(as.character(cons$sample_id), ids)))
log_check("Expression <-> Immune", "Immune landscape results: one row per sample",
          sum(duplicated(land$sample_id)), pass_fail(!anyDuplicated(land$sample_id)))
land_ok <- setequal(land$sample_id, ids[clinical$in_analysis_set])
log_check("Expression <-> Immune", "Immune landscape results cover exactly the analysis set",
          sprintf("analysis set %d; landscape results %d", sum(clinical$in_analysis_set), nrow(land)),
          pass_fail(land_ok), if (land_ok) "" else "Run phase4_immune_landscape.R before this script")
if (!land_ok) stop("Immune landscape results do not match the analysis set. Run phase4_immune_landscape.R first.")

#Clincal records VS Gene expression
lg <- function(g) if (g %in% rownames(tpm)) log2(tpm[g, ] + 1) else rep(NA_real_, ncol(tpm))

# Sex from expression (sample-swap check; Toker et al. F1000Research 2016;5:2103).
# Y-linked genes decide: XIST loss is common in breast tumours, so low XIST alone is not "male".
y_genes <- intersect(c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"), rownames(tpm))
y_expr  <- colMeans(log2(tpm[y_genes, , drop = FALSE] + 1))
xist    <- lg("XIST")
sex_expr <- ifelse(y_expr > 2, "male", ifelse(y_expr < 1, "female", "ambiguous"))
sex_mismatch <- !is.na(clinical$sex) & sex_expr != "ambiguous" & sex_expr != clinical$sex
log_check("Clinical <-> Expression", "Clinical sex matches Y-chromosome gene expression",
          sprintf("%d mismatches, %d ambiguous", sum(sex_mismatch), sum(sex_expr == "ambiguous")),
          if (sum(sex_mismatch) == 0) "PASS" else "FLAG",
          paste(head(ids[sex_mismatch], 5), collapse = ", "))

# Receptor status vs expression of its gene
auc <- function(score, label) {
  ok <- !is.na(score) & label %in% c("Positive", "Negative")
  s <- score[ok]; l <- label[ok] == "Positive"; r <- rank(s)
  (sum(r[l]) - sum(l) * (sum(l) + 1) / 2) / (sum(l) * sum(!l))
}
discordant <- function(score, label) {
  neg_med <- median(score[label %in% "Negative"], na.rm = TRUE)
  pos_med <- median(score[label %in% "Positive"], na.rm = TRUE)
  (label %in% "Positive" & score < neg_med) | (label %in% "Negative" & score > pos_med)
}
conc <- list(ER = list(gene = "ESR1", status = clinical$er_status),
             PR = list(gene = "PGR", status = clinical$pr_status),
             HER2 = list(gene = "ERBB2", status = clinical$her2_status))
conc_flags <- list(); conc_rows <- list(); conc_long <- list()
for (r in names(conc)) {
  g <- lg(conc[[r]]$gene); st <- conc[[r]]$status
  d <- discordant(g, st); a <- auc(g, st)
  conc_flags[[r]] <- d
  conc_rows[[r]] <- data.frame(receptor = r, gene = conc[[r]]$gene,
                               n_positive = sum(st %in% "Positive"), n_negative = sum(st %in% "Negative"),
                               AUC = round(a, 3), n_discordant = sum(d))
  log_check("Clinical <-> Expression", sprintf("%s status vs %s expression (AUC)", r, conc[[r]]$gene),
            sprintf("AUC %.3f; %d discordant", a, sum(d)), if (isTRUE(a >= 0.85)) "PASS" else "FLAG",
            "Discordant = positive but below negative-group median, or the reverse")
  conc_long[[r]] <- data.frame(receptor = paste0(r, " (", conc[[r]]$gene, ")"),
                               status = ifelse(is.na(st), "Unknown", st), expr = g)
}
conc_tab <- do.call(rbind, conc_rows)
write.csv(conc_tab, file.path(out_dir, "receptor_expression_concordance.csv"), row.names = FALSE)

# Stricter HER2 (sensitivity definition, data-driven, not a guideline):
# a recorded positive is kept only if FISH-positive or ERBB2 mRNA is above the 95th
erbb2 <- lg("ERBB2")
neg95 <- quantile(erbb2[clinical$her2_status %in% "Negative"], 0.95, na.rm = TRUE)
her2_supported <- clinical$her2_status %in% "Positive" &
  (clinical$her2_fish %in% "Positive" | erbb2 > neg95)
her2_status_strict <- ifelse(clinical$her2_status %in% "Negative", "Negative",
                             ifelse(her2_supported, "Positive", NA))
rs <- ifelse(is.na(clinical$hr_status) | is.na(her2_status_strict), NA,
             ifelse(clinical$hr_status == "Positive" & her2_status_strict == "Negative", "HR+/HER2-",
                    ifelse(clinical$hr_status == "Positive" & her2_status_strict == "Positive", "HR+/HER2+",
                           ifelse(clinical$hr_status == "Negative" & her2_status_strict == "Positive", "HR-/HER2+", "Triple-negative"))))
receptor_subtype_strict <- factor(rs, levels = levels(clinical$receptor_subtype))
log_check("Clinical <-> Expression", "HER2-positive calls supported by FISH or high ERBB2 mRNA",
          sprintf("%d of %d", sum(her2_supported), sum(clinical$her2_status %in% "Positive")), "INFO",
          "Unsupported calls are missing in her2_status_strict (sensitivity analysis)")

# PAM50 vs clinical receptor subtype (analysis set)
a_rows <- clinical$in_analysis_set
pam_rec <- table(PAM50 = clinical$pam50[a_rows], receptor_subtype = clinical$receptor_subtype[a_rows],
                 useNA = "ifany")
write.csv(as.data.frame.matrix(pam_rec), file.path(out_dir, "PAM50_vs_receptor_subtype.csv"))

#Master Table
li <- match(ids, land$sample_id)
master <- data.frame(
  clinical,
  sex_expression = sex_expr, sex_mismatch = sex_mismatch,
  her2_status_strict = her2_status_strict, receptor_subtype_strict = receptor_subtype_strict,
  ERBB2_log2tpm = round(erbb2, 3),
  er_expression_discordant = conc_flags$ER, pr_expression_discordant = conc_flags$PR,
  her2_expression_discordant = conc_flags$HER2,
  cluster = land$cluster[li], cluster_label = land$cluster_label[li],
  immune_score = land$immune_score[li], immune_group = land$immune_group[li],
  check.names = FALSE, stringsAsFactors = FALSE)
log_check("Clinical <-> Immune", "Every analysis-set sample has an immune cluster",
          sprintf("%d / %d", sum(!is.na(master$cluster[a_rows])), sum(a_rows)),
          pass_fail(all(!is.na(master$cluster[a_rows]))))

write.csv(master, file.path(meta_dir, "clinical_master_table.csv"), row.names = FALSE, na = "")
saveRDS(master, file.path(meta_dir, "clinical_master_table.rds"))

flag_cols <- c("sex_mismatch", "er_expression_discordant", "pr_expression_discordant",
               "her2_expression_discordant")
any_flag <- Reduce(`|`, lapply(flag_cols, function(cn) master[[cn]] %in% TRUE))
flagged <- master[any_flag, c("sample_id", "in_analysis_set", flag_cols, "pam50",
                              "er_status", "pr_status", "her2_status", "absolute_purity")]
write.csv(flagged, file.path(out_dir, "expression_discordant_samples.csv"), row.names = FALSE, na = "")

# Completeness of each variable in the analysis set
vars <- c("age", "race", "menopause_status", "pam50", "er_status", "pr_status", "her2_status",
          "receptor_subtype", "her2_status_strict", "receptor_subtype_strict",
          "stage", "pT", "pN", "pM", "histology_group",
          "OS", "PFI", "DFI", "DSS", "leukocyte_fraction_methylation", "absolute_purity")
a <- master[a_rows, ]
completeness <- data.frame(
  variable = vars,
  n_available = sapply(vars, function(v) sum(!is.na(a[[v]]))),
  pct_available = sapply(vars, function(v) round(100 * mean(!is.na(a[[v]])), 1)),
  levels_or_range = sapply(vars, function(v) {
    x <- a[[v]]
    if (all(is.na(x))) return("no values")
    if (is.numeric(x)) sprintf("%.2f to %.2f", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
    else paste(names(table(x)), table(x), sep = " = ", collapse = "; ")
  }), row.names = NULL)
write.csv(completeness, file.path(out_dir, "variable_completeness.csv"), row.names = FALSE)

checks_df <- do.call(rbind, checks)
write.csv(checks_df, file.path(out_dir, "sample_matching_qc_checks.csv"), row.names = FALSE)

#QC Figures
p_sex <- ggplot(data.frame(xist = xist, y = y_expr,
                           sex = ifelse(is.na(clinical$sex), "unknown", clinical$sex)),
                aes(xist, y, colour = sex)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_hline(yintercept = c(1, 2), linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(female = "#C2185B", male = "#1565C0", unknown = "grey60")) +
  labs(x = "XIST, log2(TPM + 1)", y = "Y-chromosome genes, mean log2(TPM + 1)",
       colour = "Clinical sex", title = "Sex check: clinical record vs expression") + theme_bw()
ggsave(file.path(out_dir, "sex_check_expression.png"), p_sex, width = 6.5, height = 5, dpi = 200)

cl_df <- do.call(rbind, conc_long)
cl_df$status <- factor(cl_df$status, levels = c("Positive", "Negative", "Unknown"))
p_rec <- ggplot(cl_df, aes(status, expr, fill = status)) +
  geom_boxplot(outlier.size = 0.4) + facet_wrap(~ receptor, nrow = 1) +
  scale_fill_manual(values = c(Positive = "#D95F02", Negative = "#1B9E77", Unknown = "grey75")) +
  labs(x = "Clinical status", y = "Gene expression, log2(TPM + 1)",
       title = "Receptor status vs expression of its gene") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_dir, "receptor_status_vs_expression.png"), p_rec, width = 10, height = 4.5, dpi = 200)

op <- options(width = 200)

print(checks_df[, c("link", "check", "result", "status")], row.names = FALSE)
cat("\nReceptor status vs expression:\n"); print(conc_tab, row.names = FALSE)
cat("\nPAM50 vs receptor subtype, analysis set (counts):\n"); print(pam_rec)
cat("\nVariable completeness in the analysis set (n =", nrow(a), "):\n")
print(completeness[, 1:3], row.names = FALSE)
cat("\nSamples with an expression discordance flag:", nrow(flagged),
    "(", sum(flagged$in_analysis_set), "in the analysis set )\n")
cat("\nSaved: metadata/clinical_master_table.csv and results/clinical_integration/\n")
options(op)
