#Clinical Immune Integration

for (p in c("ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(library(ggplot2))
set.seed(2026)

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
meta_dir <- file.path(project_dir, "metadata")
out_dir  <- file.path(project_dir, "results", "clinical_immune_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

m    <- readRDS(file.path(meta_dir, "clinical_master_table.rds"))
land <- readRDS(file.path(project_dir, "results", "immune_landscape", "sample_level_results.rds"))
a <- m[m$in_analysis_set, ]
if (!setequal(a$sample_id, land$sample_id))
  stop("Master table and immune landscape results differ. Run scripts 1-3 first.")
li <- match(a$sample_id, land$sample_id)
for (f in c("CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial")) a[[f]] <- land[[f]][li]

# Stromal score built the same way as the immune score (mean of standardised features)
a$stroma_score <- rowMeans(scale(cbind(a$CAFs, a$Endothelial)))

lab_order <- unique(a$cluster_label[order(a$cluster)])
a$cluster_name <- factor(sub("^C[0-9]+: ", "", a$cluster_label), levels = sub("^C[0-9]+: ", "", lab_order))
cl_cols <- setNames(hcl.colors(nlevels(a$cluster_name), "Dark 3"), levels(a$cluster_name))

# Ordered categories
a$pT <- factor(a$pT, levels = c("T1", "T2", "T3", "T4"))
a$pN <- factor(a$pN, levels = c("N0", "N1", "N2", "N3"))
a$pM <- factor(a$pM, levels = c("M0", "M1"))
a$histology_group <- factor(a$histology_group, levels = c("Ductal", "Lobular", "Other/mixed"))
a$menopause_status <- factor(a$menopause_status, levels = c("Pre/peri-menopausal", "Post-menopausal"))
for (v in c("er_status", "pr_status", "her2_status", "her2_status_strict"))
  a[[v]] <- factor(a[[v]], levels = c("Negative", "Positive"))
a$race <- factor(a$race)

cat_vars <- c("PAM50" = "pam50",
              "Receptor subtype" = "receptor_subtype",
              "Receptor subtype (strict HER2)" = "receptor_subtype_strict",
              "ER status" = "er_status", "PR status" = "pr_status",
              "HER2 status" = "her2_status", "HER2 status (strict)" = "her2_status_strict",
              "Stage" = "stage", "T stage" = "pT", "N stage" = "pN", "M stage" = "pM",
              "Histology" = "histology_group", "Menopausal status" = "menopause_status",
              "Race" = "race")
step_of <- c("PAM50" = "33 Subtype", "Receptor subtype" = "33 Subtype",
             "Receptor subtype (strict HER2)" = "33 Subtype",
             "ER status" = "35 Receptors", "PR status" = "35 Receptors",
             "HER2 status" = "35 Receptors", "HER2 status (strict)" = "35 Receptors",
             "Stage" = "34 Stage", "T stage" = "34 Stage", "N stage" = "34 Stage",
             "M stage" = "34 Stage", "Histology" = "36 Other", "Menopausal status" = "36 Other",
             "Race" = "36 Other", "Age" = "36 Other")

cramers_v <- function(tb) {
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tb) * (min(dim(tb)) - 1))))
}
eps2 <- function(H, n) as.numeric(H / (n - 1))
med_iqr <- function(x) sprintf("%.2f (%.2f to %.2f)", median(x), quantile(x, .25), quantile(x, .75))

#Independent Validation
val_rows <- list()
refs <-c(leukocyte_fraction_methylation = "Leukocyte fraction (DNA methylation)",
         absolute_purity = "Tumour purity (ABSOLUTE)")
for (r in names(refs)) {
  for (f in c("immune_score", "CD8_T_cells", "B_cells", "Macrophages", "stroma_score", "CAFs", "Endothelial")) {
    ok <- !is.na(a[[r]])
    ct <- suppressWarnings(cor.test(a[[f]][ok], a[[r]][ok], method = "spearman", exact = FALSE))
    val_rows[[length(val_rows) + 1]] <- data.frame(
      reference = refs[[r]], rna_feature = f, n = sum(ok),
      spearman_rho = round(unname(ct$estimate), 3), p_value = signif(ct$p.value, 3))
  }
}
validation <- do.call(rbind, val_rows)
write.csv(validation, file.path(out_dir, "rna_scores_vs_methylation_and_purity.csv"), row.names = FALSE)

val_cluster <- do.call(rbind, lapply(names(refs), function(r) {
  ok <- !is.na(a[[r]])
  kw <- kruskal.test(a[[r]][ok] ~ a$cluster_name[ok])
  meds <- tapply(a[[r]][ok], a$cluster_name[ok], med_iqr)
  data.frame(reference = refs[[r]], t(meds), kruskal_p = signif(kw$p.value, 3),
             epsilon_squared = round(eps2(kw$statistic, sum(ok)), 3), check.names = FALSE)
}))
write.csv(val_cluster, file.path(out_dir, "methylation_and_purity_by_cluster.csv"), row.names = FALSE)

assoc <- list(); by_level <- list(); tab1 <- list()
for (lab in names(cat_vars)) {
  v <- a[[cat_vars[[lab]]]]
  v <- droplevels(if (is.factor(v)) v else factor(v))
  ok <- !is.na(v)
  if (nlevels(v[ok]) < 2) next
  vv <- v[ok]; cl <- a$cluster_name[ok]
  
  tb <- table(vv, cl)
  p_cl <- suppressWarnings(chisq.test(tb, simulate.p.value = TRUE, B = 10000)$p.value)
  kw_i <- kruskal.test(a$immune_score[ok] ~ vv)
  kw_s <- kruskal.test(a$stroma_score[ok] ~ vv)
  assoc[[lab]] <- data.frame(
    step = step_of[[lab]], variable = lab, n = sum(ok), n_categories = nlevels(vv),
    outcome = c("Immune cluster", "Immune score", "Stromal score"),
    test = c("Chi-square (Monte-Carlo)", "Kruskal-Wallis", "Kruskal-Wallis"),
    effect_measure = c("Cramer's V", "Epsilon-squared", "Epsilon-squared"),
    effect_size = round(c(cramers_v(tb), eps2(kw_i$statistic, sum(ok)), eps2(kw_s$statistic, sum(ok))), 3),
    p_value = c(p_cl, kw_i$p.value, kw_s$p.value))
  
  by_level[[lab]] <- do.call(rbind, lapply(levels(vv), function(lv) {
    s <- ok & v %in% lv
    data.frame(variable = lab, category = lv, n = sum(s),
               immune_score_median_iqr = med_iqr(a$immune_score[s]),
               stroma_score_median_iqr = med_iqr(a$stroma_score[s]),
               t(setNames(sprintf("%.1f%%", 100 * as.numeric(tb[lv, ]) / sum(tb[lv, ])),
                          paste0("in ", colnames(tb)))),
               check.names = FALSE)
  }))
  
  # Cohort characteristics by cluster: n (% of column)
  tab1[[lab]] <- do.call(rbind, c(
    list(data.frame(variable = lab, category = "", Overall = "",
                    t(setNames(rep("", nlevels(a$cluster_name)), levels(a$cluster_name))),
                    check.names = FALSE)),
    lapply(c(levels(vv), "Missing"), function(lv) {
      s <- if (lv == "Missing") is.na(v) else (v %in% lv)
      cell <- function(idx) sprintf("%d (%.1f%%)", sum(s & idx), 100 * sum(s & idx) / sum(idx))
      data.frame(variable = "", category = lv, Overall = cell(rep(TRUE, nrow(a))),
                 t(sapply(levels(a$cluster_name), function(cn) cell(a$cluster_name == cn))),
                 check.names = FALSE)
    })))
}

# Age (continuous)
kw_age <- kruskal.test(a$age ~ a$cluster_name)
sp_i <- suppressWarnings(cor.test(a$age, a$immune_score, method = "spearman", exact = FALSE))
sp_s <- suppressWarnings(cor.test(a$age, a$stroma_score, method = "spearman", exact = FALSE))
assoc[["Age"]] <- data.frame(
  step = step_of[["Age"]], variable = "Age", n = sum(!is.na(a$age)), n_categories = NA,
  outcome = c("Immune cluster", "Immune score", "Stromal score"),
  test = c("Kruskal-Wallis (age across clusters)", "Spearman", "Spearman"),
  effect_measure = c("Epsilon-squared", "Spearman rho", "Spearman rho"),
  effect_size = round(c(eps2(kw_age$statistic, sum(!is.na(a$age))),
                        unname(sp_i$estimate), unname(sp_s$estimate)), 3),
  p_value = c(kw_age$p.value, sp_i$p.value, sp_s$p.value))
age_row <- data.frame(variable = "Age, years", category = "median (IQR)",
                      Overall = med_iqr(a$age),
                      t(sapply(levels(a$cluster_name), function(cn) med_iqr(a$age[a$cluster_name == cn]))),
                      check.names = FALSE)

assoc_df <- do.call(rbind, assoc)
assoc_df$q_value_fdr <- p.adjust(assoc_df$p_value, method = "BH")
assoc_df$p_value <- signif(assoc_df$p_value, 3)
assoc_df$q_value_fdr <- signif(assoc_df$q_value_fdr, 3)
rownames(assoc_df) <- NULL
write.csv(assoc_df, file.path(out_dir, "clinical_variable_associations.csv"), row.names = FALSE)

by_level_df <- do.call(rbind, by_level); rownames(by_level_df) <- NULL
write.csv(by_level_df, file.path(out_dir, "immune_and_stromal_scores_by_clinical_category.csv"), row.names = FALSE)

n_row <- data.frame(variable = "Tumours", category = "n", Overall = as.character(nrow(a)),
                    t(sapply(levels(a$cluster_name), function(cn) as.character(sum(a$cluster_name == cn)))),
                    check.names = FALSE)
tab1_df <- rbind(n_row, age_row, do.call(rbind, tab1)); rownames(tab1_df) <- NULL
write.csv(tab1_df, file.path(out_dir, "cohort_characteristics_by_cluster.csv"), row.names = FALSE)

# Adjusted Models
# Do immune/stromal scores differ by stage, age and histology after accounting for PAM50?
mdat <- a[, c("immune_score", "stroma_score", "pam50", "stage", "age", "histology_group")]
mdat <- mdat[complete.cases(mdat), ]
mdat$pam50 <- relevel(droplevels(mdat$pam50), ref = "Luminal A")
mdat$stage <- droplevels(mdat$stage)
mdat$histology_group <- droplevels(mdat$histology_group)
mdat$age_per_10y <- mdat$age / 10

adj_terms <- list(); adj_coef <- list()
for (y in c("immune_score", "stroma_score")) {
  fit <- lm(as.formula(paste(y, "~ pam50 + stage + age_per_10y + histology_group")), data = mdat)
  d1 <- drop1(fit, test = "F")
  adj_terms[[y]] <- data.frame(outcome = y, term = rownames(d1)[-1],
                               df = d1$Df[-1], F_value = round(d1$`F value`[-1], 2),
                               p_value = signif(d1$`Pr(>F)`[-1], 3))
  ci <- confint(fit)
  adj_coef[[y]] <- data.frame(outcome = y, coefficient = names(coef(fit)),
                              estimate = round(coef(fit), 3),
                              ci_low = round(ci[, 1], 3), ci_high = round(ci[, 2], 3),
                              p_value = signif(summary(fit)$coefficients[, 4], 3))
}
adj_terms_df <- do.call(rbind, adj_terms); rownames(adj_terms_df) <- NULL
adj_coef_df  <- do.call(rbind, adj_coef);  rownames(adj_coef_df) <- NULL
write.csv(adj_terms_df, file.path(out_dir, "adjusted_model_terms.csv"), row.names = FALSE)
write.csv(adj_coef_df, file.path(out_dir, "adjusted_model_coefficients.csv"), row.names = FALSE)

#Figures
ok <- !is.na(a$leukocyte_fraction_methylation)
rho_lf <- validation$spearman_rho[validation$reference == refs[[1]] & validation$rna_feature == "immune_score"]
p_val1 <- ggplot(a[ok, ], aes(immune_score, leukocyte_fraction_methylation, colour = cluster_name)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_colour_manual(values = cl_cols) +
  labs(x = "Immune score (RNA, consensus)", y = "Leukocyte fraction (DNA methylation)", colour = NULL,
       title = sprintf("RNA immune score vs methylation leukocyte fraction (Spearman rho = %.2f, n = %d)",
                       rho_lf, sum(ok))) +
  theme_bw()
ggsave(file.path(out_dir, "immune_score_vs_methylation_leukocyte_fraction.png"), p_val1,
       width = 8, height = 5.5, dpi = 200)

vlong <- rbind(
  data.frame(measure = refs[[1]], cluster = a$cluster_name, value = a$leukocyte_fraction_methylation),
  data.frame(measure = refs[[2]], cluster = a$cluster_name, value = a$absolute_purity))
vlong <- vlong[!is.na(vlong$value), ]
p_val2 <- ggplot(vlong, aes(cluster, value, fill = cluster)) +
  geom_boxplot(outlier.size = 0.4) + facet_wrap(~ measure, scales = "free_y") +
  scale_fill_manual(values = cl_cols) +
  labs(x = NULL, y = NULL, title = "DNA-based measures by immune cluster") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_dir, "methylation_and_purity_by_cluster_boxplot.png"), p_val2,
       width = 9, height = 4.5, dpi = 200)

# Scores by subtype and stage
sv <- c("PAM50" = "pam50", "Receptor subtype" = "receptor_subtype", "Stage" = "stage")
slong <- do.call(rbind, lapply(names(sv), function(lab) {
  v <- a[[sv[[lab]]]]
  rbind(data.frame(variable = lab, category = as.character(v), score = "Immune score", value = a$immune_score),
        data.frame(variable = lab, category = as.character(v), score = "Stromal score", value = a$stroma_score))
}))
slong <- slong[!is.na(slong$category), ]
slong$category <- factor(slong$category, levels = unique(c(levels(a$pam50), levels(a$receptor_subtype),
                                                           levels(a$stage))))
p_sc <- ggplot(slong, aes(category, value)) +
  geom_boxplot(outlier.size = 0.3, fill = "grey85") +
  facet_grid(score ~ variable, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Score (mean standardised features)",
       title = "Immune and stromal scores by subtype and stage") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "immune_and_stromal_scores_by_subtype_and_stage.png"), p_sc,
       width = 12, height = 6.5, dpi = 200)

# Cluster mix within each clinical group
cv <- c("PAM50" = "pam50", "Receptor subtype" = "receptor_subtype", "Stage" = "stage",
        "Histology" = "histology_group")
clong <- do.call(rbind, lapply(names(cv), function(lab)
  data.frame(variable = lab, category = as.character(a[[cv[[lab]]]]), cluster = a$cluster_name)))
clong <- clong[!is.na(clong$category), ]
clong$category <- factor(clong$category, levels = unique(c(levels(a$pam50), levels(a$receptor_subtype),
                                                           levels(a$stage), levels(a$histology_group))))
clong$variable <- factor(clong$variable, levels = names(cv))
p_mix <- ggplot(clong, aes(category, fill = cluster)) +
  geom_bar(position = "fill") +
  facet_grid(. ~ variable, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cl_cols) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  labs(x = NULL, y = "Tumours", fill = "Immune cluster",
       title = "Immune cluster composition within clinical groups") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "cluster_composition_by_clinical_group.png"), p_mix,
       width = 12, height = 5, dpi = 200)

# Effect-size overview
eff <- assoc_df[!grepl("strict", assoc_df$variable), ]
eff$abs_effect <- abs(eff$effect_size)
eff$significant <- ifelse(eff$q_value_fdr < 0.05, "FDR < 0.05", "not significant")
eff$variable <- factor(eff$variable, levels = rev(unique(eff$variable)))
p_eff <- ggplot(eff, aes(outcome, variable, size = abs_effect, colour = significant)) +
  geom_point() +
  scale_colour_manual(values = c("FDR < 0.05" = "#B2182B", "not significant" = "grey65")) +
  scale_size_area(max_size = 9) +
  labs(x = NULL, y = NULL, size = "Effect size (|value|)", colour = NULL,
       title = "Associations of clinical variables with immune outcomes",
       caption = "Cramer's V for cluster; epsilon-squared for scores; |Spearman rho| for age vs scores") +
  theme_bw()
ggsave(file.path(out_dir, "clinical_association_effect_sizes.png"), p_eff, width = 8, height = 6, dpi = 200)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_session_info.txt"))

op <- options(width = 200)
cat("Tumours:", nrow(a), "\n")
cat("\nA. RNA scores vs DNA-based measures (Spearman):\n")
print(validation, row.names = FALSE)
cat("\nDNA-based measures by cluster, median (IQR):\n"); print(val_cluster, row.names = FALSE)
cat("\nB-E. Associations (all outcomes, FDR-adjusted):\n")
print(assoc_df[, c("step", "variable", "n", "outcome", "effect_measure", "effect_size", "p_value", "q_value_fdr")],
      row.names = FALSE)
cat("\nScores and cluster mix by category:\n")
print(by_level_df, row.names = FALSE)
cat("\nF. Adjusted models (each term adjusted for the others; n =", nrow(mdat), "):\n")
print(adj_terms_df, row.names = FALSE)
cat("\nSaved to:", out_dir, "\n"); print(sort(list.files(out_dir)))
options(op)
