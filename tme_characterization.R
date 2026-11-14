# ============================================================
# TUMOUR MICROENVIRONMENT CHARACTERISATION  (steps 37–43)
# TCGA-BRCA | 1,069 tumours in the analysis set
#
# Run order:  1. clinical_data_and_analysis_set.R
#             2. phase4_immune_landscape.R
#             3. clinical_master_table_qc.R
#             4. clinical_immune_integration.R
#             5. tme_characterization.R   (this script)
#
# Uses the existing three immune clusters (no re-clustering) and adds what the
# cells are doing:
#   37. Immune composition: CIBERSORT LM22 in relative mode = share of each cell type
#       within the immune infiltrate (Newman et al. Nat Methods 2015;12:453-457)
#   38. Relationships between immune, stromal, functional and tumour signals
#   39. Functional and tumour-related signals (gene signatures, sources below)
#   40–41. Signals by cluster; basal-like "hot" vs "cold" tumours
#   42. Signals by PAM50 (other clinical associations: clinical_immune_integration.R;
#       histologic grade is not available in TCGA-BRCA)
#   43. Integrated TME figures
#
# Gene signatures (gene lists checked against the cited sources):
#   - Cytolytic activity: geometric mean of GZMA and PRF1 in TPM, 0.01 offset
#     (Rooney et al. Cell 2015;160:48-61)
#   - IFN-gamma response: 10-gene preliminary IFN-gamma set
#     (Ayers et al. J Clin Invest 2017;127:2930-2940)
#   - Antigen presentation: 18 antigen processing and presentation genes
#     (Wang et al. eLife 2019;8:e49020; gene list from Leone et al. 2013).
#     Wang et al. scored with GSVA; here the mean z-score is used (stated in methods).
#   - TGF-beta signalling and proliferation (G2M checkpoint): MSigDB Hallmark sets
#     (Liberzon et al. Cell Syst 2015;1:417-425), read with the msigdbr package
#   - Immune checkpoints: single genes (CD274/PD-L1, PDCD1/PD-1, PDCD1LG2/PD-L2,
#     CTLA4, LAG3, HAVCR2/TIM-3, TIGIT, IDO1)
# Signature score = mean of per-gene z-scores of log2(TPM + 1) across the analysis set.
# Statistics: Kruskal-Wallis (epsilon-squared), Wilcoxon rank-sum, Spearman;
# Benjamini-Hochberg FDR within each table.
# ============================================================

for (p in c("ggplot2", "pheatmap", "msigdbr"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(library(ggplot2))
set.seed(2026)

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
meta_dir <- file.path(project_dir, "metadata")
out_dir  <- file.path(project_dir, "results", "tme_characterization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nn <- function(s) tolower(gsub("[^A-Za-z0-9]", "", s))
eps2 <- function(H, n) as.numeric(H / (n - 1))
med_iqr <- function(x) sprintf("%.2f (%.2f to %.2f)", median(x, na.rm = TRUE),
                               quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))

## 1. DATA =======================================================
m    <- readRDS(file.path(meta_dir, "clinical_master_table.rds"))
land <- readRDS(file.path(project_dir, "results", "immune_landscape", "sample_level_results.rds"))
a <- m[m$in_analysis_set, ]
if (!setequal(a$sample_id, land$sample_id)) stop("Run scripts 1-3 first: analysis set and immune results differ.")
li <- match(a$sample_id, land$sample_id)
for (f in c("CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial")) a[[f]] <- land[[f]][li]
a$stroma_score <- rowMeans(scale(cbind(a$CAFs, a$Endothelial)))
lab_order <- unique(a$cluster_label[order(a$cluster)])
a$cluster_name <- factor(sub("^C[0-9]+: ", "", a$cluster_label), levels = sub("^C[0-9]+: ", "", lab_order))
cl_cols <- setNames(hcl.colors(nlevels(a$cluster_name), "Dark 3"), levels(a$cluster_name))
ids <- a$sample_id

tpm <- read.delim(file.path(project_dir, "data", "CIBERSORT", "BRCA_TPM_1106.txt"),
                  check.names = FALSE, stringsAsFactors = FALSE)
genes <- tpm[[1]]; tpm <- as.matrix(tpm[, -1]); storage.mode(tpm) <- "numeric"
rownames(tpm) <- genes
colnames(tpm) <- gsub(".", "-", colnames(tpm), fixed = TRUE)
if (!all(ids %in% colnames(tpm))) stop("Some analysis-set samples are missing from the expression matrix.")
tpm <- tpm[, ids]
lt  <- log2(tpm + 1)

## 2. GENE SETS ==================================================
hallmark <- NULL
hallmark <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                     error = function(e) tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
                                                  error = function(e2) NULL))
hall <- function(name) if (is.null(hallmark)) character(0) else unique(hallmark$gene_symbol[hallmark$gs_name == name])

# Fallback for TGF-beta if msigdbr is unavailable: MSigDB v7.0 Hallmark list (53 genes)
tgf_fallback <- c("APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1", "CDK9", "CDKN1C", "CTNNB1",
                  "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB",
                  "KLF10", "LEFTY2", "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA",
                  "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1", "SMAD1", "SMAD3",
                  "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1", "TGFB1", "TGFBR1", "TGIF1", "THBS1",
                  "TJP1", "TRIM33", "UBE2D3", "WWTR1", "XIAP")
tgf <- hall("HALLMARK_TGF_BETA_SIGNALING"); if (!length(tgf)) tgf <- tgf_fallback
g2m <- hall("HALLMARK_G2M_CHECKPOINT")
if (!length(g2m)) message("msigdbr unavailable: proliferation (G2M) score skipped. ",
                          "Install msigdbr (and msigdbdf if prompted) to include it.")

gene_sets <- list(
  "IFN-gamma response"   = c("IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "PRF1", "GZMA", "HLA-DRA"),
  "Antigen presentation" = c("PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10", "TAP1", "TAP2", "ERAP1",
                             "ERAP2", "CANX", "CALR", "PDIA3", "TAPBP", "B2M", "HLA-A", "HLA-B", "HLA-C"),
  "TGF-beta signalling"  = tgf,
  "Proliferation (G2M)"  = g2m)
gene_sets <- gene_sets[lengths(gene_sets) > 0]
set_source <- c("IFN-gamma response" = "Ayers et al. 2017 (10-gene preliminary IFN-gamma set)",
                "Antigen presentation" = "Wang et al. 2019 eLife (from Leone et al. 2013)",
                "TGF-beta signalling" = if (is.null(hallmark)) "MSigDB Hallmark v7.0 (fallback list)" else
                  paste("MSigDB Hallmark via msigdbr", as.character(packageVersion("msigdbr"))),
                "Proliferation (G2M)" = paste("MSigDB Hallmark via msigdbr",
                                              if (requireNamespace("msigdbr", quietly = TRUE)) as.character(packageVersion("msigdbr")) else ""))
checkpoints <- c("PD-L1 (CD274)" = "CD274", "PD-1 (PDCD1)" = "PDCD1", "PD-L2 (PDCD1LG2)" = "PDCD1LG2",
                 "CTLA4" = "CTLA4", "LAG3" = "LAG3", "TIM-3 (HAVCR2)" = "HAVCR2", "TIGIT" = "TIGIT",
                 "IDO1" = "IDO1")
checkpoints <- checkpoints[checkpoints %in% rownames(lt)]

score_mean_z <- function(g) {
  x <- lt[intersect(g, rownames(lt)), , drop = FALSE]
  x <- x[apply(x, 1, sd) > 0, , drop = FALSE]
  colMeans(t(scale(t(x))))
}

sig <- data.frame(row.names = ids, check.names = FALSE)
sig[["Cytolytic activity"]] <- log2(sqrt((tpm["GZMA", ] + 0.01) * (tpm["PRF1", ] + 0.01)))
for (s in names(gene_sets)) sig[[s]] <- score_mean_z(gene_sets[[s]])
for (cp in names(checkpoints)) sig[[cp]] <- lt[checkpoints[[cp]], ]
sig_names <- c("Cytolytic activity", names(gene_sets))
cp_names  <- names(checkpoints)

genes_used <- rbind(
  data.frame(signature = "Cytolytic activity", source = "Rooney et al. 2015 Cell",
             n_genes_listed = 2, n_genes_found = sum(c("GZMA", "PRF1") %in% rownames(lt)),
             genes_missing = paste(setdiff(c("GZMA", "PRF1"), rownames(lt)), collapse = ", "),
             genes = "GZMA, PRF1"),
  do.call(rbind, lapply(names(gene_sets), function(s) data.frame(
    signature = s, source = set_source[[s]], n_genes_listed = length(gene_sets[[s]]),
    n_genes_found = sum(gene_sets[[s]] %in% rownames(lt)),
    genes_missing = paste(setdiff(gene_sets[[s]], rownames(lt)), collapse = ", "),
    genes = paste(gene_sets[[s]], collapse = ", ")))))
write.csv(genes_used, file.path(out_dir, "signature_gene_sets_used.csv"), row.names = FALSE)

a <- cbind(a, sig[ids, , drop = FALSE])

## 3. STEP 37: IMMUNE COMPOSITION (CIBERSORT, RELATIVE) ===========
cib <- readRDS(file.path(project_dir, "results", "deconvolution_QC", "aligned_CIBERSORT_scores.rds"))
cib <- cib[match(ids, cib$sample_id), ]
cib_col <- function(lab) { i <- match(nn(lab), nn(names(cib))); if (is.na(i)) stop("CIBERSORT column not found: ", lab); names(cib)[i] }
lineages <- list(
  "CD8 T cells" = "T cells CD8",
  "CD4 T cells" = c("T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
                    "T cells follicular helper"),
  "Tregs" = "T cells regulatory (Tregs)", "Gamma-delta T cells" = "T cells gamma delta",
  "B cells" = c("B cells naive", "B cells memory"), "Plasma cells" = "Plasma cells",
  "NK cells" = c("NK cells resting", "NK cells activated"), "Monocytes" = "Monocytes",
  "Macrophages M0" = "Macrophages M0", "Macrophages M1" = "Macrophages M1", "Macrophages M2" = "Macrophages M2",
  "Dendritic cells" = c("Dendritic cells resting", "Dendritic cells activated"),
  "Mast cells" = c("Mast cells resting", "Mast cells activated"),
  "Eosinophils" = "Eosinophils", "Neutrophils" = "Neutrophils")
comp <- sapply(lineages, function(cols) rowSums(cib[, sapply(cols, cib_col), drop = FALSE]))
comp <- 100 * comp                                    # % of immune infiltrate
mac <- comp[, "Macrophages M0"] + comp[, "Macrophages M1"] + comp[, "Macrophages M2"]
a$M2_share_of_macrophages <- ifelse(mac > 0, 100 * comp[, "Macrophages M2"] / mac, NA)

cib_diag <- readRDS(file.path(project_dir, "results", "deconvolution_QC", "aligned_CIBERSORT_diagnostics.rds"))
cib_sig <- cib_diag[["P-value"]][match(ids, cib_diag$sample_id)] < 0.05

comp_tab <- do.call(rbind, lapply(colnames(comp), function(cn) {
  kw <- kruskal.test(comp[, cn] ~ a$cluster_name)
  data.frame(cell_type = cn,
             t(setNames(sprintf("%.1f", tapply(comp[, cn], a$cluster_name, mean)), paste0("mean % ", levels(a$cluster_name)))),
             t(setNames(sprintf("%.1f", tapply(comp[cib_sig, cn], a$cluster_name[cib_sig], mean)),
                        paste0("mean % (P<0.05) ", levels(a$cluster_name)))),
             epsilon_squared = round(eps2(kw$statistic, nrow(a)), 3), p_value = kw$p.value, check.names = FALSE)
}))
comp_tab$q_value_fdr <- signif(p.adjust(comp_tab$p_value, "BH"), 3); comp_tab$p_value <- signif(comp_tab$p_value, 3)
write.csv(comp_tab, file.path(out_dir, "immune_composition_by_cluster_cibersort.csv"), row.names = FALSE)

comp_long <- do.call(rbind, lapply(levels(a$cluster_name), function(cl)
  data.frame(cluster = cl, cell_type = colnames(comp), pct = colMeans(comp[a$cluster_name == cl, , drop = FALSE]))))
comp_long$cluster <- factor(comp_long$cluster, levels = levels(a$cluster_name))
comp_long$cell_type <- factor(comp_long$cell_type, levels = rev(colnames(comp)))
p_comp <- ggplot(comp_long, aes(cluster, pct, fill = cell_type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = setNames(hcl.colors(ncol(comp), "Set 3"), colnames(comp))) +
  labs(x = NULL, y = "Mean share of immune infiltrate (%)", fill = "Cell type",
       title = "Immune composition by cluster (CIBERSORT, relative)",
       caption = "Shares of the immune infiltrate, not of the whole tumour") +
  theme_bw()
ggsave(file.path(out_dir, "immune_composition_by_cluster_barplot.png"), p_comp, width = 8, height = 6, dpi = 200)

## 4. STEP 38: RELATIONSHIPS =====================================
rel_vars <- c("immune_score", "stroma_score", "CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial",
              sig_names, cp_names, "leukocyte_fraction_methylation", "absolute_purity")
rel_lab <- c(immune_score = "Immune score", stroma_score = "Stromal score",
             leukocyte_fraction_methylation = "Leukocyte fraction (methylation)",
             absolute_purity = "Tumour purity (ABSOLUTE)")
rel <- sapply(rel_vars, function(v) a[[v]])
colnames(rel) <- ifelse(colnames(rel) %in% names(rel_lab), rel_lab[colnames(rel)], colnames(rel))
cm <- cor(rel, method = "spearman", use = "pairwise.complete.obs")
write.csv(round(cm, 3), file.path(out_dir, "tme_feature_spearman_correlations.csv"))
pheatmap::pheatmap(cm, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6,
                   breaks = seq(-1, 1, length.out = 101),
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   main = "Relationships between TME features (Spearman)",
                   filename = file.path(out_dir, "tme_feature_correlation_heatmap.png"), width = 11, height = 10)

## 5. STEPS 39–41: SIGNALS BY CLUSTER ============================
feat <- c(sig_names, cp_names, "M2_share_of_macrophages")
by_cluster <- do.call(rbind, lapply(feat, function(f) {
  ok <- !is.na(a[[f]])
  kw <- kruskal.test(a[[f]][ok] ~ a$cluster_name[ok])
  data.frame(feature = f, t(tapply(a[[f]][ok], a$cluster_name[ok], med_iqr)),
             epsilon_squared = round(eps2(kw$statistic, sum(ok)), 3), p_value = kw$p.value, check.names = FALSE)
}))
by_cluster$q_value_fdr <- signif(p.adjust(by_cluster$p_value, "BH"), 3); by_cluster$p_value <- signif(by_cluster$p_value, 3)
write.csv(by_cluster, file.path(out_dir, "functional_signals_by_cluster.csv"), row.names = FALSE)

blong <- do.call(rbind, lapply(c(sig_names, cp_names), function(f)
  data.frame(feature = f, cluster = a$cluster_name, value = a[[f]])))
blong$feature <- factor(blong$feature, levels = c(sig_names, cp_names))
p_box <- ggplot(blong, aes(cluster, value, fill = cluster)) +
  geom_boxplot(outlier.size = 0.3) + facet_wrap(~ feature, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = cl_cols) +
  labs(x = NULL, y = "Score (signatures: mean z; checkpoints and CYT: log2 expression)",
       title = "Functional and tumour-related signals by immune cluster") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, "functional_signals_by_cluster_boxplot.png"), p_box, width = 13, height = 8, dpi = 200)

# Functional signals relative to immune content: high cytotoxicity per unit of immune infiltrate?
resid_tab <- do.call(rbind, lapply(c(sig_names, cp_names), function(f) {
  r <- resid(lm(a[[f]] ~ a$immune_score))
  kw <- kruskal.test(r ~ a$cluster_name)
  data.frame(feature = f, t(round(tapply(r, a$cluster_name, median), 3)),
             epsilon_squared = round(eps2(kw$statistic, nrow(a)), 3), p_value = kw$p.value, check.names = FALSE)
}))
resid_tab$q_value_fdr <- signif(p.adjust(resid_tab$p_value, "BH"), 3); resid_tab$p_value <- signif(resid_tab$p_value, 3)
write.csv(resid_tab, file.path(out_dir, "functional_signals_adjusted_for_immune_score.csv"), row.names = FALSE)

## 6. BASAL-LIKE: HOT vs COLD ====================================
bas <- a$pam50 %in% "Basal-like" & a$cluster_name %in% c("Immune-enriched", "Depleted")
b <- a[bas, ]; b$group <- droplevels(b$cluster_name)
std <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)   # z within the whole analysis set
hc_feat <- c("immune_score", "stroma_score", sig_names, cp_names, "M2_share_of_macrophages",
             "leukocyte_fraction_methylation", "absolute_purity", "age")
hot_cold <- do.call(rbind, lapply(hc_feat, function(f) {
  z <- std(a[[f]])[bas]
  g1 <- z[b$group == "Immune-enriched"]; g2 <- z[b$group == "Depleted"]
  wt <- suppressWarnings(wilcox.test(g1, g2, exact = FALSE))
  data.frame(feature = ifelse(f %in% names(rel_lab), rel_lab[f], f),
             n_hot = sum(!is.na(g1)), n_cold = sum(!is.na(g2)),
             median_z_hot = round(median(g1, na.rm = TRUE), 2), median_z_cold = round(median(g2, na.rm = TRUE), 2),
             difference = round(median(g1, na.rm = TRUE) - median(g2, na.rm = TRUE), 2), p_value = wt$p.value)
}))
hot_cold$q_value_fdr <- signif(p.adjust(hot_cold$p_value, "BH"), 3); hot_cold$p_value <- signif(hot_cold$p_value, 3)
write.csv(hot_cold, file.path(out_dir, "basal_like_hot_vs_cold.csv"), row.names = FALSE)

write.csv(as.data.frame.matrix(table(b$group, b$receptor_subtype, useNA = "ifany")),
          file.path(out_dir, "basal_like_hot_vs_cold_receptor_subtype.csv"))
write.csv(as.data.frame.matrix(table(b$group, b$stage, useNA = "ifany")),
          file.path(out_dir, "basal_like_hot_vs_cold_stage.csv"))

hc_plot <- hot_cold
hc_plot$feature <- factor(hc_plot$feature, levels = hc_plot$feature[order(hc_plot$difference)])
hc_plot$sig <- ifelse(hc_plot$q_value_fdr < 0.05, "FDR < 0.05", "not significant")
p_hc <- ggplot(hc_plot, aes(difference, feature, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey60") + geom_point(size = 3) +
  scale_colour_manual(values = c("FDR < 0.05" = "#B2182B", "not significant" = "grey60")) +
  labs(x = "Median difference, hot minus cold (z units)", y = NULL, colour = NULL,
       title = sprintf("Basal-like tumours: Immune-enriched (n = %d) vs Depleted (n = %d)",
                       sum(b$group == "Immune-enriched"), sum(b$group == "Depleted"))) +
  theme_bw()
ggsave(file.path(out_dir, "basal_like_hot_vs_cold.png"), p_hc, width = 8.5, height = 6.5, dpi = 200)

## 7. STEP 42: SIGNALS BY PAM50 ==================================
by_pam <- do.call(rbind, lapply(c(sig_names, cp_names), function(f) {
  ok <- !is.na(a$pam50)
  kw <- kruskal.test(a[[f]][ok] ~ droplevels(a$pam50[ok]))
  data.frame(feature = f, t(round(tapply(a[[f]][ok], droplevels(a$pam50[ok]), median), 2)),
             epsilon_squared = round(eps2(kw$statistic, sum(ok)), 3), p_value = kw$p.value, check.names = FALSE)
}))
by_pam$q_value_fdr <- signif(p.adjust(by_pam$p_value, "BH"), 3); by_pam$p_value <- signif(by_pam$p_value, 3)
write.csv(by_pam, file.path(out_dir, "functional_signals_by_PAM50.csv"), row.names = FALSE)

## 8. STEP 43: INTEGRATED TME FIGURES ============================
zfeat <- c("CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial", sig_names, cp_names)
Z <- sapply(zfeat, function(f) std(a[[f]]))
rownames(Z) <- ids
colnames(Z) <- gsub("_", " ", colnames(Z))

# 8a. Cluster-level summary
cl_mean <- t(sapply(levels(a$cluster_name), function(cl) colMeans(Z[a$cluster_name == cl, , drop = FALSE])))
rownames(cl_mean) <- paste0(levels(a$cluster_name), " (n=", table(a$cluster_name), ")")
pheatmap::pheatmap(t(cl_mean), cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.2f",
                   breaks = seq(-1.2, 1.2, length.out = 101),
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   gaps_row = c(5, 5 + length(sig_names)),
                   main = "TME profile by cluster (mean z-score)",
                   filename = file.path(out_dir, "tme_profile_by_cluster_heatmap.png"), width = 6.5, height = 7.5)

# 8b. Tumour-level integrated heatmap
o <- order(a$cluster_name, -a$immune_score)
hm <- t(Z[o, ]); hm[hm > 2.5] <- 2.5; hm[hm < -2.5] <- -2.5
pam_cols <- c("Luminal A" = "#1F4E9C", "Luminal B" = "#7FB2E5", "HER2-enriched" = "#E78AC3",
              "Basal-like" = "#8C2D04", "Normal-like" = "#33A02C")
ann <- data.frame(Cluster = a$cluster_name[o], PAM50 = droplevels(a$pam50[o]),
                  `Leukocyte fraction` = a$leukocyte_fraction_methylation[o],
                  Purity = a$absolute_purity[o], row.names = ids[o], check.names = FALSE)
ann_cols <- list(Cluster = cl_cols, PAM50 = pam_cols[levels(ann$PAM50)],
                 `Leukocyte fraction` = colorRampPalette(c("white", "#54278F"))(50),
                 Purity = colorRampPalette(c("white", "#252525"))(50))
pheatmap::pheatmap(hm, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
                   annotation_col = ann, annotation_colors = ann_cols,
                   gaps_col = cumsum(table(a$cluster_name[o]))[-nlevels(a$cluster_name)],
                   gaps_row = c(5, 5 + length(sig_names)),
                   breaks = seq(-2.5, 2.5, length.out = 101),
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   main = sprintf("Integrated tumour microenvironment (n = %d)", nrow(a)),
                   filename = file.path(out_dir, "integrated_tme_heatmap.png"), width = 13, height = 7.5)

out_scores <- data.frame(sample_id = ids, cluster = a$cluster_name, pam50 = a$pam50,
                         a[, c(sig_names, cp_names, "M2_share_of_macrophages")], comp, check.names = FALSE)
write.csv(out_scores, file.path(out_dir, "tme_sample_level_scores.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_session_info.txt"))

## 9. REPORT ====================================================
op <- options(width = 220)
cat("\n================ TME CHARACTERISATION ================\n")
cat("Tumours:", nrow(a), "| CIBERSORT P<0.05 subset:", sum(cib_sig, na.rm = TRUE), "\n")
cat("\nGene sets used:\n"); print(genes_used[, c("signature", "source", "n_genes_listed", "n_genes_found", "genes_missing")], row.names = FALSE)
cat("\n37. Immune composition by cluster (CIBERSORT, mean % of immune infiltrate):\n"); print(comp_tab, row.names = FALSE)
cat("\n39-41. Functional signals by cluster, median (IQR):\n"); print(by_cluster, row.names = FALSE)
cat("\nSame signals after removing the effect of overall immune content (median residual):\n"); print(resid_tab, row.names = FALSE)
cat("\nBasal-like hot vs cold:\n"); print(hot_cold, row.names = FALSE)
cat("\nBasal-like hot vs cold, receptor subtype:\n"); print(table(b$group, b$receptor_subtype, useNA = "ifany"))
cat("\n42. Signals by PAM50 (median):\n"); print(by_pam, row.names = FALSE)
cat("\n38. Key correlations (Spearman):\n")
key <- c("Immune score", "Stromal score", sig_names, "Tumour purity (ABSOLUTE)")
print(round(cm[key, key], 2))
cat("\nSaved to:", out_dir, "\n"); print(sort(list.files(out_dir)))
options(op)