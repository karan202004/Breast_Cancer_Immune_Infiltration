#Immune landscape

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
qc_dir   <- file.path(project_dir, "results", "deconvolution_QC")
cmp_dir  <- file.path(project_dir, "results", "method_comparison_consensus")
meta_dir <- file.path(project_dir, "metadata")
out_dir  <- file.path(project_dir, "results", "immune_landscape")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

# Features used for clustering (all five consensus features = immune + stroma).
# For immune-only clustering use: c("CD8_T_cells", "B_cells", "Macrophages")
cluster_features <- c("CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial")
immune_features  <- c("CD8_T_cells", "B_cells", "Macrophages")
stroma_features  <- c("CAFs", "Endothelial")

k_range       <- 2:6    # numbers of clusters tested
n_boot        <- 100    # subsampling runs for stability
boot_frac     <- 0.8    # fraction of samples per run
min_stability <- 0.75   # minimum mean adjusted Rand index to accept a k
k_override    <- 3      # k = 3 chosen after review: most stable (ARI 0.967) and separates
# immune and stromal axes; k = 2 mainly reflected tumour purity

pkgs <- c("ggplot2", "pheatmap", "cluster", "mclust")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(ggplot2); library(mclust) })

nn  <- function(s) tolower(gsub("[^A-Za-z0-9]", "", s))
int <- function(v) qnorm((rank(v, ties.method = "average") - 0.5) / length(v))

#Load and align data
cons <- readRDS(file.path(cmp_dir, "06_consensus_immune_features.rds"))
clin <- readRDS(file.path(meta_dir, "clinical_variables.rds"))   # from clinical_data_and_analysis_set.R
if (!identical(as.character(cons$sample_id), as.character(clin$sample_id)))
  stop("Consensus and clinical files are not in the same sample order.")

read_aligned <- function(m) {
  x <- readRDS(file.path(qc_dir, paste0("aligned_", m, "_scores.rds")))
  if (!identical(as.character(x$sample_id), as.character(cons$sample_id)))
    stop(m, ": sample order differs from consensus file.")
  x
}
q   <- read_aligned("quantiseq")
ep  <- read_aligned("EPIC")
mcp <- read_aligned("MCPcounter")
xc  <- read_aligned("xCell")

getcol <- function(df, name) {
  i <- match(nn(name), nn(names(df)))
  if (is.na(i)) stop("Column not found: ", name)
  as.numeric(df[[i]])
}

#Analysis set
n_all    <- nrow(clin)
analysis <- clin$in_analysis_set
if (is.null(analysis)) stop("Run clinical_data_and_analysis_set.R first.")

reason_tab <- table(clin$exclusion_reason)
exclusions <- data.frame(
  reason = c("Total samples", names(reason_tab), "Analysis set"),
  n = c(n_all, as.integer(reason_tab), sum(analysis)))
write.csv(exclusions, file.path(out_dir, "sample_exclusions.csv"), row.names = FALSE)

ids <- as.character(clin$sample_id[analysis])
cl  <- clin[analysis, ]
n   <- length(ids)

#Feature table
df <- data.frame(
  CD8_T_cells = cons$CD8_T_cells, B_cells = cons$B_cells, Macrophages = cons$Macrophages,
  CAFs = cons$CAFs, Endothelial = cons$Endothelial,
  Tregs_quanTIseq       = getcol(q,   "T cell regulatory (Tregs)"),
  Tregs_xCell           = getcol(xc,  "T cell regulatory (Tregs)"),
  MyeloidDC_MCPcounter  = getcol(mcp, "Myeloid dendritic cell"),
  MyeloidDC_xCell       = getcol(xc,  "Myeloid dendritic cell"),
  Monocytes_MCPcounter  = getcol(mcp, "Monocyte"),
  Monocytes_xCell       = getcol(xc,  "Monocyte"),
  ImmuneScore_xCell     = getcol(xc,  "immune score"),
  StromaScore_xCell     = getcol(xc,  "stroma score"),
  check.names = FALSE)[analysis, ]
rownames(df) <- ids

na_feat <- colSums(is.na(df))
if (any(na_feat > 0)) stop("Missing values in features: ",
                           paste(names(na_feat)[na_feat > 0], collapse = ", "))

overall_cols <- c("ImmuneScore_xCell", "StromaScore_xCell")
consensus_cols <- c("CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial")
tier <- setNames(ifelse(names(df) %in% consensus_cols, "Consensus",
                        ifelse(names(df) %in% overall_cols, "Overall score", "Exploratory")),
                 names(df))

# Put every feature on a comparable rank-normal scale for plots/heatmaps
df_int <- as.data.frame(lapply(df, int), check.names = FALSE)
rownames(df_int) <- ids

dist_summary <- do.call(rbind, lapply(names(df), function(f) {
  v <- df[[f]]
  data.frame(feature = f, tier = tier[[f]], n = length(v),
             median = round(median(v), 4), q1 = round(quantile(v, .25), 4),
             q3 = round(quantile(v, .75), 4), min = round(min(v), 4),
             max = round(max(v), 4), zero_fraction = round(mean(v == 0), 3),
             row.names = NULL)
}))
write.csv(dist_summary, file.path(out_dir, "feature_distribution_summary.csv"), row.names = FALSE)

long <- do.call(rbind, lapply(names(df_int), function(f)
  data.frame(feature = f, tier = tier[[f]], value = df_int[[f]])))
long$feature <- factor(long$feature, levels = names(df))
long$tier <- factor(long$tier, levels = c("Consensus", "Exploratory", "Overall score"))

p24a <- ggplot(long, aes(feature, value, fill = tier)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_grid(. ~ tier, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Relative level (rank-normal scale)",
       title = "Immune and stromal features across tumours") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, "feature_distribution_violin_plot.png"), p24a, width = 11, height = 5, dpi = 200)

#Interpretable fractions (quanTIseq and EPIC estimate % of all cells)
frac_tier <- function(method, cell) {
  c0 <- tolower(cell)
  if (grepl("uncharacterized", c0)) return("Residual (non-immune/tumour)")
  if (method == "EPIC" && grepl("cd8", c0)) return("Unreliable")
  if (grepl("cd8|^b cell|macrophage|fibroblast|endothelial", c0)) return("Consensus")
  if (grepl("treg|monocyte|dendritic", c0)) return("Exploratory")
  "Unreliable"
}
frac_long <- do.call(rbind, c(
  lapply(setdiff(names(q), "sample_id"), function(cn)
    data.frame(method = "quanTIseq", cell = cn, tier = frac_tier("quanTIseq", cn),
               pct = 100 * as.numeric(q[[cn]])[analysis])),
  lapply(setdiff(names(ep), "sample_id"), function(cn)
    data.frame(method = "EPIC", cell = cn, tier = frac_tier("EPIC", cn),
               pct = 100 * as.numeric(ep[[cn]])[analysis]))))

frac_summary <- aggregate(pct ~ method + cell + tier, data = frac_long,
                          FUN = function(v) round(c(median = median(v), q1 = quantile(v, .25),
                                                    q3 = quantile(v, .75)), 2))
frac_summary <- do.call(data.frame, frac_summary)
names(frac_summary) <- c("method", "cell", "tier", "median_pct", "q1_pct", "q3_pct")
write.csv(frac_summary, file.path(out_dir, "cell_fractions_quanTIseq_EPIC_summary.csv"), row.names = FALSE)

p24b <- ggplot(subset(frac_long, !grepl("Residual", tier)),
               aes(reorder(cell, pct, median), pct, fill = tier)) +
  geom_boxplot(outlier.size = 0.3) + coord_flip() +
  facet_wrap(~ method, scales = "free_y") +
  scale_fill_manual(values = c(Consensus = "#2B8CBE", Exploratory = "#FDB863", Unreliable = "grey75")) +
  labs(x = NULL, y = "% of all cells in tumour", fill = "Reliability",
       title = "Estimated cell fractions (non-immune/tumour fraction not shown)") +
  theme_bw()
ggsave(file.path(out_dir, "cell_fractions_quanTIseq_EPIC_boxplot.png"), p24b, width = 11, height = 5, dpi = 200)

##Correlations
cm <- cor(df, method = "spearman")
write.csv(round(cm, 3), file.path(out_dir, "feature_spearman_correlation_matrix.csv"))
pheatmap::pheatmap(
  cm, display_numbers = TRUE, number_format = "%.2f", fontsize_number = 7,
  breaks = seq(-1, 1, length.out = 101),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Spearman correlation between immune/stromal features",
  filename = file.path(out_dir, "feature_spearman_correlation_heatmap.png"), width = 9, height = 8)

#Clusturing
X <- scale(as.matrix(df[, cluster_features]))
D <- dist(X)

ari <- function(a, b) {
  tb <- table(a, b); s <- function(x) sum(choose(x, 2))
  sij <- s(tb); sa <- s(rowSums(tb)); sb <- s(colSums(tb))
  e <- sa * sb / choose(sum(tb), 2); mx <- (sa + sb) / 2
  (sij - e) / (mx - e)
}
run_km <- function(Xm, k, nstart = 50) kmeans(Xm, centers = k, nstart = nstart, iter.max = 100)$cluster

k_sel <- do.call(rbind, lapply(k_range, function(k) {
  full <- run_km(X, k)
  sil  <- mean(cluster::silhouette(full, D)[, 3])
  stab <- replicate(n_boot, {
    idx <- sample(n, floor(boot_frac * n))
    ari(full[idx], run_km(X[idx, , drop = FALSE], k, nstart = 20))
  })
  data.frame(k = k, mean_silhouette = round(sil, 3),
             stability_ARI = round(mean(stab), 3),
             smallest_cluster = min(table(full)))
}))

ok <- k_sel$stability_ARI >= min_stability
k_best <- if (!is.na(k_override)) k_override else
  if (any(ok)) k_sel$k[ok][which.max(k_sel$mean_silhouette[ok])] else
    k_sel$k[which.max(k_sel$mean_silhouette)]
k_sel$chosen <- k_sel$k == k_best
write.csv(k_sel, file.path(out_dir, "cluster_number_selection.csv"), row.names = FALSE)

png(file.path(out_dir, "cluster_number_selection_plot.png"), width = 1400, height = 600, res = 150)
par(mfrow = c(1, 2))
plot(k_sel$k, k_sel$mean_silhouette, type = "b", pch = 19, xlab = "Number of clusters (k)",
     ylab = "Mean silhouette", main = "Separation"); abline(v = k_best, lty = 2, col = "red")
plot(k_sel$k, k_sel$stability_ARI, type = "b", pch = 19, xlab = "Number of clusters (k)",
     ylab = "Mean ARI (subsampling)", main = "Stability", ylim = c(0, 1))
abline(h = min_stability, lty = 3); abline(v = k_best, lty = 2, col = "red")
dev.off()

raw_cl <- run_km(X, k_best, nstart = 100)

#Order clusters from most to least immune-rich
imm_by_cl <- tapply(rowMeans(X[, immune_features, drop = FALSE]), raw_cl, mean)
cl_order  <- order(imm_by_cl, decreasing = TRUE)
cluster   <- factor(paste0("C", match(raw_cl, as.integer(names(imm_by_cl))[cl_order])),
                    levels = paste0("C", seq_len(k_best)))

#Phenotypes
profile <- t(sapply(levels(cluster), function(c) colMeans(X[cluster == c, , drop = FALSE])))
lvl <- function(z) ifelse(z > 0.5, "high", ifelse(z < -0.5, "low", "intermediate"))
imm_z <- rowMeans(profile[, intersect(immune_features, colnames(profile)), drop = FALSE])
str_z <- if (all(stroma_features %in% colnames(profile)))
  rowMeans(profile[, stroma_features, drop = FALSE]) else rep(NA, nrow(profile))
cl_label <- paste0(rownames(profile), ": Immune-", lvl(imm_z),
                   if (!all(is.na(str_z))) paste0(" / Stroma-", lvl(str_z)) else "")
names(cl_label) <- rownames(profile)

# Final names for k = 3, assigned from each cluster's profile (not its number):
# highest immune -> Immune-enriched; of the rest, highest stroma -> Stroma-rich; last -> Depleted
if (k_best == 3) {
  final <- setNames(character(3), rownames(profile))
  i_imm <- which.max(imm_z)
  rest  <- setdiff(seq_len(3), i_imm)
  i_str <- rest[which.max(str_z[rest])]
  final[i_imm] <- "Immune-enriched"
  final[i_str] <- "Stroma-rich"
  final[setdiff(rest, i_str)] <- "Depleted"
  cl_label <- paste0(rownames(profile), ": ", final)
  names(cl_label) <- rownames(profile)
}

prof_out <- data.frame(cluster = rownames(profile), label = cl_label,
                       n = as.integer(table(cluster)[rownames(profile)]),
                       round(profile, 2), immune_z = round(imm_z, 2), stroma_z = round(str_z, 2),
                       check.names = FALSE)
write.csv(prof_out, file.path(out_dir, "cluster_profiles.csv"), row.names = FALSE)

pm <- profile; rownames(pm) <- paste0(cl_label, " (n=", table(cluster)[rownames(profile)], ")")
pheatmap::pheatmap(
  pm, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f",
  breaks = seq(-1.5, 1.5, length.out = 101),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Cluster profiles (mean standardised score)",
  filename = file.path(out_dir, "cluster_profiles_heatmap.png"), width = 8, height = 1.2 + 0.6 * k_best)

# Clinical associations
h <- cl$histology
hist_grp <- ifelse(is.na(h), NA,
                   ifelse(grepl("lobular", h, TRUE) & !grepl("duct", h, TRUE), "Lobular",
                          ifelse(grepl("duct", h, TRUE) & !grepl("lobular", h, TRUE), "Ductal", "Other/mixed")))

xtab <- function(group, var) {
  ok <- !is.na(group) & !is.na(var)
  g <- droplevels(as.factor(group[ok])); v <- droplevels(as.factor(var[ok]))
  tb <- table(g, v)
  pct <- round(100 * prop.table(tb, 1), 1)
  p <- suppressWarnings(chisq.test(tb, simulate.p.value = TRUE, B = 10000)$p.value)
  res <- data.frame(group = rownames(tb),
                    setNames(as.data.frame.matrix(tb), paste0(colnames(tb), " (n)")),
                    setNames(as.data.frame.matrix(pct), paste0(colnames(pct), " (%)")),
                    check.names = FALSE)
  list(table = res, p = p)
}

assoc <- list(PAM50 = xtab(cluster, cl$pam50), Stage = xtab(cluster, cl$stage),
              Histology = xtab(cluster, hist_grp))
for (a in names(assoc))
  write.csv(assoc[[a]]$table, file.path(out_dir, paste0("cluster_vs_", a, ".csv")), row.names = FALSE)
age_p <- kruskal.test(cl$age ~ cluster)$p.value
assoc_p <- data.frame(variable = c(names(assoc), "Age"),
                      test = c(rep("Chi-square (simulated p, B=10000)", 3), "Kruskal-Wallis"),
                      p_value = signif(c(sapply(assoc, `[[`, "p"), age_p), 3))
write.csv(assoc_p, file.path(out_dir, "cluster_clinical_association_tests.csv"), row.names = FALSE)

pam_df <- data.frame(cluster = factor(cl_label[as.character(cluster)], levels = cl_label), pam50 = cl$pam50)
pam_df <- pam_df[!is.na(pam_df$pam50), ]
pam_cols <- c("Luminal A" = "#1F4E9C", "Luminal B" = "#7FB2E5", "HER2-enriched" = "#E78AC3",
              "Basal-like" = "#D7301F", "Normal-like" = "#33A02C", "Unknown" = "grey85")
p27 <- ggplot(pam_df, aes(cluster, fill = pam50)) +
  geom_bar(position = "fill") + scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  scale_fill_manual(values = pam_cols) +
  labs(x = NULL, y = "Tumours", fill = "PAM50", title = "PAM50 subtype composition by cluster") +
  theme_bw()
ggsave(file.path(out_dir, "cluster_PAM50_composition_barplot.png"), p27, width = 7, height = 5, dpi = 200)

# Consensus features by PAM50 subtype
by_pam <- do.call(rbind, lapply(cluster_features, function(f) {
  ok <- !is.na(cl$pam50)
  med <- tapply(df[[f]][ok], droplevels(cl$pam50[ok]), median)
  data.frame(feature = f, t(round(med, 3)),
             kruskal_p = signif(kruskal.test(df[[f]][ok] ~ droplevels(cl$pam50[ok]))$p.value, 3),
             check.names = FALSE)
}))
write.csv(by_pam, file.path(out_dir, "consensus_features_by_PAM50.csv"), row.names = FALSE)

pam_long <- do.call(rbind, lapply(cluster_features, function(f)
  data.frame(feature = f, pam50 = cl$pam50, value = df[[f]])))
pam_long <- pam_long[!is.na(pam_long$pam50), ]
pam_long$feature <- factor(pam_long$feature, levels = cluster_features)
p27b <- ggplot(pam_long, aes(pam50, value, fill = pam50)) +
  geom_boxplot(outlier.size = 0.3) + facet_wrap(~ feature, nrow = 1) +
  scale_fill_manual(values = pam_cols) +
  labs(x = NULL, y = "Consensus score", title = "Consensus features by PAM50 subtype") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, "consensus_features_by_PAM50_boxplot.png"), p27b, width = 13, height = 4.5, dpi = 200)

#Immune high/low
# Mean of standardised CD8, B-cell and macrophage scores.
immune_score <- rowMeans(scale(as.matrix(df[, immune_features])))

mc <- Mclust(immune_score, G = 1:3, modelNames = "V", verbose = FALSE)
if (is.null(mc)) mc <- Mclust(immune_score, G = 1:3, verbose = FALSE)

if (!is.null(mc) && mc$G >= 2) {
  comp_order <- order(mc$parameters$mean)
  grp_names  <- if (mc$G == 2) c("Immune-low", "Immune-high") else
    c("Immune-low", "Immune-intermediate", "Immune-high")
  immune_group <- factor(grp_names[match(mc$classification, comp_order)], levels = grp_names)
  grp_method <- sprintf("Gaussian mixture on immune score: %d components (BIC-selected)", mc$G)
} else {
  immune_group <- factor(paste0("Immune-", lvl(imm_z))[match(as.character(cluster), rownames(profile))])
  grp_method <- "Mixture model found one component (no natural immune-high/low split); groups taken from each cluster's immune level"
}

grp_summary <- data.frame(
  group = levels(immune_group), n = as.integer(table(immune_group)),
  immune_score_min = round(tapply(immune_score, immune_group, min), 2),
  immune_score_max = round(tapply(immune_score, immune_group, max), 2))
write.csv(grp_summary, file.path(out_dir, "immune_group_summary.csv"), row.names = FALSE)

grp_assoc <- list(PAM50 = xtab(immune_group, cl$pam50), Cluster = xtab(immune_group, cluster),
                  Stage = xtab(immune_group, cl$stage))
for (a in names(grp_assoc))
  write.csv(grp_assoc[[a]]$table, file.path(out_dir, paste0("immune_group_vs_", a, ".csv")), row.names = FALSE)

grp_cols <- setNames(colorRampPalette(c("#4575B4", "#D73027"))(nlevels(immune_group)), levels(immune_group))
p28 <- ggplot(data.frame(s = immune_score, g = immune_group), aes(s, fill = g)) +
  geom_histogram(bins = 50, colour = NA) + scale_fill_manual(values = grp_cols) +
  labs(x = "Immune score (mean standardised CD8, B-cell and macrophage scores)", y = "Tumours",
       fill = NULL, title = "Data-driven immune groups") + theme_bw()
ggsave(file.path(out_dir, "immune_group_histogram.png"), p28, width = 8, height = 4.5, dpi = 200)

#Landscape heatmap + Sample result
o <- order(cluster, -immune_score)
hm <- t(as.matrix(df_int[o, ]))
hm[hm > 2.5] <- 2.5; hm[hm < -2.5] <- -2.5

na_unknown <- function(x) { x <- as.character(x); x[is.na(x)] <- "Unknown"; x }
ann <- data.frame(
  Cluster = factor(cl_label[as.character(cluster[o])], levels = cl_label),
  Immune  = factor(as.character(immune_group[o]), levels = levels(immune_group)),
  PAM50   = factor(na_unknown(cl$pam50[o]), levels = names(pam_cols)),
  Stage   = factor(na_unknown(cl$stage[o]), levels = c("I", "II", "III", "IV", "Unknown")),
  row.names = colnames(hm))
ann_cols <- list(
  Cluster = setNames(hcl.colors(k_best, "Dark 3"), cl_label),
  Immune  = grp_cols,
  PAM50   = pam_cols[levels(droplevels(ann$PAM50))],
  Stage   = c(I = "#FEE5D9", II = "#FCAE91", III = "#FB6A4A", IV = "#CB181D", Unknown = "grey85")
  [levels(droplevels(ann$Stage))])
ann$PAM50 <- droplevels(ann$PAM50); ann$Stage <- droplevels(ann$Stage)

pheatmap::pheatmap(
  hm, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
  annotation_col = ann, annotation_colors = ann_cols,
  gaps_col = cumsum(table(cluster[o]))[-k_best], gaps_row = c(5, 11),
  breaks = seq(-2.5, 2.5, length.out = 101),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = sprintf("TCGA-BRCA immune landscape (n = %d)", n),
  filename = file.path(out_dir, "immune_landscape_heatmap.png"), width = 13, height = 6.5)

results <- data.frame(sample_id = ids, patient = cl$patient,
                      cluster = cluster, cluster_label = cl_label[as.character(cluster)],
                      immune_score = round(immune_score, 4), immune_group = immune_group,
                      pam50 = cl$pam50, stage = cl$stage, age = cl$age, histology_group = hist_grp,
                      df, check.names = FALSE)
write.csv(results, file.path(out_dir, "sample_level_results.csv"), row.names = FALSE)
saveRDS(results, file.path(out_dir, "sample_level_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_session_info.txt"))

cat("\nAnalysis set:\n"); print(exclusions, row.names = FALSE)

cat("\nconsensus feature distributions saved; cell fractions (median %):\n")
print(subset(frac_summary, !grepl("Residual", tier))[, c("method", "cell", "tier", "median_pct")],
      row.names = FALSE)

cat("\ncorrelations among consensus features (Spearman):\n")
print(round(cm[cluster_features, cluster_features], 2))

cat("\nchoosing the number of clusters:\n"); print(k_sel, row.names = FALSE)
cat("Chosen k:", k_best, if (!is.na(k_override)) "(manual override)" else "", "\n")

cat("\ncluster profiles (mean standardised score):\n")
print(prof_out, row.names = FALSE)
cat("\nPAM50 composition by cluster (% of cluster):\n")
print(assoc$PAM50$table[, c(1, grep("\\(%\\)", names(assoc$PAM50$table)))], row.names = FALSE)
cat("\nClinical association tests:\n"); print(assoc_p, row.names = FALSE)
cat("\nConsensus features by PAM50 (median):\n"); print(by_pam, row.names = FALSE)

cat("\n grp_method, "\n"); print(grp_summary, row.names = FALSE)
cat("\nPAM50 by immune group (% of group):\n")
print(grp_assoc$PAM50$table[, c(1, grep("\\(%\\)", names(grp_assoc$PAM50$table)))], row.names = FALSE)

cat("\n figures and per-sample results saved to:\n", out_dir, "\n")
cat("\nFiles created:\n"); print(sort(list.files(out_dir)))
