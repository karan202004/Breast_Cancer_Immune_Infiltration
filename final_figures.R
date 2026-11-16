# Final Figures

suppressPackageStartupMessages({
  library(survival); library(splines); library(ggplot2)
})
for (p in c("patchwork", "pheatmap")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(patchwork)

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
res_dir  <- file.path(project_dir, "results")
main_dir <- file.path(res_dir, "figures", "main")
supp_dir <- file.path(res_dir, "figures", "supplementary")
for (d in c(main_dir, supp_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


m    <- readRDS(file.path(project_dir, "metadata", "clinical_master_table.rds"))
land <- readRDS(file.path(res_dir, "immune_landscape", "sample_level_results.rds"))
tme  <- read.csv(file.path(res_dir, "tme_characterization", "tme_sample_level_scores.csv"), check.names = FALSE)
mv   <- read.csv(file.path(res_dir, "survival_analysis", "multivariable_cox_results.csv"), check.names = FALSE)

a <- m[m$in_analysis_set, ]
if (!setequal(a$sample_id, land$sample_id) || !setequal(a$sample_id, tme$sample_id))
  stop("Saved results do not cover the same tumours. Re-run scripts 1-7 in order.")
land <- land[match(a$sample_id, land$sample_id), ]
tme  <- tme[match(a$sample_id, tme$sample_id), ]

landscape_features <- c("CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial",
                        "ImmuneScore_xCell", "StromaScore_xCell",
                        "Tregs_quanTIseq", "Tregs_xCell", "MyeloidDC_MCPcounter", "MyeloidDC_xCell",
                        "Monocytes_MCPcounter", "Monocytes_xCell")
for (f in landscape_features) a[[f]] <- land[[f]]
a$stroma_score <- rowMeans(scale(cbind(a$CAFs, a$Endothelial)))

lab_order <- unique(a$cluster_label[order(a$cluster)])
a$cluster_name <- factor(sub("^C[0-9]+: ", "", a$cluster_label), levels = sub("^C[0-9]+: ", "", lab_order))
a$immune_group <- factor(a$immune_group, levels = c("Immune-low", "Immune-high"))
a$pam50 <- factor(a$pam50, levels = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like"))
a$stage <- factor(a$stage, levels = c("I", "II", "III", "IV"))
a$histology_group <- factor(a$histology_group, levels = c("Ductal", "Lobular", "Other/mixed"))
a$receptor_subtype <- factor(a$receptor_subtype)

sig_names <- intersect(c("Cytolytic activity", "IFN-gamma response", "Antigen presentation",
                         "TGF-beta signalling", "Proliferation (G2M)"), names(tme))
cp_names  <- intersect(c("PD-L1 (CD274)", "PD-1 (PDCD1)", "PD-L2 (PDCD1LG2)", "CTLA4", "LAG3",
                         "TIM-3 (HAVCR2)", "TIGIT", "IDO1"), names(tme))
cib_names <- intersect(c("CD8 T cells", "CD4 T cells", "Tregs", "Gamma-delta T cells", "NK cells",
                         "B cells", "Plasma cells",
                         "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells",
                         "Mast cells", "Eosinophils", "Neutrophils"), names(tme))
for (f in c(sig_names, cp_names, "M2_share_of_macrophages")) a[[f]] <- tme[[f]]


cl_cols    <- c("Immune-enriched" = "#E69F00", "Stroma-rich" = "#009E73", "Depleted" = "#7F7F7F")
pam_cols   <- c("Luminal A" = "#1F4E9C", "Luminal B" = "#7FB2E5", "HER2-enriched" = "#E78AC3",
                "Basal-like" = "#D7301F", "Normal-like" = "#B8E186")
grp_cols   <- c("Immune-low" = "#D9D9D9", "Immune-high" = "#252525")
stage_cols <- c(I = "#F2F0F7", II = "#CBC9E2", III = "#9E9AC8", IV = "#54278F", Unknown = "white")
div_cols   <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
cib_cols <- c("CD8 T cells" = "#08306B", "CD4 T cells" = "#2171B5", "Tregs" = "#6BAED6",
              "Gamma-delta T cells" = "#C6DBEF", "NK cells" = "#6A51A3",
              "B cells" = "#006D2C", "Plasma cells" = "#74C476",
              "Monocytes" = "#FDD0A2", "Macrophages M0" = "#FD8D3C", "Macrophages M1" = "#D94801",
              "Macrophages M2" = "#7F2704", "Dendritic cells" = "#FEE391",
              "Mast cells" = "#525252", "Eosinophils" = "#969696", "Neutrophils" = "#D9D9D9")

theme_pub <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10), strip.background = element_blank(),
        strip.text = element_text(face = "bold"), legend.key.size = unit(0.35, "cm"),
        plot.tag = element_text(face = "bold", size = 12))

# Legend inside the panel (ggplot2 >= 3.5 uses legend.position.inside)
legend_inside <- function(x, y) {
  if (packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = c(x, y))
  else theme(legend.position = c(x, y), legend.justification = c(x, y))
}

clean_lab <- function(x) {
  map <- c(CD8_T_cells = "CD8 T cells", B_cells = "B cells", Macrophages = "Macrophages", CAFs = "CAFs",
           Endothelial = "Endothelial cells", immune_score = "Immune score", stroma_score = "Stromal score",
           M2_share_of_macrophages = "M2 share of macrophages (%)", age = "Age (years)",
           leukocyte_fraction_methylation = "Leukocyte fraction (methylation)",
           absolute_purity = "Tumour purity (ABSOLUTE)",
           ImmuneScore_xCell = "Immune score (xCell)", StromaScore_xCell = "Stroma score (xCell)",
           Tregs_quanTIseq = "Tregs (quanTIseq)", Tregs_xCell = "Tregs (xCell)",
           MyeloidDC_MCPcounter = "Myeloid DCs (MCP-counter)", MyeloidDC_xCell = "Myeloid DCs (xCell)",
           Monocytes_MCPcounter = "Monocytes (MCP-counter)", Monocytes_xCell = "Monocytes (xCell)")
  ifelse(x %in% names(map), map[x], x)
}
int <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))   # rank-based inverse normal
std <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
n_lab <- function(f) paste0(levels(f), "\n(n = ", as.integer(table(f)), ")")

save_gg <- function(p, dir, name, w, h) {
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}
save_ph <- function(dir, name, w, h, ...) {
  for (ext in c(".png", ".pdf"))
    pheatmap::pheatmap(..., filename = file.path(dir, paste0(name, ext)), width = w, height = h)
}
plan <- list()
add_plan <- function(placement, figure, file, panels, source)
  plan[[length(plan) + 1]] <<- data.frame(placement, figure, file, panels, source_script = source)

#Survival Helpers
surv_data <- function(e, var) {
  t <- a[[paste0(e, ".time")]]; s <- a[[e]]; g <- a[[var]]
  ok <- !is.na(t) & t > 0 & !is.na(s) & !is.na(g)
  data.frame(time = t[ok] / 365.25, status = s[ok], grp = droplevels(factor(g[ok], levels = levels(g))))
}

# Kaplan-Meier with numbers at risk; dashed where fewer than `low_n` remain at risk
km_plot <- function(dd, cols, title, note = NULL, ymin = 0.4, tmax = 10, low_n = 20) {
  fit <- survfit(Surv(time, status) ~ grp, data = dd)
  lv <- levels(dd$grp); cols <- cols[lv]
  st <- factor(rep(sub("^grp=", "", names(fit$strata)), fit$strata), levels = lv)
  df <- data.frame(grp = st, time = fit$time, surv = fit$surv, n.risk = fit$n.risk, n.censor = fit$n.censor)
  df <- rbind(data.frame(grp = factor(lv, levels = lv), time = 0, surv = 1,
                         n.risk = as.integer(table(dd$grp)), n.censor = 0), df)
  df <- df[order(df$grp, df$time), ]
  df <- df[df$time <= tmax, ]
  solid <- do.call(rbind, lapply(split(df, df$grp), function(x) {
    cut <- which(x$n.risk < low_n)[1]
    if (is.na(cut)) return(x)
    if (cut == 1) return(x[0, ])
    rbind(x[seq_len(cut - 1), ], transform(x[cut - 1, ], time = x$time[cut]))
  }))
  lr <- survdiff(Surv(time, status) ~ grp, data = dd)
  p_lr <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  lab <- sprintf("%s (n = %d, events = %d)", lv, as.integer(table(dd$grp)),
                 as.integer(tapply(dd$status, dd$grp, sum)))
  txt <- paste0("Log-rank p = ", format.pval(p_lr, digits = 2, eps = 1e-4), if (!is.null(note)) paste0("\n", note))
  
  p <- ggplot(df, aes(time, surv, colour = grp)) +
    geom_step(linetype = "22", linewidth = 0.5) +
    geom_step(data = solid, linewidth = 0.8) +
    geom_point(data = df[df$n.censor > 0, ], shape = 3, size = 0.9, stroke = 0.4) +
    annotate("text", x = tmax, y = 1, label = txt, hjust = 1, vjust = 1, size = 2.8) +
    scale_colour_manual(values = cols, labels = lab, name = NULL) +
    scale_x_continuous(breaks = seq(0, tmax, 2), limits = c(0, tmax)) +
    scale_y_continuous(breaks = seq(ymin, 1, 0.1)) +
    coord_cartesian(ylim = c(ymin, 1)) +
    labs(x = NULL, y = "Event-free probability", title = title) +
    theme_pub + legend_inside(0.02, 0.02) +
    theme(legend.background = element_blank(), legend.text = element_text(size = 7))
  
  rs <- summary(fit, times = seq(0, tmax, 2), extend = TRUE)
  rt <- data.frame(grp = factor(sub("^grp=", "", as.character(rs$strata)), levels = rev(lv)),
                   time = rs$time, n = rs$n.risk)
  r <- ggplot(rt, aes(time, grp, label = n, colour = grp)) + geom_text(size = 2.6) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_continuous(breaks = seq(0, tmax, 2), limits = c(0, tmax)) +
    labs(x = "Years", y = NULL, subtitle = "Number at risk") +
    theme_pub + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      plot.subtitle = element_text(size = 7.5, face = "bold"),
                      axis.text.y = element_text(size = 7))
  p / r + plot_layout(heights = c(3.4, 1))
}

# Adjusted spline of the immune score (same model as survival_analysis.R)
spline_plot <- function(e) {
  d <- a
  d$immune_score_sd <- as.numeric(scale(d$immune_score)); d$age10 <- d$age / 10
  d$pam50r <- relevel(droplevels(d$pam50), ref = "Luminal A")
  d$time <- d[[paste0(e, ".time")]] / 365.25; d$status <- d[[e]]
  d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$status) &
           complete.cases(d[, c("age10", "stage", "pam50r", "immune_score_sd")]), ]
  d <- droplevels(d)
  basis <- ns(d$immune_score_sd, df = 3); d$imm_ns <- basis
  f_ns  <- coxph(Surv(time, status) ~ imm_ns + age10 + stage + pam50r, data = d)
  f_lin <- coxph(Surv(time, status) ~ immune_score_sd + age10 + stage + pam50r, data = d)
  p_nl <- pchisq(2 * (f_ns$loglik[2] - f_lin$loglik[2]), 2, lower.tail = FALSE)
  grid <- seq(quantile(d$immune_score_sd, 0.02), quantile(d$immune_score_sd, 0.98), length.out = 200)
  B <- predict(basis, grid) - matrix(predict(basis, 0), length(grid), 3, byrow = TRUE)
  k <- startsWith(names(coef(f_ns)), "imm_ns")
  lp <- drop(B %*% coef(f_ns)[k]); se <- sqrt(rowSums((B %*% vcov(f_ns)[k, k]) * B))
  cv <- data.frame(x = grid, hr = exp(lp), lo = exp(lp - 1.96 * se), hi = exp(lp + 1.96 * se),
                   lin = exp(coef(f_lin)["immune_score_sd"] * grid))
  hr_lin <- summary(f_lin)$conf.int["immune_score_sd", ]
  ggplot(cv, aes(x)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#E69F00", alpha = 0.25) +
    geom_line(aes(y = hr), linewidth = 0.8) +
    geom_line(aes(y = lin), linetype = "dotted", colour = "grey30") +
    geom_rug(data = d[d$immune_score_sd >= min(grid) & d$immune_score_sd <= max(grid), ],
             aes(x = immune_score_sd), inherit.aes = FALSE, alpha = 0.3, length = unit(0.02, "npc")) +
    scale_y_log10() +
    annotate("text", x = max(grid), y = max(cv$hi), hjust = 1, vjust = 1, size = 2.7,
             label = sprintf("Linear HR per SD %.2f (%.2f-%.2f)\nNon-linearity p = %s",
                             hr_lin[1], hr_lin[3], hr_lin[4], format.pval(p_nl, digits = 2))) +
    labs(x = "Immune score (SD from cohort mean)", y = "Adjusted hazard ratio (log scale)",
         title = sprintf("%s: immune score (spline, 3 df)", e),
         caption = "Adjusted for age, stage and PAM50; reference = cohort mean; dotted = linear fit") +
    theme_pub
}

#Main figure 1: immune landscape heatmap
o <- order(a$cluster_name, -a$immune_score)
na_unknown <- function(x, lv) factor(ifelse(is.na(x), "Unknown", as.character(x)), levels = c(lv, "Unknown"))
ann <- data.frame(Cluster = a$cluster_name[o], `Immune group` = a$immune_group[o],
                  PAM50 = droplevels(na_unknown(a$pam50[o], levels(a$pam50))),
                  Stage = droplevels(na_unknown(a$stage[o], levels(a$stage))),
                  row.names = a$sample_id[o], check.names = FALSE)
ann_cols <- list(Cluster = cl_cols, `Immune group` = grp_cols,
                 PAM50 = c(pam_cols, Unknown = "white")[levels(ann$PAM50)],
                 Stage = stage_cols[levels(ann$Stage)])
make_hm <- function(feats) {
  hm <- t(sapply(feats, function(f) int(a[[f]])[o])); hm[hm > 2.5] <- 2.5; hm[hm < -2.5] <- -2.5
  rownames(hm) <- clean_lab(feats); colnames(hm) <- a$sample_id[o]; hm
}
gaps_c <- cumsum(table(a$cluster_name))[-nlevels(a$cluster_name)]
save_ph(main_dir, "immune_landscape_heatmap", 11, 4.6,
        mat = make_hm(landscape_features[1:7]), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
        annotation_col = ann, annotation_colors = ann_cols, gaps_col = gaps_c, gaps_row = 5,
        breaks = seq(-2.5, 2.5, length.out = 101), color = div_cols, fontsize = 8, border_color = NA,
        main = sprintf("Immune landscape of TCGA-BRCA (n = %d); values: rank-normalised scores", nrow(a)))
add_plan("Main", "Figure 1", "immune_landscape_heatmap",
         "Heatmap of the five consensus features and xCell overall scores; tumours ordered by cluster then immune score; annotations: cluster, immune group, PAM50, stage",
         "phase4_immune_landscape.R")

save_ph(supp_dir, "immune_landscape_heatmap_all_features", 11, 6,
        mat = make_hm(landscape_features), cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
        annotation_col = ann, annotation_colors = ann_cols, gaps_col = gaps_c, gaps_row = c(5, 7),
        breaks = seq(-2.5, 2.5, length.out = 101), color = div_cols, fontsize = 8, border_color = NA,
        main = "Immune landscape including exploratory features (Tregs, myeloid DCs, monocytes)")
add_plan("Supplementary", "", "immune_landscape_heatmap_all_features",
         "Figure 1 with exploratory single-method features added", "phase4_immune_landscape.R")

#Main figure 2: clusters and clinical context
# A. cluster profiles
prof <- do.call(rbind, lapply(landscape_features[1:5], function(f) {
  z <- int(a[[f]])
  data.frame(feature = clean_lab(f), cluster = levels(a$cluster_name), z = as.numeric(tapply(z, a$cluster_name, mean)))
}))
prof$feature <- factor(prof$feature, levels = rev(clean_lab(landscape_features[1:5])))
prof$cluster <- factor(prof$cluster, levels = levels(a$cluster_name), labels = n_lab(a$cluster_name))
p2a <- ggplot(prof, aes(cluster, feature, fill = z)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", z)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1.3, 1.3), name = "Mean z") +
  labs(x = NULL, y = NULL, title = "Cluster profiles") + theme_pub + theme(axis.line = element_blank())

# B. RNA immune score vs DNA methylation leukocyte fraction
ok <- !is.na(a$leukocyte_fraction_methylation)
rho <- cor(a$immune_score[ok], a$leukocyte_fraction_methylation[ok], method = "spearman")
p2b <- ggplot(a[ok, ], aes(immune_score, leukocyte_fraction_methylation, colour = cluster_name)) +
  geom_point(size = 0.7, alpha = 0.7) + scale_colour_manual(values = cl_cols, name = NULL) +
  labs(x = "Immune score (RNA consensus)", y = "Leukocyte fraction (DNA methylation)",
       title = "Independent validation by DNA methylation",
       subtitle = sprintf("Spearman rho = %.2f, n = %d", rho, sum(ok))) +
  theme_pub + legend_inside(0.02, 0.98)

# C. cluster composition within clinical groups, with n per bar
cv <- c("PAM50" = "pam50", "Receptor subtype" = "receptor_subtype", "Stage" = "stage", "Histology" = "histology_group")
comp_cl <- do.call(rbind, lapply(names(cv), function(lab) {
  v <- a[[cv[[lab]]]]; ok <- !is.na(v)
  tb <- table(droplevels(factor(v[ok])), a$cluster_name[ok])
  data.frame(variable = lab, category = rep(rownames(tb), ncol(tb)), cluster = rep(colnames(tb), each = nrow(tb)),
             count = as.vector(tb), n = rep(rowSums(tb), ncol(tb)))
}))
comp_cl$prop <- comp_cl$count / comp_cl$n
comp_cl$variable <- factor(comp_cl$variable, levels = names(cv))
comp_cl$category <- factor(comp_cl$category, levels = unique(c(levels(a$pam50), levels(a$receptor_subtype),
                                                               levels(a$stage), levels(a$histology_group))))
comp_cl$cluster <- factor(comp_cl$cluster, levels = rev(levels(a$cluster_name)))
n_bar <- unique(comp_cl[, c("variable", "category", "n")])
p2c <- ggplot(comp_cl, aes(category, prop, fill = cluster)) + geom_col(width = 0.8) +
  geom_text(data = n_bar, aes(category, 1.04, label = paste0("n=", n)), inherit.aes = FALSE, size = 2.3) +
  facet_grid(. ~ variable, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cl_cols, breaks = levels(a$cluster_name), name = "Immune cluster") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), breaks = seq(0, 1, 0.25), limits = c(0, 1.08)) +
  labs(x = NULL, y = "Tumours", title = "Immune cluster composition within clinical groups") +
  theme_pub + theme(axis.text.x = element_text(angle = 40, hjust = 1))

fig2 <- (p2a | p2b) / p2c + plot_layout(heights = c(1, 1.15)) + plot_annotation(tag_levels = "A")
save_gg(fig2, main_dir, "immune_clusters_and_clinical_context", 10, 8.5)
add_plan("Main", "Figure 2", "immune_clusters_and_clinical_context",
         "A cluster profiles (mean z of consensus features); B RNA immune score vs methylation leukocyte fraction; C cluster composition within PAM50, receptor subtype, stage and histology (n per bar)",
         "phase4_immune_landscape.R; clinical_immune_integration.R")

#Main figure 3: tumour microenvironment by cluster
# A. CIBERSORT composition coloured by lineage
cib_long <- do.call(rbind, lapply(levels(a$cluster_name), function(cl)
  data.frame(cluster = cl, cell_type = cib_names, pct = colMeans(tme[a$cluster_name == cl, cib_names, drop = FALSE], na.rm = TRUE))))
cib_long$cluster <- factor(cib_long$cluster, levels = levels(a$cluster_name), labels = n_lab(a$cluster_name))
cib_long$cell_type <- factor(cib_long$cell_type, levels = rev(cib_names))
p3a <- ggplot(cib_long, aes(cluster, pct, fill = cell_type)) + geom_col(width = 0.7, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = cib_cols, breaks = cib_names, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Mean share of immune infiltrate (%)", title = "Immune composition (CIBERSORT, relative)",
       subtitle = "Blue: T/NK; green: B lineage; orange: myeloid; grey: granulocytes") +
  theme_pub

# B. cluster mean z of immune features, functional signatures and checkpoints
zf <- c(setNames(landscape_features[1:5], rep("Immune and stromal cells", 5)),
        setNames(sig_names, rep("Functional signatures", length(sig_names))),
        setNames(cp_names, rep("Immune checkpoints", length(cp_names))))
tp <- do.call(rbind, lapply(seq_along(zf), function(i) {
  z <- std(a[[zf[i]]])
  data.frame(block = names(zf)[i], feature = clean_lab(zf[i]), cluster = levels(a$cluster_name),
             z = as.numeric(tapply(z, a$cluster_name, mean, na.rm = TRUE)))
}))
tp$block <- factor(tp$block, levels = unique(names(zf)))
tp$feature <- factor(tp$feature, levels = rev(clean_lab(zf)))
tp$cluster <- factor(tp$cluster, levels = levels(a$cluster_name))
p3b <- ggplot(tp, aes(cluster, feature, fill = z)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", z)), size = 2.4) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1.3, 1.3), name = "Mean z") +
  labs(x = NULL, y = NULL, title = "TME profile by cluster") +
  theme_pub + theme(axis.line = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
                    strip.text.y = element_text(angle = 0, hjust = 0))

fig3 <- (p3a | p3b) + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "A")
save_gg(fig3, main_dir, "tumour_microenvironment_by_cluster", 11, 6.5)
add_plan("Main", "Figure 3", "tumour_microenvironment_by_cluster",
         "A CIBERSORT immune composition by cluster coloured by lineage; B cluster mean z-scores of cell features, functional signatures and checkpoint genes",
         "tme_characterization.R")

#Main figure 4: Survival
hr_ie <- mv[mv$endpoint == "PFI" & mv$exposure == "Immune cluster (ref Depleted)" & mv$level == "Immune-enriched", ]
note_ie <- if (nrow(hr_ie)) sprintf("Adjusted HR Immune-enriched vs Depleted\n%.2f (%.2f-%.2f), p = %s",
                                    hr_ie$HR, hr_ie$CI_lower, hr_ie$CI_upper, format.pval(hr_ie$p, digits = 2)) else NULL
p4a <- km_plot(surv_data("PFI", "cluster_name"), cl_cols, "Progression-free interval by immune cluster", note_ie)
p4b <- spline_plot("PFI")

fp <- mv[mv$endpoint %in% c("PFI", "OS"), ]
fp$label <- ifelse(fp$exposure == "Immune cluster (ref Depleted)", paste(fp$level, "vs Depleted"),
                   sub(" \\(per SD\\)", "", fp$exposure))
fp$label <- factor(fp$label, levels = rev(unique(fp$label)))
fp$endpoint <- factor(fp$endpoint, levels = c("PFI", "OS"))
fp$sig <- ifelse(!is.na(fp$FDR_immune_features) & fp$FDR_immune_features < 0.05, "*", "")
p4c <- ggplot(fp, aes(HR, label, colour = endpoint, shape = endpoint)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), width = 0, orientation = "y", position = position_dodge(0.55)) +
  geom_point(size = 1.8, position = position_dodge(0.55), fill = "white") +
  geom_text(aes(x = CI_upper, label = sig), position = position_dodge(0.55), hjust = -0.4, size = 3.5,
            show.legend = FALSE) +
  scale_x_log10() + scale_colour_manual(values = c(PFI = "black", OS = "grey55"), name = NULL) +
  scale_shape_manual(values = c(PFI = 16, OS = 21), name = NULL) +
  labs(x = "Adjusted hazard ratio (95% CI, log scale)", y = NULL,
       title = "Immune features and outcome",
       caption = "Continuous features per SD; adjusted for age, stage, PAM50; * FDR < 0.05") +
  theme_pub + theme(legend.position = "top")

fig4 <- wrap_elements(p4a) | (p4b / p4c + plot_layout(heights = c(1, 1.3)))
fig4 <- fig4 + plot_layout(widths = c(1.1, 1)) + plot_annotation(tag_levels = "A")
save_gg(fig4, main_dir, "survival_by_immune_infiltration", 11.5, 7.5)
add_plan("Main", "Figure 4", "survival_by_immune_infiltration",
         "A Kaplan-Meier PFI by cluster (dashed where <20 at risk; y-axis from 0.4); B adjusted spline of immune score for PFI; C adjusted hazard ratios for PFI and OS",
         "survival_analysis.R")

#Supplementary figures (redrawn)
# Kaplan-Meier, all four endpoints
ep <- c(PFI = "Progression-free interval", OS = "Overall survival", DFI = "Disease-free interval",
        DSS = "Disease-specific survival")
km_grid <- function(var, cols, name, title) {
  pl <- lapply(names(ep), function(e) wrap_elements(km_plot(surv_data(e, var), cols, ep[[e]])))
  g <- wrap_plots(pl, ncol = 2) + plot_annotation(title = title,
                                                  caption = "Dashed: fewer than 20 patients at risk. Y-axis starts at 0.4.",
                                                  theme = theme(plot.title = element_text(face = "bold")))
  save_gg(g, supp_dir, name, 11, 10)
  add_plan("Supplementary", "", name, paste("Kaplan-Meier, four endpoints:", title), "survival_analysis.R")
}
km_grid("cluster_name", cl_cols, "km_by_immune_cluster_all_endpoints", "By immune cluster")
km_grid("immune_group", c("Immune-low" = "#9E9E9E", "Immune-high" = "#252525"),
        "km_by_immune_group_all_endpoints", "By immune group (mixture-model split)")
a$basal_hot_cold <- factor(ifelse(a$pam50 %in% "Basal-like" & a$cluster_name %in% c("Immune-enriched", "Depleted"),
                                  ifelse(a$cluster_name == "Immune-enriched", "Hot (Immune-enriched)", "Cold (Depleted)"), NA),
                           levels = c("Cold (Depleted)", "Hot (Immune-enriched)"))
km_grid("basal_hot_cold", c("Cold (Depleted)" = "#7F7F7F", "Hot (Immune-enriched)" = "#E69F00"),
        "km_basal_like_hot_vs_cold_all_endpoints", "Basal-like tumours: hot vs cold")

save_gg(spline_plot("OS"), supp_dir, "immune_score_spline_hazard_ratio_OS", 6, 4.5)
add_plan("Supplementary", "", "immune_score_spline_hazard_ratio_OS",
         "Adjusted spline of immune score for OS (lower tail rests on few patients)", "survival_analysis.R")

# DNA-based measures by cluster
vlong <- rbind(data.frame(measure = "Leukocyte fraction (DNA methylation)", cluster = a$cluster_name, value = a$leukocyte_fraction_methylation),
               data.frame(measure = "Tumour purity (ABSOLUTE)", cluster = a$cluster_name, value = a$absolute_purity))
vlong <- vlong[!is.na(vlong$value), ]
ps <- ggplot(vlong, aes(cluster, value, fill = cluster)) + geom_boxplot(outlier.size = 0.4, width = 0.6) +
  facet_wrap(~ measure, scales = "free_y") + scale_fill_manual(values = cl_cols, guide = "none") +
  labs(x = NULL, y = NULL, title = "DNA-based measures by immune cluster") + theme_pub
save_gg(ps, supp_dir, "methylation_and_purity_by_cluster", 7, 3.8)
add_plan("Supplementary", "", "methylation_and_purity_by_cluster",
         "Methylation leukocyte fraction and ABSOLUTE purity by cluster", "clinical_immune_integration.R")

# Functional signals by cluster
blong <- do.call(rbind, lapply(c(sig_names, cp_names), function(f)
  data.frame(feature = f, cluster = a$cluster_name, value = a[[f]])))
blong$feature <- factor(blong$feature, levels = c(sig_names, cp_names))
ps <- ggplot(blong, aes(cluster, value, fill = cluster)) + geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~ feature, scales = "free_y", ncol = 5) + scale_fill_manual(values = cl_cols, guide = "none") +
  labs(x = NULL, y = "Score (signatures: mean z; checkpoints and cytolytic activity: log2 expression)",
       title = "Functional signals and checkpoint genes by immune cluster") +
  theme_pub + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_gg(ps, supp_dir, "functional_signals_by_cluster", 11, 6.5)
add_plan("Supplementary", "", "functional_signals_by_cluster",
         "Distributions of signatures and checkpoint genes by cluster", "tme_characterization.R")

# TME correlation heatmap with readable labels
cm <- as.matrix(read.csv(file.path(res_dir, "tme_characterization", "tme_feature_spearman_correlations.csv"),
                         row.names = 1, check.names = FALSE))
rownames(cm) <- colnames(cm) <- clean_lab(colnames(cm))
save_ph(supp_dir, "tme_feature_correlation_heatmap", 9, 8.5, mat = cm, display_numbers = TRUE,
        number_format = "%.2f", fontsize_number = 5.5, fontsize = 8, breaks = seq(-1, 1, length.out = 101),
        color = div_cols, border_color = NA, main = "Spearman correlations between TME features")
add_plan("Supplementary", "", "tme_feature_correlation_heatmap",
         "Spearman correlations between cell features, signatures, checkpoints, methylation LF and purity", "tme_characterization.R")

# Basal-like hot vs cold, TME differences
hc <- read.csv(file.path(res_dir, "tme_characterization", "basal_like_hot_vs_cold.csv"), check.names = FALSE)
hc$feature <- clean_lab(hc$feature)
hc$feature <- factor(hc$feature, levels = hc$feature[order(hc$difference)])
hc$sig <- ifelse(hc$q_value_fdr < 0.05, "FDR < 0.05", "Not significant")
ps <- ggplot(hc, aes(difference, feature, colour = sig)) + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(size = 2.2) + scale_colour_manual(values = c("FDR < 0.05" = "#E69F00", "Not significant" = "grey60"), name = NULL) +
  labs(x = "Median difference, hot minus cold (z units)", y = NULL,
       title = sprintf("Basal-like: Immune-enriched (n = %d) vs Depleted (n = %d)", hc$n_hot[1], hc$n_cold[1])) +
  theme_pub
save_gg(ps, supp_dir, "basal_like_hot_vs_cold_tme", 7, 5.5)
add_plan("Supplementary", "", "basal_like_hot_vs_cold_tme",
         "TME differences between hot and cold basal-like tumours", "tme_characterization.R")

# Integrated sample-level TME heatmap
zfe <- c(landscape_features[1:5], sig_names, cp_names)
Z <- t(sapply(zfe, function(f) std(a[[f]])[o])); Z[Z > 2.5] <- 2.5; Z[Z < -2.5] <- -2.5
rownames(Z) <- clean_lab(zfe); colnames(Z) <- a$sample_id[o]
ann2 <- data.frame(Cluster = a$cluster_name[o], PAM50 = ann$PAM50, Stage = ann$Stage,
                   `Leukocyte fraction` = a$leukocyte_fraction_methylation[o], Purity = a$absolute_purity[o],
                   row.names = a$sample_id[o], check.names = FALSE)
ann2_cols <- c(ann_cols[c("Cluster", "PAM50", "Stage")],
               list(`Leukocyte fraction` = colorRampPalette(c("white", "#7F2704"))(50),
                    Purity = colorRampPalette(c("white", "black"))(50)))
save_ph(supp_dir, "integrated_tme_heatmap", 11, 6.5, mat = Z, cluster_rows = FALSE, cluster_cols = FALSE,
        show_colnames = FALSE, annotation_col = ann2, annotation_colors = ann2_cols, gaps_col = gaps_c,
        gaps_row = c(5, 5 + length(sig_names)), breaks = seq(-2.5, 2.5, length.out = 101), color = div_cols,
        fontsize = 8, border_color = NA, main = "Integrated tumour microenvironment (z-scores)")
add_plan("Supplementary", "", "integrated_tme_heatmap",
         "Sample-level heatmap of cell features, signatures and checkpoints with methylation LF and purity", "tme_characterization.R")

#Supplementary figures (redrawn)
copies <- list(
  c("", "immune_landscape/feature_distribution_violin_plot.png", "feature_distributions", "Distributions of all deconvolution features"),
  c("", "immune_landscape/cluster_number_selection_plot.png", "cluster_number_selection", "Choice of k for k-means"),
  c("", "immune_landscape/feature_spearman_correlation_heatmap.png", "landscape_feature_correlations", "Correlations between landscape features"),
  c("", "clinical_immune_integration/immune_and_stromal_scores_by_subtype_and_stage.png", "immune_and_stromal_scores_by_subtype_and_stage", "Scores by PAM50, receptor subtype, stage"),
  c("", "survival_analysis/immune_score_by_subtype_forest_plot.png", "immune_score_effect_by_subtype", "Immune score HR within subtypes"),
  c("", "survival_analysis/schoenfeld_residuals_immune_score_PFI.png", "schoenfeld_residuals_immune_score_PFI", "Proportional hazards check for the immune score"),
  c("", "clinical_integration/receptor_status_vs_expression.png", "receptor_status_vs_expression", "QC: clinical ER/PR/HER2 vs mRNA"),
  c("", "clinical_integration/sex_check_expression.png", "sex_check_expression", "QC: sex vs Y-chromosome genes"))
for (cp in copies) {
  src <- file.path(res_dir, cp[2])
  if (file.exists(src)) {
    file.copy(src, file.path(supp_dir, paste0(cp[3], ".png")), overwrite = TRUE)
    add_plan("Supplementary", "", cp[3], cp[4], paste("copied from", dirname(cp[2])))
  } else message("Not found, skipped: ", cp[2])
}

plan_df <- do.call(rbind, plan)
supp_order <- c("feature_distributions", "landscape_feature_correlations", "cluster_number_selection",
                "immune_landscape_heatmap_all_features", "methylation_and_purity_by_cluster",
                "immune_and_stromal_scores_by_subtype_and_stage", "integrated_tme_heatmap",
                "tme_feature_correlation_heatmap", "functional_signals_by_cluster", "basal_like_hot_vs_cold_tme",
                "km_by_immune_cluster_all_endpoints", "km_by_immune_group_all_endpoints",
                "km_basal_like_hot_vs_cold_all_endpoints", "immune_score_spline_hazard_ratio_OS",
                "immune_score_effect_by_subtype", "schoenfeld_residuals_immune_score_PFI",
                "receptor_status_vs_expression", "sex_check_expression")
s_idx <- match(plan_df$file, supp_order)
plan_df$figure[plan_df$placement == "Supplementary"] <- paste0("Figure S", rank(s_idx[plan_df$placement == "Supplementary"]))
plan_df <- plan_df[order(plan_df$placement != "Main", ifelse(is.na(s_idx), 0, s_idx), plan_df$figure), ]
write.csv(plan_df, file.path(res_dir, "figures", "figure_plan.csv"), row.names = FALSE)

print(plan_df[, c("placement", "figure", "file")], row.names = FALSE)
cat("\nMain:", main_dir, "\nSupplementary:", supp_dir, "\n")
