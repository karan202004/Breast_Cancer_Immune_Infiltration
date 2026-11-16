#Survival analysis 


suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
meta_dir <- file.path(project_dir, "metadata")
out_dir  <- file.path(project_dir, "results", "survival_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


m    <- readRDS(file.path(meta_dir, "clinical_master_table.rds"))
land <- readRDS(file.path(project_dir, "results", "immune_landscape", "sample_level_results.rds"))
a <- m[m$in_analysis_set, ]
if (!setequal(a$sample_id, land$sample_id))
  stop("Master table and immune landscape results differ. Run scripts 1-3 first.")
li <- match(a$sample_id, land$sample_id)
features <- c("CD8_T_cells", "B_cells", "Macrophages", "CAFs", "Endothelial")
for (f in features) a[[f]] <- land[[f]][li]
a$stroma_score <- rowMeans(scale(cbind(a$CAFs, a$Endothelial)))   # as in clinical_immune_integration.R

# Exposures (continuous ones per 1 SD)
cont_exposures <- c(immune_score = "Immune score", stroma_score = "Stromal score",
                    CD8_T_cells = "CD8 T cells", B_cells = "B cells", Macrophages = "Macrophages",
                    CAFs = "CAFs", Endothelial = "Endothelial cells")
for (v in names(cont_exposures)) a[[paste0(v, "_sd")]] <- as.numeric(scale(a[[v]]))

lab_order <- unique(a$cluster_label[order(a$cluster)])
a$cluster_name <- factor(sub("^C[0-9]+: ", "", a$cluster_label), levels = sub("^C[0-9]+: ", "", lab_order))
a$cluster_ref  <- relevel(a$cluster_name, ref = "Depleted")
a$immune_group <- factor(a$immune_group, levels = c("Immune-low", "Immune-high"))

# Covariates
a$age10 <- a$age / 10
a$stage <- factor(a$stage, levels = c("I", "II", "III", "IV"))
a$stage_2 <- factor(ifelse(a$stage %in% c("I", "II"), "I-II", ifelse(is.na(a$stage), NA, "III-IV")),
                    levels = c("I-II", "III-IV"))
a$pam50 <- relevel(factor(a$pam50), ref = "Luminal A")
a$histology_group <- factor(a$histology_group, levels = c("Ductal", "Lobular", "Other/mixed"))
a$menopause_status <- factor(a$menopause_status, levels = c("Pre/peri-menopausal", "Post-menopausal"))
a$receptor_subtype <- factor(a$receptor_subtype)
if ("HR+/HER2-" %in% levels(a$receptor_subtype)) a$receptor_subtype <- relevel(a$receptor_subtype, ref = "HR+/HER2-")
# Three subtype groups following Denkert et al. 2018 (luminal / HER2 / triple-negative-like)
a$subtype3 <- factor(ifelse(a$pam50 %in% c("Luminal A", "Luminal B", "Normal-like"), "Luminal",
                            ifelse(is.na(a$pam50), NA, as.character(a$pam50))),
                     levels = c("Luminal", "HER2-enriched", "Basal-like"))
a$basal_hot_cold <- factor(ifelse(a$pam50 %in% "Basal-like" & a$cluster_name %in% c("Immune-enriched", "Depleted"),
                                  ifelse(a$cluster_name == "Immune-enriched", "Hot (Immune-enriched)", "Cold (Depleted)"), NA),
                           levels = c("Cold (Depleted)", "Hot (Immune-enriched)"))

# Endpoints: time in years; time <= 0 set to missing
endpoints <- c(PFI = "Progression-free interval", OS = "Overall survival",
               DFI = "Disease-free interval", DSS = "Disease-specific survival")
for (e in names(endpoints)) {
  t <- a[[paste0(e, ".time")]]
  bad <- is.na(t) | t <= 0 | is.na(a[[e]])
  a[[paste0(e, "_years")]]  <- ifelse(bad, NA, t / 365.25)
  a[[paste0(e, "_status")]] <- ifelse(bad, NA, a[[e]])
}

# Endpoint-specific data set with generic time/status columns
ep_data <- function(e, d = a, censor_at = NULL) {
  d$time <- d[[paste0(e, "_years")]]; d$status <- d[[paste0(e, "_status")]]
  d <- d[!is.na(d$time) & !is.na(d$status), ]
  if (!is.null(censor_at)) {
    d$status <- ifelse(d$time > censor_at, 0, d$status)
    d$time <- pmin(d$time, censor_at)
  }
  d
}
covars_main <- c("age10", "stage", "pam50")
cc_data <- function(d, vars) droplevels(d[complete.cases(d[, vars, drop = FALSE]), ])

# Helpers
cox_fit <- function(rhs, d) coxph(as.formula(paste("Surv(time, status) ~", rhs)), data = d, ties = "efron")

# Hazard-ratio rows for terms whose name starts with `var`
hr_rows <- function(fit, var) {
  s <- summary(fit); co <- s$coefficients; ci <- s$conf.int
  k <- startsWith(rownames(co), var)
  data.frame(term = rownames(co)[k], level = sub(paste0("^", var), "", rownames(co)[k]),
             HR = ci[k, "exp(coef)"], CI_lower = ci[k, "lower .95"], CI_upper = ci[k, "upper .95"],
             p = co[k, "Pr(>|z|)"], row.names = NULL)
}
# Likelihood-ratio test of `var` in a fitted model (refit without it on the same rows)
lrt_p <- function(fit, d, rhs_without) {
  if (is.null(rhs_without)) return(pchisq(2 * diff(fit$loglik), length(coef(fit)), lower.tail = FALSE))
  red <- cox_fit(rhs_without, d)
  pchisq(2 * (fit$loglik[2] - red$loglik[2]), length(coef(fit)) - length(coef(red)), lower.tail = FALSE)
}
fmt_hr <- function(h, l, u) sprintf("%.2f (%.2f-%.2f)", h, l, u)
# Round estimates to 3 significant digits; leave counts as they are and AIC to 1 decimal
r3 <- function(df) {
  for (v in names(df)) if (is.numeric(df[[v]])) {
    if (grepl("^AIC", v)) df[[v]] <- round(df[[v]], 1)
    else if (!all(df[[v]] == round(df[[v]]), na.rm = TRUE)) df[[v]] <- signif(df[[v]], 3)
  }
  df
}

#Kaplan-Meier 
km_groups <- list(
  list(var = "cluster_name",   title = "Immune cluster",         file = "km_by_immune_cluster.png",
       cols = setNames(hcl.colors(3, "Dark 3"), levels(a$cluster_name))),
  list(var = "immune_group",   title = "Immune group",           file = "km_by_immune_group.png",
       cols = c("Immune-low" = "#2166AC", "Immune-high" = "#B2182B")),
  list(var = "basal_hot_cold", title = "Basal-like: hot vs cold", file = "km_basal_like_hot_vs_cold.png",
       cols = c("Cold (Depleted)" = "#2166AC", "Hot (Immune-enriched)" = "#B2182B")))

km_rows <- list()
for (g in km_groups) {
  png(file.path(out_dir, g$file), width = 12, height = 10.5, units = "in", res = 300)
  par(mfrow = c(2, 2), oma = c(0, 0, 2.5, 0))
  for (e in names(endpoints)) {
    d <- ep_data(e); d <- d[!is.na(d[[g$var]]), ]; d$grp <- droplevels(d[[g$var]])
    fit <- survfit(Surv(time, status) ~ grp, data = d)
    lr  <- survdiff(Surv(time, status) ~ grp, data = d)
    p_lr <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
    lv <- levels(d$grp); cols <- g$cols[lv]
    
    par(mar = c(5 + length(lv) + 1, 10, 3, 1))
    plot(fit, col = cols, lwd = 2, xlim = c(0, 10), ylim = c(0, 1), xaxt = "n", mark.time = TRUE,
         xlab = "", ylab = "Event-free probability", main = endpoints[[e]])
    axis(1, at = seq(0, 10, 2)); mtext("Years", side = 1, line = 2.3, cex = 0.8)
    legend("bottomleft", legend = sprintf("%s (n=%d, events=%d)", lv, as.integer(table(d$grp)),
                                          as.integer(tapply(d$status, d$grp, sum))),
           col = cols, lwd = 2, bty = "n", cex = 0.8)
    legend("topright", legend = sprintf("Log-rank p = %s", format.pval(p_lr, digits = 2)), bty = "n", cex = 0.85)
    at <- seq(0, 10, 2)
    s <- summary(fit, times = at, extend = TRUE)
    nr <- split(s$n.risk, sub("^grp=", "", as.character(s$strata)))
    mtext("At risk", side = 1, line = 4, at = -0.6, adj = 1, cex = 0.7, font = 2)
    for (i in seq_along(lv)) {
      mtext(lv[i], side = 1, line = 4 + i, at = -0.6, adj = 1, cex = 0.7, col = cols[i])
      mtext(nr[[lv[i]]], side = 1, line = 4 + i, at = at, cex = 0.7, col = cols[i])
    }
    s5 <- summary(fit, times = 5, extend = TRUE)
    km_rows[[length(km_rows) + 1]] <- data.frame(
      grouping = g$title, endpoint = e, group = sub("^grp=", "", as.character(s5$strata)),
      n = as.integer(table(d$grp)), events = as.integer(tapply(d$status, d$grp, sum)),
      at_risk_5y = s5$n.risk, event_free_5y_pct = round(100 * s5$surv, 1),
      CI_lower_pct = round(100 * s5$lower, 1), CI_upper_pct = round(100 * s5$upper, 1),
      logrank_p = signif(p_lr, 3))
  }
  mtext(sprintf("Kaplan-Meier by %s (TCGA-BRCA)", tolower(g$title)), outer = TRUE, cex = 1.2, font = 2)
  dev.off()
}
km_tab <- do.call(rbind, km_rows)
write.csv(km_tab, file.path(out_dir, "km_event_free_estimates.csv"), row.names = FALSE)

#Univariate COX
uni_vars <- c(setNames(paste0(names(cont_exposures), "_sd"), paste(cont_exposures, "(per SD)")),
              "Immune cluster (ref Depleted)" = "cluster_ref",
              "Immune group (ref Immune-low)" = "immune_group",
              "Age (per 10 years)" = "age10", "Stage (ref I)" = "stage", "PAM50 (ref Luminal A)" = "pam50",
              "Histology (ref Ductal)" = "histology_group",
              "Menopausal status" = "menopause_status", "Receptor subtype" = "receptor_subtype")
uni_rows <- list()
for (e in names(endpoints)) {
  d0 <- ep_data(e)
  for (lab in names(uni_vars)) {
    v <- uni_vars[[lab]]; d <- cc_data(d0, v)
    if (is.factor(d[[v]]) && nlevels(d[[v]]) < 2) next
    fit <- cox_fit(v, d)
    r <- hr_rows(fit, v)
    uni_rows[[length(uni_rows) + 1]] <- data.frame(endpoint = e, variable = lab, r[, -1],
                                                   n = fit$n, events = fit$nevent, variable_LRT_p = lrt_p(fit, d, NULL))
  }
}
uni_tab <- do.call(rbind, uni_rows)
uni_tab$level[uni_tab$level == ""] <- "per unit"
feat_rows <- uni_tab$variable %in% paste(cont_exposures, "(per SD)")
uni_tab$FDR_immune_features <- NA
for (e in names(endpoints)) {
  k <- feat_rows & uni_tab$endpoint == e
  uni_tab$FDR_immune_features[k] <- p.adjust(uni_tab$p[k], "BH")
}
uni_tab$HR_95CI <- fmt_hr(uni_tab$HR, uni_tab$CI_lower, uni_tab$CI_upper)
write.csv(r3(uni_tab), file.path(out_dir, "univariate_cox_results.csv"), row.names = FALSE)

#Mutivariable COX
mv_exposures <- c(setNames(paste0(names(cont_exposures), "_sd"), paste(cont_exposures, "(per SD)")),
                  "Immune cluster (ref Depleted)" = "cluster_ref")
mv_rows <- list(); full_rows <- list(); ph_rows <- list(); main_fits <- list()
for (e in names(endpoints)) {
  d <- cc_data(ep_data(e), covars_main)
  adj <- paste(covars_main, collapse = " + ")
  base_terms <- length(coef(cox_fit(adj, d)))
  for (lab in names(mv_exposures)) {
    v <- mv_exposures[[lab]]
    fit <- cox_fit(paste(v, "+", adj), d)
    r <- hr_rows(fit, v)
    mv_rows[[length(mv_rows) + 1]] <- data.frame(endpoint = e, exposure = lab, r[, -1],
                                                 n = fit$n, events = fit$nevent, model_terms = length(coef(fit)),
                                                 events_per_term = round(fit$nevent / length(coef(fit)), 1),
                                                 exposure_LRT_p = lrt_p(fit, d, adj))
    if (v %in% c("immune_score_sd", "cluster_ref")) {
      main_fits[[paste(e, v)]] <- fit
      s <- summary(fit)
      full_rows[[length(full_rows) + 1]] <- data.frame(endpoint = e, model = lab, term = rownames(s$coefficients),
                                                       HR = s$conf.int[, 1], CI_lower = s$conf.int[, 3], CI_upper = s$conf.int[, 4],
                                                       p = s$coefficients[, "Pr(>|z|)"], row.names = NULL)
      z <- cox.zph(fit)
      ph_rows[[length(ph_rows) + 1]] <- data.frame(endpoint = e, model = lab, term = rownames(z$table),
                                                   chisq = z$table[, "chisq"], df = z$table[, "df"], p = z$table[, "p"], row.names = NULL)
    }
  }
}
mv_tab <- do.call(rbind, mv_rows)
mv_tab$level[mv_tab$level == ""] <- "per SD"
mv_tab$FDR_immune_features <- NA
for (e in names(endpoints)) {
  k <- mv_tab$endpoint == e & mv_tab$exposure %in% paste(cont_exposures, "(per SD)")
  mv_tab$FDR_immune_features[k] <- p.adjust(mv_tab$p[k], "BH")
}
mv_tab$HR_95CI <- fmt_hr(mv_tab$HR, mv_tab$CI_lower, mv_tab$CI_upper)
mv_tab$adjusted_for <- "age (per 10 y), stage, PAM50"
write.csv(r3(mv_tab), file.path(out_dir, "multivariable_cox_results.csv"), row.names = FALSE)
full_tab <- do.call(rbind, full_rows)
full_tab$HR_95CI <- fmt_hr(full_tab$HR, full_tab$CI_lower, full_tab$CI_upper)
write.csv(r3(full_tab), file.path(out_dir, "multivariable_cox_full_models.csv"), row.names = FALSE)

#Proportional hazards and linearity
ph_tab <- do.call(rbind, ph_rows)
ph_tab$flag <- ifelse(ph_tab$p < 0.05, "PH possibly violated", "")
write.csv(r3(ph_tab), file.path(out_dir, "proportional_hazards_checks.csv"), row.names = FALSE)

png(file.path(out_dir, "schoenfeld_residuals_immune_score_PFI.png"), width = 7, height = 5, units = "in", res = 300)
zp <- cox.zph(main_fits[["PFI immune_score_sd"]])
plot(zp, var = "immune_score_sd", resid = TRUE, se = TRUE, col = "#B2182B", lwd = 2,
     xlab = "Time (years)", ylab = "Beta(t) for immune score (per SD)",
     main = sprintf("PFI: immune score effect over time (PH test p = %s)",
                    format.pval(zp$table["immune_score_sd", "p"], digits = 2)))
abline(h = coef(main_fits[["PFI immune_score_sd"]])["immune_score_sd"], lty = 2, col = "grey40")
abline(h = 0, lty = 3)
dev.off()

lin_rows <- list()
for (e in c("PFI", "OS")) {
  d <- cc_data(ep_data(e), covars_main)
  adj <- paste(covars_main, collapse = " + ")
  basis <- ns(d$immune_score_sd, df = 3); d$imm_ns <- basis
  f_lin <- cox_fit(paste("immune_score_sd +", adj), d)
  f_ns  <- cox_fit(paste("imm_ns +", adj), d)
  lin_rows[[e]] <- data.frame(endpoint = e, n = f_lin$n, events = f_lin$nevent,
                              AIC_linear = AIC(f_lin), AIC_spline = AIC(f_ns),
                              nonlinearity_LRT_p = pchisq(2 * (f_ns$loglik[2] - f_lin$loglik[2]), 2, lower.tail = FALSE))
  {
    grid <- seq(quantile(d$immune_score_sd, 0.02), quantile(d$immune_score_sd, 0.98), length.out = 100)
    B <- predict(basis, grid) - matrix(predict(basis, 0), length(grid), 3, byrow = TRUE)
    k <- startsWith(names(coef(f_ns)), "imm_ns")
    lp <- drop(B %*% coef(f_ns)[k]); se <- sqrt(rowSums((B %*% vcov(f_ns)[k, k]) * B))
    png(file.path(out_dir, paste0("immune_score_spline_hazard_ratio_", e, ".png")), width = 7, height = 5, units = "in", res = 300)
    plot(grid, exp(lp), type = "l", lwd = 2, log = "y", ylim = range(exp(c(lp - 1.96 * se, lp + 1.96 * se))),
         xlab = "Immune score (SD units; 0 = cohort mean)", ylab = "Adjusted hazard ratio (log scale)",
         main = paste0(e, ": immune score as a spline (adjusted for age, stage, PAM50)"))
    polygon(c(grid, rev(grid)), exp(c(lp - 1.96 * se, rev(lp + 1.96 * se))), col = adjustcolor("#B2182B", 0.2), border = NA)
    abline(h = 1, lty = 2); rug(d$immune_score_sd[d$immune_score_sd >= min(grid) & d$immune_score_sd <= max(grid)], col = "grey60")
    lines(grid, exp(coef(f_lin)["immune_score_sd"] * grid), lty = 3, col = "grey30")
    legend("topright", c("Spline (95% CI)", "Linear fit"), lty = c(1, 3), col = c("black", "grey30"), bty = "n")
    dev.off()
  }
}
lin_tab <- do.call(rbind, lin_rows)
write.csv(r3(lin_tab), file.path(out_dir, "immune_score_linearity_check.csv"), row.names = FALSE)

#Sensitivity Analysis (PFI and OS)
sens_specs <- list(
  "Main model"                          = list(adj = "age10 + stage + pam50", sub = NULL, censor = NULL),
  "Plus histology"                      = list(adj = "age10 + stage + pam50 + histology_group", sub = NULL, censor = NULL),
  "Censored at 10 years"                = list(adj = "age10 + stage + pam50", sub = NULL, censor = 10),
  "Excluding stage IV"                  = list(adj = "age10 + stage + pam50", sub = quote(stage != "IV"), censor = NULL),
  "Stratified by PAM50 (not adjusted)"  = list(adj = "age10 + stage + strata(pam50)", sub = NULL, censor = NULL),
  "Receptor subtype instead of PAM50"   = list(adj = "age10 + stage + receptor_subtype", sub = NULL, censor = NULL),
  "Stratified by stage"                 = list(adj = "age10 + strata(stage) + pam50", sub = NULL, censor = NULL),
  # Relaxes every covariate flagged by cox.zph: age effect allowed to change with log(time),
  # stage and PAM50 as strata (separate baseline hazards) instead of adjusted terms
  "Age time-varying; stratified by stage and PAM50" =
    list(adj = "age10 + tt(age10) + strata(stage) + strata(pam50)", sub = NULL, censor = NULL))
sens_rows <- list()
for (e in c("PFI", "OS")) for (sn in names(sens_specs)) {
  sp <- sens_specs[[sn]]
  vars <- all.vars(as.formula(paste("~", sp$adj)))
  d <- ep_data(e, censor_at = sp$censor)
  if (!is.null(sp$sub)) d <- d[which(eval(sp$sub, d)), ]
  d <- cc_data(d, vars)
  for (v in c("immune_score_sd", "cluster_ref")) {
    fit <- if (grepl("tt\\(", sp$adj))
      coxph(as.formula(paste("Surv(time, status) ~", v, "+", sp$adj)), data = d, ties = "efron",
            tt = function(x, t, ...) x * log(t)) else cox_fit(paste(v, "+", sp$adj), d)
    r <- hr_rows(fit, v)
    sens_rows[[length(sens_rows) + 1]] <- data.frame(endpoint = e, analysis = sn,
                                                     exposure = ifelse(v == "immune_score_sd", "Immune score (per SD)", "Immune cluster (ref Depleted)"),
                                                     r[, -1], n = fit$n, events = fit$nevent)
  }
}
sens_tab <- do.call(rbind, sens_rows)
sens_tab$level[sens_tab$level == ""] <- "per SD"
sens_tab$HR_95CI <- fmt_hr(sens_tab$HR, sens_tab$CI_lower, sens_tab$CI_upper)
write.csv(r3(sens_tab), file.path(out_dir, "sensitivity_cox_results.csv"), row.names = FALSE)

#Effect by subtype and Basel-like Hot vs Cold
# Subgroup models use stage I-II vs III-IV to keep events per term acceptable.
sub_rows <- list(); int_rows <- list()
for (e in names(endpoints)) {
  d <- cc_data(ep_data(e), c("age10", "stage_2", "subtype3"))
  f_int <- cox_fit("immune_score_sd * subtype3 + age10 + stage_2", d)
  f_add <- cox_fit("immune_score_sd + subtype3 + age10 + stage_2", d)
  int_rows[[e]] <- data.frame(endpoint = e, n = f_int$n, events = f_int$nevent,
                              interaction_LRT_p = pchisq(2 * (f_int$loglik[2] - f_add$loglik[2]), 2, lower.tail = FALSE))
  for (s in levels(d$subtype3)) {
    ds <- droplevels(d[d$subtype3 == s, ])
    rhs <- if (nlevels(ds$stage_2) > 1) "immune_score_sd + age10 + stage_2" else "immune_score_sd + age10"
    fit <- cox_fit(rhs, ds)
    r <- hr_rows(fit, "immune_score_sd")
    sub_rows[[length(sub_rows) + 1]] <- data.frame(endpoint = e, subtype = s, r[, c("HR", "CI_lower", "CI_upper", "p")],
                                                   n = fit$n, events = fit$nevent, interaction_LRT_p = int_rows[[e]]$interaction_LRT_p)
  }
}
sub_tab <- do.call(rbind, sub_rows)
sub_tab$HR_95CI <- fmt_hr(sub_tab$HR, sub_tab$CI_lower, sub_tab$CI_upper)
sub_tab$adjusted_for <- "age (per 10 y), stage I-II vs III-IV"
write.csv(r3(sub_tab), file.path(out_dir, "immune_score_effect_by_subtype.csv"), row.names = FALSE)

bas_rows <- list()
for (e in names(endpoints)) {
  d <- cc_data(ep_data(e), c("basal_hot_cold", "age10", "stage_2"))
  for (mdl in c("Unadjusted", "Adjusted (age, stage I-II vs III-IV)")) {
    rhs <- if (mdl == "Unadjusted") "basal_hot_cold" else "basal_hot_cold + age10 + stage_2"
    fit <- cox_fit(rhs, d)
    r <- hr_rows(fit, "basal_hot_cold")
    bas_rows[[length(bas_rows) + 1]] <- data.frame(endpoint = e, model = mdl,
                                                   comparison = "Hot (Immune-enriched) vs Cold (Depleted)", r[, c("HR", "CI_lower", "CI_upper", "p")],
                                                   n = fit$n, events = fit$nevent)
  }
}
bas_tab <- do.call(rbind, bas_rows)
bas_tab$HR_95CI <- fmt_hr(bas_tab$HR, bas_tab$CI_lower, bas_tab$CI_upper)
write.csv(r3(bas_tab), file.path(out_dir, "basal_like_hot_vs_cold_cox.csv"), row.names = FALSE)

#Forest plot
forest <- function(df, labels, file, main, cols, xlab = "Adjusted hazard ratio (95% CI, log scale)") {
  eps <- unique(df$endpoint); nl <- length(labels)
  png(file.path(out_dir, file), width = 9, height = 1.2 + 0.45 * nl * length(eps) / 1.6, units = "in", res = 300)
  par(mar = c(4.5, 15, 3, 1))
  xr <- range(c(df$CI_lower, df$CI_upper), na.rm = TRUE); xr <- c(min(xr[1], 0.5), max(xr[2], 2))
  plot(NA, xlim = xr, ylim = c(0.5, nl + 0.5), log = "x", yaxt = "n", ylab = "", xlab = xlab)
  title(main, adj = 1, cex.main = 0.95)
  abline(v = 1, lty = 2, col = "grey50")
  axis(2, at = nl:1, labels = labels, las = 1, cex.axis = 0.8)
  off <- seq(0.25, -0.25, length.out = length(eps))   # first endpoint on top
  for (j in seq_along(eps)) {
    dj <- df[df$endpoint == eps[j], ]
    y <- (nl:1)[match(dj$label, labels)] + off[j]
    segments(dj$CI_lower, y, dj$CI_upper, y, col = cols[eps[j]], lwd = 1.8)
    points(dj$HR, y, pch = 19, col = cols[eps[j]], cex = 0.9)
  }
  legend("topright", legend = eps, col = cols[eps], pch = 19, lwd = 1.8, bty = "n", cex = 0.85)
  dev.off()
}
ep_cols <- c(PFI = "#B2182B", OS = "#2166AC", DFI = "#4DAF4A", DSS = "#984EA3")

fd <- mv_tab
fd$label <- ifelse(fd$exposure == "Immune cluster (ref Depleted)", paste(fd$level, "vs Depleted"), sub(" \\(per SD\\)", " (per SD)", fd$exposure))
forest(fd, labels = unique(fd$label), file = "multivariable_hazard_ratio_forest_plot.png", cols = ep_cols,
       main = "Immune features, adjusted for age, stage and PAM50")

fs <- sub_tab; fs$label <- fs$subtype
forest(fs, labels = levels(a$subtype3), file = "immune_score_by_subtype_forest_plot.png", cols = ep_cols,
       main = "Immune score (per SD) by subtype, adjusted for age and stage",
       xlab = "Hazard ratio (95% CI, log scale)")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_session_info_survival.txt"))

op <- options(width = 200)
show <- function(df, cols) print(df[, cols], row.names = FALSE)
cat("\n Kaplan-Meier: 5-year event-free % and log-rank p \n")
show(km_tab, c("grouping", "endpoint", "group", "n", "events", "event_free_5y_pct", "CI_lower_pct", "CI_upper_pct", "logrank_p"))
cat("\n Univariate Cox: immune features and clusters \n")
u <- uni_tab[uni_tab$variable %in% c(paste(cont_exposures, "(per SD)"), "Immune cluster (ref Depleted)", "Immune group (ref Immune-low)"), ]
show(r3(u), c("endpoint", "variable", "level", "HR_95CI", "p", "FDR_immune_features", "n", "events"))
cat("\nUnivariate Cox: clinical covariates (PFI and OS) \n")
u <- uni_tab[!feat_rows & !grepl("Immune", uni_tab$variable) & uni_tab$endpoint %in% c("PFI", "OS"), ]
show(r3(u), c("endpoint", "variable", "level", "HR_95CI", "p", "variable_LRT_p", "n", "events"))
cat("\n Multivariable Cox (age, stage, PAM50) \n")
show(r3(mv_tab), c("endpoint", "exposure", "level", "HR_95CI", "p", "FDR_immune_features", "exposure_LRT_p", "n", "events", "events_per_term"))
cat("\n Full main models (PFI and OS) \n")
show(r3(full_tab[full_tab$endpoint %in% c("PFI", "OS"), ]), c("endpoint", "model", "term", "HR_95CI", "p"))
cat("\n Proportional hazards (cox.zph) \n")
show(r3(ph_tab), c("endpoint", "model", "term", "chisq", "df", "p", "flag"))
cat("\n Linearity of immune score \n"); print(r3(lin_tab), row.names = FALSE)
cat("\n Sensitivity analyses \n")
show(r3(sens_tab), c("endpoint", "analysis", "exposure", "level", "HR_95CI", "p", "n", "events"))
cat("\n Immune score by subtype \n")
show(r3(sub_tab), c("endpoint", "subtype", "HR_95CI", "p", "n", "events", "interaction_LRT_p"))
cat("\n--- Basal-like hot vs cold ---\n")
show(r3(bas_tab), c("endpoint", "model", "HR_95CI", "p", "n", "events"))
cat("\nOutputs written to:", out_dir, "\n")
options(op)
