# Survival outcome data check


suppressPackageStartupMessages(library(survival))

project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
out_dir <- file.path(project_dir, "results", "survival_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

m <- readRDS(file.path(project_dir, "metadata", "clinical_master_table.rds"))
a <- m[m$in_analysis_set, ]
lab_order <- unique(a$cluster_label[order(a$cluster)])
a$cluster_name <- factor(sub("^C[0-9]+: ", "", a$cluster_label), levels = sub("^C[0-9]+: ", "", lab_order))

endpoints <- c(OS = "Overall survival", PFI = "Progression-free interval",
               DFI = "Disease-free interval", DSS = "Disease-specific survival")


ep_summary <- do.call(rbind, lapply(names(endpoints), function(e) {
  st <- a[[e]]; tm <- a[[paste0(e, ".time")]]
  ok <- !is.na(st) & !is.na(tm)
  bad_time <- ok & tm <= 0
  use <- ok & tm > 0
  rev_km <- survfit(Surv(tm[use] / 365.25, 1 - st[use]) ~ 1)      # reverse KM: censoring as the event
  km <- survfit(Surv(tm[use] / 365.25, st[use]) ~ 1)
  s5 <- summary(km, times = 5, extend = TRUE)
  data.frame(endpoint = e, description = endpoints[[e]],
             n_available = sum(ok), pct_available = round(100 * mean(ok), 1),
             time_zero_or_negative = sum(bad_time),
             n_usable = sum(use), events = sum(st[use] == 1), pct_events = round(100 * mean(st[use] == 1), 1),
             median_follow_up_years = round(unname(summary(rev_km)$table["median"]), 1),
             max_follow_up_years = round(max(tm[use]) / 365.25, 1),
             events_after_5y = sum(st[use] == 1 & tm[use] > 5 * 365.25),
             km_5y_event_free_pct = round(100 * s5$surv, 1),
             model_terms_allowed_10_epv = floor(sum(st[use] == 1) / 10))
}))
write.csv(ep_summary, file.path(out_dir, "survival_endpoint_completeness.csv"), row.names = FALSE)

#Events by group (are there enough events in each group?)
event_tab <- function(group, gname) {
  do.call(rbind, lapply(names(endpoints), function(e) {
    st <- a[[e]]; tm <- a[[paste0(e, ".time")]]
    use <- !is.na(st) & !is.na(tm) & tm > 0 & !is.na(group)
    g <- droplevels(as.factor(group[use]))
    data.frame(grouping = gname, endpoint = e, group = levels(g),
               n = as.integer(table(g)), events = as.integer(tapply(st[use], g, sum)))
  }))
}
by_group <- rbind(event_tab(a$cluster_name, "Immune cluster"),
                  event_tab(a$pam50, "PAM50"),
                  event_tab(a$stage, "Stage"),
                  event_tab(a$histology_group, "Histology"),
                  event_tab(ifelse(a$pam50 %in% "Basal-like" & a$cluster_name %in% c("Immune-enriched", "Depleted"),
                                   as.character(a$cluster_name), NA), "Basal-like hot vs cold"))
write.csv(by_group, file.path(out_dir, "survival_events_by_group.csv"), row.names = FALSE)

#Candidate covariates: completeness among tumours with PFI
pfi_ok <- !is.na(a$PFI) & !is.na(a$PFI.time) & a$PFI.time > 0
covars <- c(age = "Age", stage = "Stage", pT = "T stage", pN = "N stage", pam50 = "PAM50",
            receptor_subtype = "Receptor subtype", histology_group = "Histology",
            menopause_status = "Menopausal status", absolute_purity = "Tumour purity")
cov_tab <- data.frame(
  variable = unname(covars),
  n_available_with_PFI = sapply(names(covars), function(v) sum(!is.na(a[[v]][pfi_ok]))),
  pct_available = sapply(names(covars), function(v) round(100 * mean(!is.na(a[[v]][pfi_ok])), 1)),
  model_terms = sapply(names(covars), function(v) {
    x <- a[[v]][pfi_ok]; if (is.numeric(x)) 1 else length(unique(na.omit(x))) - 1 }),
  row.names = NULL)
cc <- complete.cases(a[pfi_ok, c("age", "stage", "pam50")])
write.csv(cov_tab, file.path(out_dir, "survival_covariate_completeness.csv"), row.names = FALSE)

op <- options(width = 200)
cat("Tumours in analysis set:", nrow(a), "\n\n")
print(ep_summary, row.names = FALSE)
cat("\nEvents by group:\n")
print(reshape(by_group, idvar = c("grouping", "group"), timevar = "endpoint", direction = "wide"), row.names = FALSE)
cat("\nCandidate covariates (tumours with usable PFI):\n"); print(cov_tab, row.names = FALSE)
cat("\nTumours with PFI and complete age + stage + PAM50:", sum(cc), "of", sum(pfi_ok), "\n")
options(op)
