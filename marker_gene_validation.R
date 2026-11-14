out_dir <- file.path(project_dir, "results", "method_comparison_consensus")

primary_methods <- c("quantiseq", "EPIC", "MCPcounter", "xCell")
min_methods     <- 3
min_median_rho  <- 0.4

#Population-specific exclusions, each with a documented reason
consensus_exclude <- list(
  CD8_T_cells = c(EPIC = "negative correlation with all other methods (rho -0.08 to -0.23) and with CD8A/CD8B expression (rho -0.23)")
)

int_transform <- function(v) qnorm((rank(v, ties.method = "average") - 0.5) / length(v))

consensus <- list(); cons_summary <- list()
for (f in names(feat_vec)) {
  excl <- consensus_exclude[[f]]
  ms <- setdiff(intersect(primary_methods, names(feat_vec[[f]])), names(excl))
  prim <- pairs_df[pairs_df$feature == f &
                     pairs_df$method_1 %in% ms & pairs_df$method_2 %in% ms, ]
  med_rho <- if (nrow(prim)) median(prim$spearman_rho) else NA
  min_rho <- if (nrow(prim)) min(prim$spearman_rho) else NA
  
  decision <- if (length(ms) < min_methods) {
    sprintf("not created: only %d primary method(s)", length(ms))
  } else if (med_rho < min_median_rho) {
    sprintf("not created: median rho %.2f < %.2f", med_rho, min_median_rho)
  } else {
    consensus[[f]] <- rowMeans(sapply(ms, function(m) int_transform(feat_vec[[f]][[m]])))
    "created"
  }
  
  cons_summary[[length(cons_summary) + 1]] <- data.frame(
    feature = f, methods_used = paste(ms, collapse = "; "), n_methods = length(ms),
    excluded = if (length(excl)) paste0(names(excl), ": ", excl, collapse = " | ") else "",
    median_rho = round(med_rho, 3), lowest_pair_rho = round(min_rho, 3),
    decision = decision, stringsAsFactors = FALSE)
}
cons_summary <- do.call(rbind, cons_summary)
write.csv(cons_summary, file.path(out_dir, "05_consensus_feature_summary.csv"), row.names = FALSE)

cons_df <- data.frame(sample_id = ref_ids, as.data.frame(consensus), check.names = FALSE)
write.csv(cons_df, file.path(out_dir, "06_consensus_immune_features.csv"), row.names = FALSE)
saveRDS(cons_df, file.path(out_dir, "06_consensus_immune_features.rds"))

print(cons_summary[, c("feature", "methods_used", "median_rho", "lowest_pair_rho", "decision")], row.names = FALSE)
cat("\nConsensus features saved:", paste(names(consensus), collapse = ", "), "\n")
