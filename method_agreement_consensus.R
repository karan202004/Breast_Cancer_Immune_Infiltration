project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
input_dir   <- file.path(project_dir, "results", "deconvolution_QC")
output_dir  <- file.path(project_dir, "results", "method_comparison_consensus")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

primary_methods   <- c("quantiseq", "EPIC", "MCPcounter", "xCell")  # consensus methods
min_methods       <- 3     # consensus needs at least this many primary methods
min_median_rho    <- 0.4   # ...and at least this median pairwise Spearman rho
cib_p_cutoff      <- 0.05  # CIBERSORT sensitivity subset

method_files <- c(
  quantiseq  = "aligned_quantiseq_scores.rds",
  EPIC       = "aligned_EPIC_scores.rds",
  MCPcounter = "aligned_MCPcounter_scores.rds",
  xCell      = "aligned_xCell_scores.rds",
  CIBERSORT  = "aligned_CIBERSORT_scores.rds"
)

scores <- list()
for (m in names(method_files)) {
  p <- file.path(input_dir, method_files[[m]])
  if (!file.exists(p)) stop("Missing file: ", p)
  x <- as.data.frame(readRDS(p), check.names = FALSE)
  if (!"sample_id" %in% names(x)) stop(m, ": no sample_id column")
  ids <- as.character(x$sample_id)
  x$sample_id <- NULL
  if (anyDuplicated(ids)) stop(m, ": duplicated sample IDs")
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  scores[[m]] <- mat
}

ref_ids <- rownames(scores[[1]])
for (m in names(scores)) {
  if (!identical(rownames(scores[[m]]), ref_ids))
    stop(m, ": sample order differs from ", names(scores)[1], " — rerun the QC step")
}
n_samples <- length(ref_ids)
if (n_samples != 1106) stop("Expected 1106 samples, found ", n_samples)

# CIBERSORT P-values for the sensitivity subset
cib_diag_path <- file.path(input_dir, "aligned_CIBERSORT_diagnostics.rds")
cib_sig <- NULL
if (file.exists(cib_diag_path)) {
  cd <- readRDS(cib_diag_path)
  if (identical(as.character(cd$sample_id), ref_ids) && "P-value" %in% names(cd))
    cib_sig <- cd[["P-value"]] < cib_p_cutoff
}

#Record every available column
available <- do.call(rbind, lapply(names(scores), function(m)
  data.frame(method = m, column = colnames(scores[[m]]), stringsAsFactors = FALSE)))
write.csv(available, file.path(output_dir, "00_available_columns.csv"), row.names = FALSE)

cat("Columns available per method:\n")
for (m in names(scores)) cat("\n", m, ":\n  ", paste(colnames(scores[[m]]), collapse = " | "), "\n")

#Cell-type mapping
map <- list(
  CD8_T_cells = list(
    quantiseq  = list(c("T cell CD8+", "T.cells.CD8", "T cells CD8")),
    EPIC       = list(c("T cell CD8+", "CD8_Tcells", "CD8 T cells")),
    MCPcounter = list(c("T cell CD8+", "CD8 T cells")),
    xCell      = list(c("T cell CD8+", "CD8+ T-cells")),
    CIBERSORT  = list("T cells CD8")
  ),
  CD4_T_cells = list(   # conventional (non-regulatory) CD4 where possible
    quantiseq  = list(c("T cell CD4+ (non-regulatory)", "T.cells.CD4", "T cells CD4")),
    EPIC       = list(c("T cell CD4+", "CD4_Tcells", "CD4 T cells")),
    xCell      = list(c("T cell CD4+ (non-regulatory)", "CD4+ T-cells")),
    CIBERSORT  = list("T cells CD4 naive", "T cells CD4 memory resting",
                      "T cells CD4 memory activated", "T cells follicular helper")
  ),
  Tregs = list(
    quantiseq  = list(c("T cell regulatory (Tregs)", "Tregs")),
    xCell      = list(c("T cell regulatory (Tregs)", "Tregs")),
    CIBERSORT  = list("T cells regulatory (Tregs)")
  ),
  B_cells = list(
    quantiseq  = list(c("B cell", "B.cells", "B cells")),
    EPIC       = list(c("B cell", "Bcells", "B cells")),
    MCPcounter = list(c("B cell", "B lineage")),
    xCell      = list(c("B cell", "B-cells")),
    CIBERSORT  = list("B cells naive", "B cells memory")
  ),
  NK_cells = list(
    quantiseq  = list(c("NK cell", "NK.cells", "NK cells")),
    EPIC       = list(c("NK cell", "NKcells", "NK cells")),
    MCPcounter = list(c("NK cell", "NK cells")),
    xCell      = list(c("NK cell", "NK cells")),
    CIBERSORT  = list("NK cells resting", "NK cells activated")
  ),
  Macrophages = list(
    quantiseq  = list(c("Macrophage M1", "Macrophages.M1", "Macrophages M1"),
                      c("Macrophage M2", "Macrophages.M2", "Macrophages M2")),
    EPIC       = list(c("Macrophage", "Macrophages")),
    MCPcounter = list(c("Macrophage/Monocyte", "Monocytic lineage")),   # approximate
    xCell      = list(c("Macrophage", "Macrophages")),
    CIBERSORT  = list("Macrophages M0", "Macrophages M1", "Macrophages M2")
  ),
  Monocytes = list(
    quantiseq  = list(c("Monocyte", "Monocytes")),
    MCPcounter = list(c("Monocyte", "Monocytic lineage")),              # approximate
    xCell      = list(c("Monocyte", "Monocytes")),
    CIBERSORT  = list("Monocytes")
  ),
  Myeloid_DC = list(
    quantiseq  = list(c("Myeloid dendritic cell", "Dendritic.cells", "Dendritic cells")),
    MCPcounter = list(c("Myeloid dendritic cell", "Myeloid dendritic cells")),
    xCell      = list(c("Myeloid dendritic cell", "DC")),
    CIBERSORT  = list("Dendritic cells resting", "Dendritic cells activated")
  ),
  Neutrophils = list(
    quantiseq  = list(c("Neutrophil", "Neutrophils")),
    MCPcounter = list(c("Neutrophil", "Neutrophils")),
    xCell      = list(c("Neutrophil", "Neutrophils")),
    CIBERSORT  = list("Neutrophils")
  ),
  # Stromal populations: sanity check (methods usually agree well here)
  CAFs = list(
    EPIC       = list(c("Cancer associated fibroblast", "CAFs")),
    MCPcounter = list(c("Cancer associated fibroblast", "Fibroblasts")),
    xCell      = list(c("Cancer associated fibroblast", "Fibroblasts"))
  ),
  Endothelial = list(
    EPIC       = list(c("Endothelial cell", "Endothelial")),
    MCPcounter = list(c("Endothelial cell", "Endothelial cells")),
    xCell      = list(c("Endothelial cell", "Endothelial cells"))
  )
)

approximate <- list(Macrophages = "MCPcounter", Monocytes = "MCPcounter")

nn <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

resolve <- function(method, comps) {
  cn  <- colnames(scores[[method]])
  key <- nn(cn)
  used <- character(0)
  for (alts in comps) {
    hit <- match(nn(alts), key)
    hit <- hit[!is.na(hit)]
    if (!length(hit)) return(list(ok = FALSE, missing = alts[1], used = used))
    used <- c(used, cn[hit[1]])
  }
  list(ok = TRUE, v = rowSums(scores[[method]][, used, drop = FALSE]), used = used)
}

#Resolve the mapping
feat_vec <- list()
mapping  <- list()

for (f in names(map)) {
  for (m in names(scores)) {
    if (is.null(map[[f]][[m]])) {
      mapping[[length(mapping) + 1]] <- data.frame(
        feature = f, method = m, status = "no equivalent population",
        columns_used = "", stringsAsFactors = FALSE)
      next
    }
    r <- resolve(m, map[[f]][[m]])
    if (!r$ok) {
      status <- paste0("NOT FOUND (missing '", r$missing, "')")
    } else if (sd(r$v) == 0) {
      status <- "constant — excluded"
    } else {
      feat_vec[[f]][[m]] <- r$v
      status <- if (m %in% approximate[[f]]) "matched (approximate)" else "matched"
    }
    mapping[[length(mapping) + 1]] <- data.frame(
      feature = f, method = m, status = status,
      columns_used = paste(r$used, collapse = " + "), stringsAsFactors = FALSE)
  }
}
mapping <- do.call(rbind, mapping)
write.csv(mapping, file.path(output_dir, "01_cell_type_mapping.csv"), row.names = FALSE)

not_found <- mapping[grepl("NOT FOUND", mapping$status), ]
if (nrow(not_found)) {
  cat("\n*** WARNING: expected columns not found — check 00_available_columns.csv ***\n")
  print(not_found, row.names = FALSE)
}

##Pairwise Spearman agreement
pairs_out <- list()
for (f in names(feat_vec)) {
  ms <- names(feat_vec[[f]])
  if (length(ms) < 2) next
  for (pr in combn(ms, 2, simplify = FALSE)) {
    a <- feat_vec[[f]][[pr[1]]]
    b <- feat_vec[[f]][[pr[2]]]
    ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE))
    
    rho_sig <- NA_real_
    if ("CIBERSORT" %in% pr && !is.null(cib_sig)) {
      rho_sig <- suppressWarnings(cor(a[cib_sig], b[cib_sig], method = "spearman"))
    }
    
    pairs_out[[length(pairs_out) + 1]] <- data.frame(
      feature = f, method_1 = pr[1], method_2 = pr[2],
      n = length(a), spearman_rho = round(unname(ct$estimate), 3),
      p_value = signif(ct$p.value, 3),
      rho_CIBERSORT_P_lt_0.05 = round(rho_sig, 3),
      stringsAsFactors = FALSE)
  }
}
pairs_df <- if (length(pairs_out)) do.call(rbind, pairs_out) else data.frame()
write.csv(pairs_df, file.path(output_dir, "02_pairwise_method_correlations.csv"), row.names = FALSE)

##Heatmaps (one per population)
for (f in unique(pairs_df$feature)) {
  tmp <- pairs_df[pairs_df$feature == f, ]
  ms  <- unique(c(tmp$method_1, tmp$method_2))
  mat <- matrix(NA_real_, length(ms), length(ms), dimnames = list(ms, ms))
  diag(mat) <- 1
  for (i in seq_len(nrow(tmp))) {
    mat[tmp$method_1[i], tmp$method_2[i]] <- tmp$spearman_rho[i]
    mat[tmp$method_2[i], tmp$method_1[i]] <- tmp$spearman_rho[i]
  }
  pheatmap::pheatmap(
    mat, cluster_rows = FALSE, cluster_cols = FALSE,
    display_numbers = TRUE, number_format = "%.2f",
    breaks = seq(-1, 1, length.out = 101),
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = paste("Method agreement (Spearman):", f),
    filename = file.path(output_dir, paste0("03_heatmap_", f, ".png")),
    width = 5.5, height = 4.5)
}

##Agreement overview per population
overview <- do.call(rbind, lapply(unique(pairs_df$feature), function(f) {
  tmp  <- pairs_df[pairs_df$feature == f, ]
  prim <- tmp[tmp$method_1 %in% primary_methods & tmp$method_2 %in% primary_methods, ]
  data.frame(
    feature = f,
    methods = length(unique(c(tmp$method_1, tmp$method_2))),
    comparisons = nrow(tmp),
    min_rho = min(tmp$spearman_rho),
    median_rho = median(tmp$spearman_rho),
    max_rho = max(tmp$spearman_rho),
    median_rho_primary = if (nrow(prim)) median(prim$spearman_rho) else NA,
    stringsAsFactors = FALSE)
}))
write.csv(overview, file.path(output_dir, "04_method_agreement_overview.csv"), row.names = FALSE)

int_transform <- function(v) qnorm((rank(v, ties.method = "average") - 0.5) / length(v))

consensus <- list()
cons_summary <- list()

for (f in names(map)) {
  ms <- intersect(primary_methods, names(feat_vec[[f]]))
  prim <- pairs_df[pairs_df$feature == f &
                     pairs_df$method_1 %in% ms & pairs_df$method_2 %in% ms, ]
  med_rho <- if (nrow(prim)) median(prim$spearman_rho) else NA
  min_rho <- if (nrow(prim)) min(prim$spearman_rho) else NA
  
  if (length(ms) < min_methods) {
    decision <- sprintf("not created: only %d primary method(s)", length(ms))
  } else if (med_rho < min_median_rho) {
    decision <- sprintf("not created: median rho %.2f < %.2f", med_rho, min_median_rho)
  } else {
    z <- sapply(ms, function(m) int_transform(feat_vec[[f]][[m]]))
    consensus[[f]] <- rowMeans(z)
    decision <- "created"
  }
  
  cons_summary[[length(cons_summary) + 1]] <- data.frame(
    feature = f, primary_methods_matched = paste(ms, collapse = "; "),
    n_methods = length(ms),
    median_rho_primary = round(med_rho, 3), lowest_pair_rho = round(min_rho, 3),
    decision = decision, stringsAsFactors = FALSE)
}
cons_summary <- do.call(rbind, cons_summary)
write.csv(cons_summary, file.path(output_dir, "05_consensus_feature_summary.csv"), row.names = FALSE)

if (length(consensus)) {
  cons_df <- data.frame(sample_id = ref_ids, as.data.frame(consensus), check.names = FALSE)
  write.csv(cons_df, file.path(output_dir, "06_consensus_immune_features.csv"), row.names = FALSE)
  saveRDS(cons_df, file.path(output_dir, "06_consensus_immune_features.rds"))
}

cat("Samples:", n_samples, "\n")
cat("CIBERSORT sensitivity subset (P <", cib_p_cutoff, "):",
    if (is.null(cib_sig)) "not available" else sum(cib_sig), "\n\n")

cat("Mapping:\n");            print(mapping[mapping$status != "no equivalent population", ], row.names = FALSE)
cat("\nAgreement overview:\n"); print(overview, row.names = FALSE)
cat("\nAll pairwise correlations:\n"); print(pairs_df, row.names = FALSE)
cat("\nConsensus decisions:\n"); print(cons_summary, row.names = FALSE)
cat("\nOutputs saved to:\n", output_dir, "\n")
