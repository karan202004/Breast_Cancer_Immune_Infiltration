# ============================================================
# INTEGRATED QC — FIVE IMMUNE DECONVOLUTION METHODS
# TCGA-BRCA | 1,106 primary tumour samples
# Uses existing result files only. Does NOT rerun deconvolution.
#
# Status levels in the report:
#   FAIL   = structural problem (missing/duplicated samples, NA, Inf,
#            wrong number of CIBERSORT cell types, unreadable file)
#   REVIEW = needs interpretation, not automatically wrong
#            (negative values, row sums off 1, zero-variance cell types,
#             all-zero samples, extra sample IDs)
#   INFO   = context only
# ============================================================

## 0. Paths and metadata --------------------------------------
project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
result_dir  <- file.path(project_dir, "results", "deconvolution_results")
qc_dir      <- file.path(project_dir, "results", "deconvolution_QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

meta_path <- file.path(project_dir, "metadata", "primary_tumor_1106_metadata.csv")
if (!file.exists(meta_path)) stop("Metadata not found: ", meta_path)

metadata <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"sample_id" %in% colnames(metadata)) stop("Metadata has no 'sample_id' column.")
expected_samples <- as.character(metadata$sample_id)

if (anyNA(expected_samples) || any(expected_samples == ""))
  stop("Empty or NA sample IDs in metadata.")
if (anyDuplicated(expected_samples))
  stop("Duplicate sample IDs in metadata.")
if (length(expected_samples) != 1106)
  stop("Expected 1106 samples in metadata, found ", length(expected_samples))

## 1. Explicit result files (no pattern matching) -------------
result_files <- c(
  quantiseq  = "quantiseq_results_1106.rds",
  EPIC       = "epic_results_1106.rds",
  MCPcounter = "mcp_counter_results_1106.rds",
  xCell      = "xcell_results_1106.rds",
  CIBERSORT  = "cibersort_LM22_results_1106.rds"
)
cib_fit_file <- "cibersort_fit_metrics_1106.rds"   # optional, used only if
# diagnostics aren't embedded

missing_files <- result_files[!file.exists(file.path(result_dir, result_files))]
if (length(missing_files) > 0) {
  stop("Missing result files:\n",
       paste0("  ", names(missing_files), ": ", missing_files, collapse = "\n"),
       "\n\nFiles actually in ", result_dir, ":\n",
       paste0("  ", list.files(result_dir), collapse = "\n"))
}
cat("Result files used:\n"); print(result_files); cat("\n")

## 2. Helpers --------------------------------------------------
flags <- list()
add_flag <- function(method, check, level, detail) {
  flags[[length(flags) + 1]] <<- data.frame(
    method = method, check = check, level = level, detail = detail,
    stringsAsFactors = FALSE)
}

norm_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

# TCGA barcodes sometimes get "-" turned into "." by check.names.
# Swap back only if it improves matching to the metadata.
fix_ids <- function(ids) {
  ids <- as.character(ids)
  alt <- gsub(".", "-", ids, fixed = TRUE)
  if (sum(expected_samples %in% alt) > sum(expected_samples %in% ids)) alt else ids
}

# CIBERSORT fit diagnostics (kept out of the cell-score matrix)
diag_map <- c(pvalue = "P-value", pval = "P-value",
              correlation = "Correlation", rmse = "RMSE",
              absolutescoresigscore = "Absolute score (sig.score)")

## 3. Loader: verifies orientation, IDs, duplicates -----------
load_method <- function(method, filename) {
  
  x <- readRDS(file.path(result_dir, filename))
  raw_class <- paste(class(x), collapse = "/")
  if (is.null(dim(x))) stop("not a matrix/data frame (class: ", raw_class, ")")
  raw_dim <- paste(dim(x), collapse = " x ")
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  
  label_col <- intersect(c("cell_type", "Cell.Type", "celltype",
                           "CellType", "name", "label"), colnames(x))
  id_col    <- intersect(c("Mixture", "sample_id", "SampleID",
                           "Sample", "sample"), colnames(x))
  
  if (length(label_col) > 0) {
    # immunedeconv style: cell types in rows, samples in columns
    features <- as.character(x[[label_col[1]]])
    x[[label_col[1]]] <- NULL
    ids <- fix_ids(colnames(x))
    samples_in_rows <- FALSE
  } else {
    # Decide orientation from where the sample IDs actually are
    row_ids <- if (length(id_col)) as.character(x[[id_col[1]]]) else rownames(x)
    if (length(id_col)) x[[id_col[1]]] <- NULL
    row_ids <- fix_ids(row_ids)
    col_ids <- fix_ids(colnames(x))
    hits_r <- sum(expected_samples %in% row_ids)
    hits_c <- sum(expected_samples %in% col_ids)
    if (hits_r == 0 && hits_c == 0) {
      stop("sample IDs not found in rows or columns. First row names: ",
           paste(head(row_ids, 3), collapse = ", "), " | first column names: ",
           paste(head(col_ids, 3), collapse = ", "))
    }
    samples_in_rows <- hits_r >= hits_c
    if (samples_in_rows) { ids <- row_ids; features <- colnames(x) }
    else                 { ids <- col_ids; features <- rownames(x) }
  }
  
  # Numeric conversion, counting values that fail to convert
  coercion_na <- 0
  for (j in seq_along(x)) {
    v <- x[[j]]
    was_na <- is.na(v) | trimws(as.character(v)) %in% c("", "NA", "NaN")
    num <- suppressWarnings(as.numeric(as.character(v)))
    coercion_na <- coercion_na + sum(is.na(num) & !was_na)
    x[[j]] <- num
  }
  mat <- as.matrix(x)
  if (!samples_in_rows) mat <- t(mat)
  rownames(mat) <- ids
  colnames(mat) <- features
  
  # Structural checks — stop this method on failure
  if (anyNA(features) || any(features == ""))
    stop("empty cell-type labels present")
  if (anyDuplicated(features))
    stop("duplicated cell-type labels: ",
         paste(unique(features[duplicated(features)]), collapse = ", "))
  dup_ids <- unique(ids[duplicated(ids)])
  if (length(dup_ids))
    stop(length(dup_ids), " duplicated sample IDs in result, e.g. ",
         paste(head(dup_ids, 3), collapse = ", "))
  miss <- setdiff(expected_samples, ids)
  if (length(miss))
    stop(length(miss), " expected samples missing, e.g. ",
         paste(head(miss, 3), collapse = ", "))
  extra <- setdiff(ids, expected_samples)
  
  mat <- mat[expected_samples, , drop = FALSE]   # metadata order
  
  # Separate fit diagnostics from cell scores
  diagnostics <- NULL
  d_idx <- which(norm_name(colnames(mat)) %in% names(diag_map))
  if (length(d_idx)) {
    diagnostics <- as.data.frame(mat[, d_idx, drop = FALSE], check.names = FALSE)
    colnames(diagnostics) <- unname(diag_map[norm_name(colnames(diagnostics))])
    mat <- mat[, -d_idx, drop = FALSE]
  }
  
  # CIBERSORT: fall back to separate fit-metrics file if needed
  if (method == "CIBERSORT" && is.null(diagnostics) &&
      file.exists(file.path(result_dir, cib_fit_file))) {
    fd <- as.data.frame(readRDS(file.path(result_dir, cib_fit_file)),
                        check.names = FALSE)
    fid <- if ("sample_id" %in% colnames(fd)) fd$sample_id else rownames(fd)
    fid <- fix_ids(fid)
    if (all(expected_samples %in% fid) && !anyDuplicated(fid)) {
      rownames(fd) <- fid
      fd$sample_id <- NULL
      keep <- norm_name(colnames(fd)) %in% names(diag_map)
      fd <- fd[expected_samples, keep, drop = FALSE]
      colnames(fd) <- unname(diag_map[norm_name(colnames(fd))])
      fd[] <- lapply(fd, function(z) suppressWarnings(as.numeric(as.character(z))))
      diagnostics <- fd
    } else {
      add_flag(method, "fit metrics file", "REVIEW",
               "cibersort_fit_metrics_1106.rds sample IDs do not align; ignored")
    }
  }
  
  list(
    scores = mat,
    diagnostics = diagnostics,
    info = data.frame(
      method = method, file = filename, loaded = TRUE,
      object_class = raw_class, raw_dim = raw_dim,
      samples_were_in = if (samples_in_rows) "rows" else "columns",
      n_extra_ids_dropped = length(extra),
      values_failed_numeric_conversion = coercion_na,
      n_diagnostic_cols = if (is.null(diagnostics)) 0 else ncol(diagnostics),
      stringsAsFactors = FALSE)
  )
}

## 4. Load all methods (one failure doesn't stop the others) ---
method_data <- list()
load_log <- list()
for (m in names(result_files)) {
  res <- tryCatch(load_method(m, result_files[[m]]), error = function(e) e)
  if (inherits(res, "error")) {
    add_flag(m, "load/alignment", "FAIL", conditionMessage(res))
    load_log[[m]] <- data.frame(
      method = m, file = result_files[[m]], loaded = FALSE,
      object_class = NA, raw_dim = NA, samples_were_in = NA,
      n_extra_ids_dropped = NA, values_failed_numeric_conversion = NA,
      n_diagnostic_cols = NA, stringsAsFactors = FALSE)
  } else {
    method_data[[m]] <- res
    load_log[[m]] <- res$info
    if (res$info$n_extra_ids_dropped > 0)
      add_flag(m, "extra sample IDs", "REVIEW",
               paste(res$info$n_extra_ids_dropped,
                     "IDs not in metadata were dropped"))
    if (res$info$values_failed_numeric_conversion > 0)
      add_flag(m, "numeric conversion", "FAIL",
               paste(res$info$values_failed_numeric_conversion,
                     "values could not be converted to numbers"))
  }
}
load_log <- do.call(rbind, load_log)
write.csv(load_log, file.path(qc_dir, "00_load_log.csv"), row.names = FALSE)

## 5. Sample-level QC ------------------------------------------
sum_to_one_methods <- c("quantiseq", "EPIC", "CIBERSORT")  # diagnostic only
sum_tol <- 0.01
xcell_aggregates <- c("immunescore", "stromascore", "microenvironmentscore")

sample_qc <- do.call(rbind, lapply(names(method_data), function(m) {
  mat <- method_data[[m]]$scores
  
  n_na      <- sum(is.na(mat) & !is.nan(mat))
  n_nonfin  <- sum(is.infinite(mat) | is.nan(mat))
  n_neg     <- sum(mat < 0, na.rm = TRUE)
  all_na    <- sum(rowSums(is.na(mat)) == ncol(mat))
  all_zero  <- sum(rowSums(mat == 0, na.rm = TRUE) == ncol(mat))
  zero_var  <- colnames(mat)[apply(mat, 2, function(v) {
    f <- v[is.finite(v)]; length(f) < 2 || sd(f) == 0 })]
  
  rs <- rowSums(mat, na.rm = TRUE)
  is_prop <- m %in% sum_to_one_methods
  n_off_one <- if (is_prop) sum(abs(rs - 1) > sum_tol) else NA
  
  # Flags
  if (n_na > 0)     add_flag(m, "missing values", "FAIL", paste(n_na, "NA values"))
  if (n_nonfin > 0) add_flag(m, "non-finite values", "FAIL", paste(n_nonfin, "Inf/NaN values"))
  if (all_na > 0)   add_flag(m, "all-NA samples", "FAIL", paste(all_na, "samples"))
  if (n_neg > 0)    add_flag(m, "negative values", "REVIEW",
                             paste(n_neg, "negative values — check whether the method can return them"))
  if (all_zero > 0) add_flag(m, "all-zero samples", "REVIEW", paste(all_zero, "samples with every score = 0"))
  if (length(zero_var))
    add_flag(m, "zero-variance cell types", "REVIEW",
             paste(length(zero_var), "constant columns:",
                   paste(head(zero_var, 6), collapse = "; ")))
  if (is_prop && n_off_one > 0)
    add_flag(m, "row sums vs 1", "REVIEW",
             sprintf("%d samples outside 1 ± %.2f (range %.3f–%.3f); expected only for relative-fraction output incl. 'uncharacterized cell'",
                     n_off_one, sum_tol, min(rs), max(rs)))
  if (m == "CIBERSORT" && ncol(mat) != 22)
    add_flag(m, "LM22 cell types", "FAIL",
             paste("expected 22 cell-type columns, found", ncol(mat),
                   ":", paste(colnames(mat), collapse = "; ")))
  
  data.frame(
    method = m,
    samples = nrow(mat),
    expected_samples = length(expected_samples),
    sample_order_matches = identical(rownames(mat), expected_samples),
    score_columns = ncol(mat),
    missing_values = n_na,
    nonfinite_values = n_nonfin,
    negative_values = n_neg,
    all_NA_samples = all_na,
    all_zero_samples = all_zero,
    zero_variance_columns = length(zero_var),
    min_row_sum = if (is_prop) round(min(rs), 4) else NA,
    max_row_sum = if (is_prop) round(max(rs), 4) else NA,
    samples_row_sum_off_1 = n_off_one,
    stringsAsFactors = FALSE
  )
}))

## 6. Per-cell-type score summary ------------------------------
feature_type <- function(m, cell) {
  if (m == "xCell" && norm_name(cell) %in% xcell_aggregates) return("xCell aggregate score")
  if (grepl("uncharacteri", cell, ignore.case = TRUE)) return("non-immune / residual")
  "cell type"
}

score_summary <- do.call(rbind, lapply(names(method_data), function(m) {
  mat <- method_data[[m]]$scores
  do.call(rbind, lapply(colnames(mat), function(cell) {
    v <- mat[, cell]; f <- v[is.finite(v)]
    q <- if (length(f)) quantile(f, c(0, .25, .5, .75, 1), names = FALSE) else rep(NA, 5)
    data.frame(
      method = m, cell_type_or_score = cell,
      feature_type = feature_type(m, cell),
      n_samples = length(v), missing = sum(is.na(v)),
      nonfinite = sum(is.infinite(v) | is.nan(v)),
      minimum = q[1], q1 = q[2], median = q[3], mean = if (length(f)) mean(f) else NA,
      q3 = q[4], maximum = q[5], sd = if (length(f) > 1) sd(f) else NA,
      negative_values = sum(v < 0, na.rm = TRUE),
      zero_fraction = round(mean(v == 0, na.rm = TRUE), 4),
      stringsAsFactors = FALSE)
  }))
}))
write.csv(score_summary, file.path(qc_dir, "02_integrated_score_summary.csv"), row.names = FALSE)

## 7. xCell follow-up ------------------------------------------
# Convergence warnings from the original run are NOT stored in the RDS.
# This checks the outputs for their likely footprint instead.
if (!is.null(method_data$xCell)) {
  xm <- method_data$xCell$scores
  agg <- norm_name(colnames(xm)) %in% xcell_aggregates
  cells <- xm[, !agg, drop = FALSE]
  frac_zero <- rowMeans(cells == 0, na.rm = TRUE)
  
  xcell_samples <- data.frame(
    sample_id = rownames(xm),
    n_celltypes = ncol(cells),
    zero_fraction = round(frac_zero, 4),
    all_celltypes_zero = frac_zero == 1,
    xm[, agg, drop = FALSE],
    check.names = FALSE)
  write.csv(xcell_samples, file.path(qc_dir, "06_xCell_sample_checks.csv"), row.names = FALSE)
  
  add_flag("xCell", "aggregate columns", "INFO",
           paste(sum(agg), "aggregate scores (ImmuneScore/StromaScore/MicroenvironmentScore) kept separate from",
                 ncol(cells), "cell types in summary"))
  add_flag("xCell", "per-sample zero fraction", "INFO",
           sprintf("median %.2f, max %.2f; %d samples with >=75%% zero cell types",
                   median(frac_zero), max(frac_zero), sum(frac_zero >= 0.75)))
  add_flag("xCell", "convergence warnings", "INFO",
           "Original-run warnings cannot be checked from saved results; see 06_xCell_sample_checks.csv for affected-looking samples")
}

## 8. CIBERSORT diagnostics ------------------------------------
cib_diag <- if (!is.null(method_data$CIBERSORT)) method_data$CIBERSORT$diagnostics else NULL

if (!is.null(method_data$CIBERSORT)) {
  cmat <- method_data$CIBERSORT$scores
  add_flag("CIBERSORT", "cell-type columns", "INFO",
           paste(ncol(cmat), "columns:", paste(colnames(cmat), collapse = "; ")))
  if (is.null(cib_diag)) {
    add_flag("CIBERSORT", "fit diagnostics", "REVIEW",
             "No P-value/Correlation/RMSE found in result or fit-metrics file")
  }
}

if (!is.null(cib_diag) && ncol(cib_diag) > 0) {
  diagnostic_summary <- do.call(rbind, lapply(colnames(cib_diag), function(nm) {
    v <- cib_diag[[nm]]; f <- v[is.finite(v)]
    data.frame(diagnostic = nm, n = length(v), missing = sum(is.na(v)),
               minimum = if (length(f)) round(min(f), 4) else NA,
               median  = if (length(f)) round(median(f), 4) else NA,
               mean    = if (length(f)) round(mean(f), 4) else NA,
               maximum = if (length(f)) round(max(f), 4) else NA,
               stringsAsFactors = FALSE)
  }))
  
  if ("P-value" %in% colnames(cib_diag)) {
    p <- cib_diag[["P-value"]]
    diagnostic_summary$samples_p_below_0.05 <- NA
    diagnostic_summary$samples_p_below_0.05[diagnostic_summary$diagnostic == "P-value"] <-
      sum(p < 0.05, na.rm = TRUE)
    if (all(is.na(p)) || all(p == 9999, na.rm = TRUE)) {
      add_flag("CIBERSORT", "P-values", "REVIEW",
               "P-values absent or 9999 — run likely used perm = 0, so per-sample significance can't be assessed")
    } else {
      add_flag("CIBERSORT", "P-values", "INFO",
               sprintf("%d of %d samples with P < 0.05", sum(p < 0.05, na.rm = TRUE), length(p)))
    }
  }
  if ("Correlation" %in% colnames(cib_diag)) {
    r <- cib_diag[["Correlation"]]
    add_flag("CIBERSORT", "fit correlation", "INFO",
             sprintf("median %.3f; %d samples below 0.5", median(r, na.rm = TRUE),
                     sum(r < 0.5, na.rm = TRUE)))
  }
  write.csv(diagnostic_summary, file.path(qc_dir, "03_CIBERSORT_diagnostics.csv"), row.names = FALSE)
}

## 9. Plots ----------------------------------------------------
if (length(method_data)) {
  png(file.path(qc_dir, "04_score_distributions.png"), width = 2200, height = 1700, res = 150)
  op <- par(no.readonly = TRUE)
  par(mfrow = c(3, 2), mar = c(11, 4, 3, 1))
  for (m in names(method_data)) {
    boxplot(as.data.frame(method_data[[m]]$scores), main = paste(m, "distributions"),
            ylab = "Abundance / score", las = 2, cex.axis = 0.5,
            outline = FALSE, col = "#2B5C8F")
  }
  par(op); dev.off()
  
  props <- intersect(sum_to_one_methods, names(method_data))
  if (length(props)) {
    png(file.path(qc_dir, "05_row_sums.png"), width = 1800, height = 550, res = 150)
    op <- par(no.readonly = TRUE)
    par(mfrow = c(1, length(props)), mar = c(5, 4, 3, 1))
    for (m in props) {
      hist(rowSums(method_data[[m]]$scores, na.rm = TRUE), breaks = 40,
           main = paste(m, "row sums"), xlab = "Sum of fractions per sample",
           col = "grey60", border = "white")
      abline(v = 1, col = "red", lty = 2)
    }
    par(op); dev.off()
  }
}

if (!is.null(cib_diag)) {
  plot_cols <- colnames(cib_diag)[sapply(cib_diag, function(v) {
    f <- v[is.finite(v)]; length(f) > 1 && !all(f == 9999) })]
  if (length(plot_cols)) {
    png(file.path(qc_dir, "07_CIBERSORT_diagnostics.png"),
        width = 550 * length(plot_cols), height = 550, res = 150)
    op <- par(no.readonly = TRUE)
    par(mfrow = c(1, length(plot_cols)), mar = c(5, 4, 3, 1))
    for (nm in plot_cols) {
      v <- cib_diag[[nm]]
      hist(v[is.finite(v)], breaks = 30, main = paste("CIBERSORT", nm), xlab = nm,
           col = "#E41A1C", border = "white")
    }
    par(op); dev.off()
  }
}

## 10. Aligned matrices for downstream analysis ----------------
for (m in names(method_data)) {
  out <- data.frame(sample_id = rownames(method_data[[m]]$scores),
                    method_data[[m]]$scores, check.names = FALSE)
  write.csv(out, file.path(qc_dir, paste0("aligned_", m, "_scores.csv")), row.names = FALSE)
  saveRDS(out, file.path(qc_dir, paste0("aligned_", m, "_scores.rds")))
}
if (!is.null(cib_diag)) {
  saveRDS(data.frame(sample_id = expected_samples, cib_diag, check.names = FALSE),
          file.path(qc_dir, "aligned_CIBERSORT_diagnostics.rds"))
}

## 11. Overall status + report ---------------------------------
flag_df <- if (length(flags)) do.call(rbind, flags) else
  data.frame(method = character(), check = character(),
             level = character(), detail = character())
flag_df <- flag_df[order(match(flag_df$level, c("FAIL", "REVIEW", "INFO")), flag_df$method), ]
write.csv(flag_df, file.path(qc_dir, "00_QC_flags.csv"), row.names = FALSE)

status <- sapply(names(result_files), function(m) {
  lv <- flag_df$level[flag_df$method == m]
  if ("FAIL" %in% lv) "FAIL" else if ("REVIEW" %in% lv) "REVIEW" else "PASS"
})

if (!is.null(sample_qc)) {
  sample_qc$QC_status <- status[sample_qc$method]
} else {
  sample_qc <- data.frame(method = names(result_files))
  sample_qc$QC_status <- status
}
write.csv(sample_qc, file.path(qc_dir, "01_integrated_sample_QC.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(qc_dir, "00_sessionInfo.txt"))

cat("\n======================================================================\n")
cat("INTEGRATED IMMUNE DECONVOLUTION QC\n")
cat("======================================================================\n\n")
cat("Load log:\n");     print(load_log, row.names = FALSE)
cat("\nSample QC:\n");  print(sample_qc, row.names = FALSE)
cat("\nFlags (FAIL / REVIEW / INFO):\n")
op <- options(width = 200); print(flag_df, row.names = FALSE); options(op)
cat("\nOverall status per method:\n"); print(status)
cat("\nOutputs saved to:\n", qc_dir, "\n")