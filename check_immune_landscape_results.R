project_dir <- "C:/Users/karan/OneDrive/Documents/Breast_Cancer_Immune_Infiltration"
out_dir <- file.path(project_dir, "results", "immune_landscape")

# 1. Histology by cluster
hist_tab <- read.csv(file.path(out_dir, "cluster_vs_Histology.csv"), check.names = FALSE)
cat("===== Histology by cluster =====\n")
print(hist_tab, row.names = FALSE)

# 2. Stage by cluster
stage_tab <- read.csv(file.path(out_dir, "cluster_vs_Stage.csv"), check.names = FALSE)
cat("\n===== Stage by cluster =====\n")
print(stage_tab, row.names = FALSE)

# 3. Test p-values for reference
cat("\n===== Association tests =====\n")
print(read.csv(file.path(out_dir, "cluster_clinical_association_tests.csv")), row.names = FALSE)

# 4. Immune score: one hump or two?
res <- readRDS(file.path(out_dir, "sample_level_results.rds"))
s <- res$immune_score

suppressPackageStartupMessages(library(mclust))
bic <- mclustBIC(s, G = 1:3, modelNames = "V", verbose = FALSE)
cat("\n===== Mixture model fit (higher BIC = better) =====\n")
print(round(bic[, "V"], 1))
cat("BIC gain, 2 groups vs 1:", round(bic[2, "V"] - bic[1, "V"], 1),
    "  (>10 = strong evidence for 2 groups)\n")

mc2 <- Mclust(s, G = 2, modelNames = "V", verbose = FALSE)
cat("\nTwo-group fit: means =", round(mc2$parameters$mean, 2),
    "| proportions =", round(mc2$parameters$pro, 2), "\n")

if (!requireNamespace("diptest", quietly = TRUE)) install.packages("diptest")
dt <- diptest::dip.test(s)
cat("\nDip test for multimodality: D =", round(dt$statistic, 4),
    "| p =", signif(dt$p.value, 3),
    "\n(p < 0.05 = evidence of more than one hump; p >= 0.05 = consistent with one hump)\n")

cat("\nImmune score summary:\n"); print(summary(s))