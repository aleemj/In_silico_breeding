############################################################
# V1 - Genomic Prediction using rrBLUP
#
# Input:
#   - AlphaSimR founderPop containing genotype data
#   - phenotype table
#
# Output:
#   - genomic prediction model
#   - marker effects
#   - GEBVs
#   - phenotype/GEBV table
############################################################

library(rrBLUP)
library(AlphaSimR)
library(data.table)


# ==========================================================
# 1. SETTINGS
# ==========================================================

genotypeFile <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/arugula_founderPop_100k.rds"
phenotypeFile <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/Genomic Mating (Actual)/phenotype.csv"
output_dir <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/Genomic Mating (Actual)/Results"
trait <- "adj_TFW"
n_folds <- 5
seed <- 12345

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==========================================================
# 2. LOAD DATA
# ==========================================================

cat("\nLoading genotype data...\n")
founderPop <- readRDS(genotypeFile)
cat("Extracting genotype matrix...\n")
geno <- pullSegSiteGeno(founderPop)
cat("\nGenotype dimensions:\n")
print(dim(geno))
cat("\nLoading phenotype data...\n")
pheno <- fread(phenotypeFile)
cat("\nPhenotype dimensions:\n")
print(dim(pheno))
cat("\nPhenotype columns:\n")
print(colnames(pheno))

# ==========================================================
# 3. CHECK GENOTYPE IDs
# ==========================================================
if (is.null(rownames(geno))) {
  stop(
    "Genotype matrix does not have row names. ",
    "We need individual IDs to match genotype and phenotype data."
  )
}
# ==========================================================
# 4. CHECK PHENOTYPE ID
# ==========================================================
if (!"Sample_ID" %in% colnames(pheno)) {
  stop(
    "Phenotype file does not contain a 'Sample_ID' column."
  )
}
if (!trait %in% colnames(pheno)) {
  stop(
    paste0(
      "Trait '", trait,
      "' was not found in phenotype file."
    )
  )
}

# ==========================================================
# 5. FIND COMMON INDIVIDUALS
# ==========================================================

geno_ids <- rownames(geno)
common_ids <- intersect(
  geno_ids,
  pheno$Sample_ID
)
cat(
  "\nNumber of common individuals:",
  length(common_ids),
  "\n"
)
if (length(common_ids) == 0) {
  stop(
    "No matching individuals found between genotype and phenotype data."
  )
}

# ==========================================================
# 6. ALIGN GENOTYPE AND PHENOTYPE
# ==========================================================

# Keep phenotype rows corresponding to genotype IDs
pheno <- pheno[
  match(common_ids, Sample_ID)
]

# Reorder genotype rows to exactly match phenotype
geno <- geno[
  common_ids,
  ,
  drop = FALSE
]

# Check alignment
if (!all(rownames(geno) == pheno$Sample_ID)) {
  stop("Genotype and phenotype IDs are not aligned.")
}

# ==========================================================
# 7. EXTRACT PHENOTYPE
# ==========================================================
y <- pheno[[trait]]
# IMPORTANT:
# Use Sample_ID, not ID
names(y) <- pheno$Sample_ID

# ==========================================================
# 8. REMOVE MISSING PHENOTYPES
# ==========================================================

keep <- !is.na(y)
y <- y[keep]
geno <- geno[
  names(y),
  ,
  drop = FALSE
]
cat(
  "\nIndividuals with phenotype:",
  length(y),
  "\n"
)

# ==========================================================
# 9. PHENOTYPE SUMMARY
# ==========================================================

cat("\nPhenotype summary:\n")
print(summary(y))

# ==========================================================
# 10. GENOTYPE QC
# ==========================================================

cat("\nStarting markers:", ncol(geno), "\n")
# Marker missingness
marker_missing <- colMeans(
  is.na(geno)
)
keep_markers <- marker_missing <= 0.10
geno <- geno[
  ,
  keep_markers,
  drop = FALSE
]
cat(
  "Markers after missingness filtering:",
  ncol(geno),
  "\n"
)

# ==========================================================
# 11. IMPUTE MISSING GENOTYPES
# ==========================================================
cat("\nImputing missing genotypes...\n")
for (j in seq_len(ncol(geno))) {
  
  missing <- is.na(geno[, j])
  
  if (any(missing)) {
    
    geno[missing, j] <- mean(
      geno[, j],
      na.rm = TRUE
    )
  }
}

# Check for remaining missing values
if (anyNA(geno)) {
  stop("Missing genotype values remain after imputation.")
}

# ==========================================================
# 12. CENTER MARKERS
# ==========================================================

marker_means <- colMeans(geno)

Z <- sweep(
  geno,
  2,
  marker_means,
  FUN = "-"
)

# ==========================================================
# 13. RRBLUP
# ==========================================================

cat("\nRunning rrBLUP...\n")
fit <- mixed.solve(
  y = y,
  Z = Z
)

# ==========================================================
# 14. EXTRACT MARKER EFFECTS
# ==========================================================

marker_effects <- fit$u

cat(
  "\nNumber of marker effects:",
  length(marker_effects),
  "\n"
)


# ==========================================================
# 15. CALCULATE GEBVs
# ==========================================================

cat("\nCalculating GEBVs...\n")

GEBV <- as.vector(
  Z %*% marker_effects
)

names(GEBV) <- rownames(Z)


# ==========================================================
# 16. CREATE RESULTS TABLE
# ==========================================================

results <- data.table(
  Sample_ID = names(GEBV),
  Phenotype = y[names(GEBV)],
  GEBV = GEBV
)


results <- results[
  order(-GEBV)
]


# ==========================================================
# 17. SAVE GEBVs
# ==========================================================

gebv_file <- file.path(
  output_dir,
  paste0(
    trait,
    "_GEBV.csv"
  )
)

fwrite(
  results,
  gebv_file
)

cat(
  "\nGEBV results saved to:\n",
  gebv_file,
  "\n"
)


# ==========================================================
# 18. SAVE MODEL
# ==========================================================

model_file <- file.path(
  output_dir,
  paste0(
    trait,
    "_rrBLUP_model.rds"
  )
)

saveRDS(
  fit,
  model_file
)

# ==========================================================
# 19. SAVE MARKER EFFECTS
# ==========================================================

marker_file <- file.path(
  output_dir,
  paste0(
    trait,
    "_marker_effects.rds"
  )
)

saveRDS(
  marker_effects,
  marker_file
)

saveRDS(
  marker_means,
  file.path(
    output_dir,
    paste0(
      trait,
      "_marker_means.rds"
    )
  )
)

# ==========================================================
# 20. SUMMARY
# ==========================================================

cat("\n========================================\n")
cat("V1 GENOMIC PREDICTION COMPLETE\n")
cat("========================================\n")

cat(
  "Trait:",
  trait,
  "\n"
)

cat(
  "Individuals:",
  nrow(geno),
  "\n"
)

cat(
  "Markers:",
  ncol(geno),
  "\n"
)

cat(
  "Phenotypic mean:",
  mean(y),
  "\n"
)

cat(
  "Phenotypic SD:",
  sd(y),
  "\n"
)

cat("\nTop 10 individuals by GEBV:\n")

print(
  head(results, 10)
)

# ==========================================================
# 21. 5-FOLD CROSS-VALIDATION
# ==========================================================

set.seed(seed)

n <- length(y)

folds <- sample(
  rep(1:n_folds, length.out = n)
)

cv_predictions <- rep(
  NA_real_,
  n
)

names(cv_predictions) <- names(y)

cv_results <- data.table(
  Fold = integer(),
  Accuracy = numeric(),
  RMSE = numeric(),
  Bias = numeric()
)


for (fold in 1:n_folds) {
  
  cat(
    "\nRunning CV fold",
    fold,
    "of",
    n_folds,
    "\n"
  )
  
  train <- folds != fold
  test  <- folds == fold
  
  Z_train <- Z[train, , drop = FALSE]
  Z_test  <- Z[test, , drop = FALSE]
  
  y_train <- y[train]
  y_test  <- y[test]
  
  # Fit model
  cv_fit <- mixed.solve(
    y = y_train,
    Z = Z_train
  )
  
  # Predict validation set
  pred <- as.vector(
    Z_test %*% cv_fit$u
  )
  
  cv_predictions[test] <- pred
  
  # Fold statistics
  fold_accuracy <- cor(
    y_test,
    pred,
    use = "complete.obs"
  )
  
  fold_rmse <- sqrt(
    mean(
      (y_test - pred)^2,
      na.rm = TRUE
    )
  )
  
  fold_bias <- mean(
    pred - y_test,
    na.rm = TRUE
  )
  
  cv_results <- rbind(
    cv_results,
    data.table(
      Fold = fold,
      Accuracy = fold_accuracy,
      RMSE = fold_rmse,
      Bias = fold_bias
    )
  )
}


# ==========================================================
# 22. OVERALL CV ACCURACY
# ==========================================================

overall_accuracy <- cor(
  y,
  cv_predictions,
  use = "complete.obs"
)

overall_RMSE <- sqrt(
  mean(
    (y - cv_predictions)^2,
    na.rm = TRUE
  )
)

overall_bias <- mean(
  cv_predictions - y,
  na.rm = TRUE
)


# ==========================================================
# 23. PRINT RESULTS
# ==========================================================

cat("\n========================================\n")
cat("5-FOLD CROSS-VALIDATION RESULTS\n")
cat("========================================\n\n")

print(cv_results)

cat("\nOverall prediction accuracy:",
    round(overall_accuracy, 4),
    "\n")

cat("Overall RMSE:",
    round(overall_RMSE, 4),
    "\n")

cat("Overall bias:",
    round(overall_bias, 4),
    "\n")


# ==========================================================
# 24. SAVE CV RESULTS
# ==========================================================

fwrite(
  cv_results,
  file.path(
    output_dir,
    paste0(
      trait,
      "_CV_fold_results.csv"
    )
  )
)

cv_predictions_table <- data.table(
  Sample_ID = names(cv_predictions),
  Observed = y[names(cv_predictions)],
  Predicted = cv_predictions,
  Fold = folds
)

fwrite(
  cv_predictions_table,
  file.path(
    output_dir,
    paste0(
      trait,
      "_CV_predictions.csv"
    )
  )
)