############################################################
# V2.1 - Simulate Progeny from a Single Cross
#
# AlphaSimR 2.1.0
#
# Input:
#   - phased haplotypes from original VCF import
#   - genetic map
#   - V1 rrBLUP marker effects
#   - V1 marker means
#
# Process:
#   Parent A x Parent B
#        ↓
#   AlphaSimR meiosis/recombination
#        ↓
#   simulated progeny
#        ↓
#   progeny genotype
#        ↓
#   V1 marker effects
#        ↓
#   progeny GEBVs
#
# Output:
#   - progeny GEBVs
#   - cross summary
#   - simulated progeny object
############################################################

library(AlphaSimR)
library(data.table)

# ==========================================================
# 1. SETTINGS
# ==========================================================

# ----------------------------------------------------------
# Input files
# ----------------------------------------------------------

haploFile <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/arugula_haplo_100k.rds"
genMapFile <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/arugula_genmap_100k.rds"
founderFile <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/arugula_founderPop_100k.rds"
markerEffectFile <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/Genomic Mating (Actual)/Results/adj_TFW_marker_effects.rds"
markerMeanFile <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/Genomic Mating (Actual)/Results/adj_TFW_marker_means.rds"

# ----------------------------------------------------------
# Output directory
# ----------------------------------------------------------

output_dir <- "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/Genomic Mating (Actual)/Results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ----------------------------------------------------------
# Trait
# ----------------------------------------------------------

trait <- "adj_TFW"

# ==========================================================
# 2. SETTINGS
# ==========================================================

parent1 <- "A10_17_1_1_A10_17_1_1"
parent2 <- "A10_17_3_2_A10_17_3_2"
nProgeny <- 1000
seed <- 12345

# ==========================================================
# 2. LOAD INPUT DATA
# ==========================================================

cat("\nLoading phased haplotypes...\n")
haplo <- readRDS(haploFile)
cat("Loading genetic map...\n")
genMap <- readRDS(genMapFile)


cat("Loading founder population...\n")

founderPop <- readRDS(
  founderFile
)


cat("Loading marker effects...\n")

marker_effects <- readRDS(
  markerEffectFile
)


cat("Loading marker means...\n")

marker_means <- readRDS(
  markerMeanFile
)


# ==========================================================
# 3. BASIC CHECKS
# ==========================================================

cat(
  "\nHaplotype dimensions:",
  nrow(haplo),
  "x",
  ncol(haplo),
  "\n"
)

cat(
  "Number of markers in genetic map:",
  nrow(genMap),
  "\n"
)

cat(
  "Number of marker effects:",
  length(marker_effects),
  "\n"
)

cat(
  "Number of marker means:",
  length(marker_means),
  "\n"
)

# ----------------------------------------------------------
# Check haplotype dimensions
# ----------------------------------------------------------

if (
  ncol(haplo) != nrow(genMap)
) {
  
  stop(
    "Number of haplotype markers does not match ",
    "number of markers in genMap."
  )
  
}

if (
  ncol(haplo) != length(marker_effects)
) {
  
  stop(
    "Number of haplotype markers does not match ",
    "number of marker effects."
  )
  
}


if (
  ncol(haplo) != length(marker_means)
) {
  
  stop(
    "Number of haplotype markers does not match ",
    "number of marker means."
  )
  
}


# ==========================================================
# 4. CHECK FOUNDERS
# ==========================================================

cat("\nFounder population:\n")

cat(
  "Class:",
  class(founderPop),
  "\n"
)

cat(
  "Individuals:",
  founderPop@nInd,
  "\n"
)

cat(
  "Chromosomes:",
  founderPop@nChr,
  "\n"
)

cat(
  "Ploidy:",
  founderPop@ploidy,
  "\n"
)

cat(
  "Markers:",
  founderPop@nLoci,
  "\n"
)

# ==========================================================
# 5. CONSTRUCT CHROMOSOME-SPECIFIC MAPS
# ==========================================================

cat("Creating chromosome-specific maps...\n")

chr_indices <- split(
  seq_len(nrow(genMap)),
  genMap$chromosome
)

# ----------------------------------------------------------
# Make AlphaSimR genetic map list
# ----------------------------------------------------------

genMap_list <- lapply(
  chr_indices,
  function(idx) {
    x <- genMap$position[idx]
    names(x) <- genMap$marker[idx]
    x
  }
)

# ==========================================================
# 6. CONSTRUCT CHROMOSOME-SPECIFIC HAPLOTYPES
# ==========================================================

cat(
  "Creating chromosome-specific haplotype matrices...\n"
)

haplo_list <- lapply(
  chr_indices,
  function(idx) {
    
    haplo[
      ,
      idx,
      drop = FALSE
    ]
    
  }
)


cat(
  "Number of chromosomes:",
  length(haplo_list),
  "\n"
)

# ==========================================================
# 7. CREATE MAPPOP
# ==========================================================
cat("Creating AlphaSimR MapPop...\n")

simMapPop <- newMapPop(
  genMap = genMap_list,
  haplotypes = haplo_list,
  ploidy = 2
)

cat(
  "\nMapPop class:",
  class(simMapPop),
  "\n"
)

cat(
  "Individuals:",
  simMapPop@nInd,
  "\n"
)

cat(
  "Chromosomes:",
  simMapPop@nChr,
  "\n"
)

cat(
  "Ploidy:",
  simMapPop@ploidy,
  "\n"
)


# ==========================================================
# 8. CREATE SIMULATION PARAMETERS
# ==========================================================

cat("\n")
cat("Creating SimParam...\n")


SP <- SimParam$new(
  simMapPop
)

# ----------------------------------------------------------
# Number of threads
# ----------------------------------------------------------

SP$nThreads <- parallel::detectCores(
  logical = FALSE
)

cat(
  "AlphaSimR threads:",
  SP$nThreads,
  "\n"
)


# ==========================================================
# 9. CONVERT MAPPOP TO POP
# ==========================================================

cat("\n")
cat("Converting MapPop to Pop...\n")


simPop <- newPop(
  simMapPop,
  simParam = SP
)


cat(
  "\nPopulation class:",
  class(simPop),
  "\n"
)


cat(
  "Individuals:",
  simPop@nInd,
  "\n"
)


# ==========================================================
# 10. RESTORE ORIGINAL SAMPLE IDs
# ==========================================================

cat("\n")
cat("Restoring original sample IDs...\n")


original_ids <- founderPop@id


if (
  length(original_ids) != simPop@nInd
) {
  
  stop(
    "Number of original IDs does not match ",
    "number of individuals in simPop."
  )
  
}

simPop@id <- original_ids

cat("\nFirst simulation population IDs:\n")

print(
  head(simPop@id)
)

# ==========================================================
# 11. CHECK REQUESTED PARENTS
# ==========================================================

cat("\n")
cat("Checking parents...\n")

p1 <- match(
  parent1,
  simPop@id
)

p2 <- match(
  parent2,
  simPop@id
)

if (is.na(p1)) {
  stop(
    "Parent 1 not found: ",
    parent1
  )
  
}

if (is.na(p2)) {
  stop(
    "Parent 2 not found: ",
    parent2
  )
}

cat(
  "\nParent 1:",
  parent1,
  "\n"
)

cat(
  "Index 1:",
  p1,
  "\n"
)

cat(
  "\nParent 2:",
  parent2,
  "\n"
)

cat(
  "Index 2:",
  p2,
  "\n"
)

# ==========================================================
# 12. CREATE CROSS PLAN
# ==========================================================

crossPlan <- matrix(
  c(
    p1,
    p2
  ),
  nrow = 1,
  ncol = 2
)

colnames(crossPlan) <- c(
  "Female",
  "Male"
)

cat("\nCross plan:\n")
print(
  crossPlan
)

# ==========================================================
# 13. SIMULATE PROGENY
# ==========================================================
cat("SIMULATING PROGENY\n")
set.seed(seed)
start_time <- Sys.time()
progeny <- makeCross(
  simPop,
  crossPlan = crossPlan,
  nProgeny = nProgeny,
  simParam = SP
)
end_time <- Sys.time()
cat(
  "\nSimulation completed.\n"
)
cat(
  "Time:",
  round(
    as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    ),
    2
  ),
  "seconds\n"
)
cat(
  "Progeny generated:",
  progeny@nInd,
  "\n"
)
# ==========================================================
# 14. EXTRACT PROGENY GENOTYPES
# ==========================================================

cat("\n")
cat("Extracting progeny genotypes...\n")
progeny_geno <- pullSegSiteGeno(
  progeny,
  simParam = SP
)
cat(
  "\nProgeny genotype dimensions:\n"
)
print(
  dim(progeny_geno)
)
# ==========================================================
# 15. CHECK MARKER ORDER
# ==========================================================

if (
  ncol(progeny_geno) != length(marker_effects)
) {
  
  stop(
    "Progeny genotype marker count does not match ",
    "V1 marker effects."
  )
  
}
# ==========================================================
# 16. CENTER PROGENY GENOTYPES
# ==========================================================

cat("Centering progeny genotypes using V1 means...\n")

Z_progeny <- sweep(
  progeny_geno,
  2,
  marker_means,
  FUN = "-"
)

# ==========================================================
# 17. CALCULATE PROGENY GEBVs
# ==========================================================
cat("Calculating progeny GEBVs...\n")

progeny_GEBV <- as.vector(
  Z_progeny %*% marker_effects
)

# ==========================================================
# 18. CREATE PROGENY RESULTS
# ==========================================================

progeny_results <- data.table(
  Progeny_ID = paste0(
    "Progeny_",
    seq_along(progeny_GEBV)
  ),
  Parent1 = parent1,
  Parent2 = parent2,
  GEBV = progeny_GEBV
)

# ----------------------------------------------------------
# Rank progeny
# ----------------------------------------------------------
progeny_results <- progeny_results[
  order(-GEBV)
]

# ==========================================================
# 19. CROSS SUMMARY
# ==========================================================

cross_summary <- data.table(
  
  Parent1 = parent1,
  
  Parent2 = parent2,
  
  N_Progeny = length(
    progeny_GEBV
  ),
  
  Mean_GEBV = mean(
    progeny_GEBV
  ),
  
  SD_GEBV = sd(
    progeny_GEBV
  ),
  
  Min_GEBV = min(
    progeny_GEBV
  ),
  
  Max_GEBV = max(
    progeny_GEBV
  ),
  
  GEBV_Q90 = quantile(
    progeny_GEBV,
    0.90
  ),
  
  GEBV_Q95 = quantile(
    progeny_GEBV,
    0.95
  ),
  
  GEBV_Q99 = quantile(
    progeny_GEBV,
    0.99
  )
  
)


# ==========================================================
# 20. PRINT RESULTS
# ==========================================================

cat("========================================\n")
cat("CROSS RESULTS\n")
cat("========================================\n")

print(
  cross_summary
)

# ==========================================================
# 21. SAVE PROGENY GEBVs
# ==========================================================

progeny_file <- file.path(
  output_dir,
  paste0(
    parent1,
    "_x_",
    parent2,
    "_progeny_GEBV.csv"
  )
)


fwrite(
  progeny_results,
  progeny_file
)

# ==========================================================
# 22. SAVE CROSS SUMMARY
# ==========================================================
summary_file <- file.path(
  output_dir,
  paste0(
    parent1,
    "_x_",
    parent2,
    "_cross_summary.csv"
  )
)

fwrite(
  cross_summary,
  summary_file
)

# ==========================================================
# 23. SAVE SIMULATED POPULATION
# ==========================================================

progeny_file_rds <- file.path(
  output_dir,
  paste0(
    parent1,
    "_x_",
    parent2,
    "_progeny.rds"
  )
)

saveRDS(
  progeny,
  progeny_file_rds
)