############################################################
# V2.2 - All Pairwise Genomic Cross Simulation
#
# AlphaSimR 2.1.0
#
# Purpose:
#   Simulate progeny for every possible pair of founders
#   and calculate the distribution of progeny GEBVs.
#
# Workflow:
#
#   Founder population
#          ↓
#   all pairwise crosses
#          ↓
#   simulate nProgeny offspring
#          ↓
#   calculate progeny GEBVs
#          ↓
#   summarize cross
#          ↓
#   discard progeny
#          ↓
#   next cross
#
# IMPORTANT:
#   This version performs a SCREENING simulation.
#
#   Recommended:
#       nProgeny = 20–50
#
#   Then perform a larger simulation on the best crosses.
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
# Output
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


# ----------------------------------------------------------
# Number of progeny per cross
#
# SCREENING VALUE
#
# Start with 20.
#
# Later:
#   50
#   100
#
# for higher precision.
# ----------------------------------------------------------

nProgeny <- 20


# ----------------------------------------------------------
# Random seed
# ----------------------------------------------------------

seed <- 12345


# ----------------------------------------------------------
# Save progress every N crosses
# ----------------------------------------------------------

checkpoint_every <- 1000


# ----------------------------------------------------------
# Starting cross
#
# Allows restarting after interruption.
#
# Example:
#
# start_cross <- 1
#
# or resume from:
#
# start_cross <- 5001
# ----------------------------------------------------------

start_cross <- 1


# ==========================================================
# 2. LOAD DATA
# ==========================================================

cat("\n")
cat("========================================\n")
cat("LOADING DATA\n")
cat("========================================\n")


cat("\nLoading haplotypes...\n")

haplo <- readRDS(
  haploFile
)


cat("Loading genetic map...\n")

genMap <- readRDS(
  genMapFile
)


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

cat("\n")
cat("========================================\n")
cat("INPUT CHECKS\n")
cat("========================================\n")


cat(
  "Individuals:",
  founderPop@nInd,
  "\n"
)


cat(
  "Markers:",
  ncol(haplo),
  "\n"
)


cat(
  "Marker effects:",
  length(marker_effects),
  "\n"
)


cat(
  "Marker means:",
  length(marker_means),
  "\n"
)


if (
  ncol(haplo) != nrow(genMap)
) {
  
  stop(
    "Haplotype markers and genetic map markers do not match."
  )
}


if (
  ncol(haplo) != length(marker_effects)
) {
  
  stop(
    "Haplotype markers and marker effects do not match."
  )
}


if (
  ncol(haplo) != length(marker_means)
) {
  
  stop(
    "Haplotype markers and marker means do not match."
  )
}


# ==========================================================
# 4. CONSTRUCT CHROMOSOME-SPECIFIC DATA
# ==========================================================

cat("\n")
cat("Creating chromosome-specific data...\n")


chr_indices <- split(
  seq_len(nrow(genMap)),
  genMap$chromosome
)


# ----------------------------------------------------------
# Genetic map
# ----------------------------------------------------------

genMap_list <- lapply(
  chr_indices,
  function(idx) {
    
    x <- genMap$position[idx]
    
    names(x) <- genMap$marker[idx]
    
    x
    
  }
)


# ----------------------------------------------------------
# Haplotypes
# ----------------------------------------------------------

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
  "Chromosomes:",
  length(haplo_list),
  "\n"
)


# ==========================================================
# 5. CREATE MAPPOP
# ==========================================================

cat("\n")
cat("Creating MapPop...\n")


simMapPop <- newMapPop(
  
  genMap = genMap_list,
  
  haplotypes = haplo_list,
  
  ploidy = 2
  
)


# ==========================================================
# 6. CREATE SIMULATION PARAMETERS
# ==========================================================

cat("\n")
cat("Creating SimParam...\n")


SP <- SimParam$new(
  simMapPop
)


# ----------------------------------------------------------
# Use available physical CPU cores
# ----------------------------------------------------------

nThreads <- parallel::detectCores(
  logical = FALSE
)


SP$nThreads <- nThreads


cat(
  "Threads:",
  SP$nThreads,
  "\n"
)


# ==========================================================
# 7. CONVERT TO POP
# ==========================================================

cat("\n")
cat("Creating simulation population...\n")


simPop <- newPop(
  simMapPop,
  simParam = SP
)


# ----------------------------------------------------------
# Restore original IDs
# ----------------------------------------------------------

simPop@id <- founderPop@id


cat(
  "Population size:",
  simPop@nInd,
  "\n"
)


# ==========================================================
# 8. CREATE ALL UNIQUE CROSS PAIRS
# ==========================================================

cat("\n")
cat("Creating cross list...\n")


nParents <- simPop@nInd


cross_pairs <- t(
  combn(
    nParents,
    2
  )
)


cross_pairs <- data.table(
  
  Cross_ID = seq_len(
    nrow(cross_pairs)
  ),
  
  Parent1_Index = cross_pairs[, 1],
  
  Parent2_Index = cross_pairs[, 2],
  
  Parent1 = simPop@id[
    cross_pairs[, 1]
  ],
  
  Parent2 = simPop@id[
    cross_pairs[, 2]
  ]
  
)


cat(
  "Number of crosses:",
  nrow(cross_pairs),
  "\n"
)


# ==========================================================
# 9. OUTPUT FILE
# ==========================================================

results_file <- file.path(
  
  output_dir,
  
  paste0(
    trait,
    "_all_crosses_",
    nProgeny,
    "progeny.csv"
  )
  
)


# ==========================================================
# 10. INITIALIZE RESULTS
# ==========================================================

results <- data.table(
  
  Cross_ID = integer(),
  
  Parent1 = character(),
  
  Parent2 = character(),
  
  N_Progeny = integer(),
  
  Mean_GEBV = numeric(),
  
  SD_GEBV = numeric(),
  
  Min_GEBV = numeric(),
  
  Max_GEBV = numeric(),
  
  GEBV_Q90 = numeric(),
  
  GEBV_Q95 = numeric(),
  
  GEBV_Q99 = numeric()
  
)


# ==========================================================
# 11. RUN ALL CROSSES
# ==========================================================

cat("\n")
cat("========================================\n")
cat("STARTING CROSS SIMULATION\n")
cat("========================================\n")


set.seed(seed)


total_crosses <- nrow(cross_pairs)


overall_start <- Sys.time()


for (
  i in start_cross:total_crosses
) {
  
  
  # --------------------------------------------------------
  # Current cross
  # --------------------------------------------------------
  
  p1 <- cross_pairs$Parent1_Index[i]
  
  p2 <- cross_pairs$Parent2_Index[i]
  
  
  parent1 <- cross_pairs$Parent1[i]
  
  parent2 <- cross_pairs$Parent2[i]
  
  
  # --------------------------------------------------------
  # Progress
  # --------------------------------------------------------
  
  if (
    i == start_cross ||
    i %% checkpoint_every == 0
  ) {
    
    elapsed <- as.numeric(
      difftime(
        Sys.time(),
        overall_start,
        units = "mins"
      )
    )
    
    completed <- i - start_cross + 1
    
    rate <- completed / max(
      elapsed,
      1e-6
    )
    
    remaining <- total_crosses - i
    
    eta <- remaining / max(
      rate,
      1e-6
    )
    
    
    cat(
      "\nCross:",
      i,
      "/",
      total_crosses,
      "\n"
    )
    
    cat(
      "Parent 1:",
      parent1,
      "\n"
    )
    
    cat(
      "Parent 2:",
      parent2,
      "\n"
    )
    
    cat(
      "Rate:",
      round(rate, 2),
      "crosses/min\n"
    )
    
    cat(
      "Estimated remaining time:",
      round(eta / 60, 2),
      "hours\n"
    )
    
  }
  
  
  # --------------------------------------------------------
  # Cross plan
  # --------------------------------------------------------
  
  crossPlan <- matrix(
    
    c(
      p1,
      p2
    ),
    
    nrow = 1,
    
    ncol = 2
    
  )
  
  
  # --------------------------------------------------------
  # Simulate progeny
  # --------------------------------------------------------
  
  progeny <- makeCross(
    
    simPop,
    
    crossPlan = crossPlan,
    
    nProgeny = nProgeny,
    
    simParam = SP
    
  )
  
  
  # --------------------------------------------------------
  # Extract progeny genotypes
  # --------------------------------------------------------
  
  progeny_geno <- pullSegSiteGeno(
    
    progeny,
    
    simParam = SP
    
  )
  
  
  # --------------------------------------------------------
  # Apply V1 marker centering
  # --------------------------------------------------------
  
  Z_progeny <- sweep(
    
    progeny_geno,
    
    2,
    
    marker_means,
    
    FUN = "-"
    
  )
  
  
  # --------------------------------------------------------
  # Calculate progeny GEBVs
  # --------------------------------------------------------
  
  progeny_GEBV <- as.vector(
    
    Z_progeny %*% marker_effects
    
  )
  
  
  # --------------------------------------------------------
  # Summarize cross
  # --------------------------------------------------------
  
  cross_result <- data.table(
    
    Cross_ID = i,
    
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
    
    GEBV_Q90 = as.numeric(
      quantile(
        progeny_GEBV,
        0.90
      )
    ),
    
    GEBV_Q95 = as.numeric(
      quantile(
        progeny_GEBV,
        0.95
      )
    ),
    
    GEBV_Q99 = as.numeric(
      quantile(
        progeny_GEBV,
        0.99
      )
    )
    
  )
  
  
  # --------------------------------------------------------
  # Add result
  # --------------------------------------------------------
  
  results <- rbind(
    
    results,
    
    cross_result
    
  )
  
  
  # --------------------------------------------------------
  # Delete progeny
  # --------------------------------------------------------
  
  rm(
    progeny,
    progeny_geno,
    Z_progeny,
    progeny_GEBV
  )
  
  
  # --------------------------------------------------------
  # Garbage collection
  # --------------------------------------------------------
  
  if (
    i %% checkpoint_every == 0
  ) {
    
    gc()
    
    
    # ------------------------------------------------------
    # Save checkpoint
    # ------------------------------------------------------
    
    fwrite(
      
      results,
      
      results_file
      
    )
    
    cat(
      "Checkpoint saved.\n"
    )
    
  }
  
}


# ==========================================================
# 12. FINAL SAVE
# ==========================================================

fwrite(
  
  results,
  
  results_file
  
)


# ==========================================================
# 13. SORT RESULTS
# ==========================================================

results_ranked <- results[
  order(
    -Mean_GEBV
  )
]


ranked_file <- file.path(
  
  output_dir,
  
  paste0(
    trait,
    "_all_crosses_",
    nProgeny,
    "progeny_RANKED.csv"
  )
  
)


fwrite(
  
  results_ranked,
  
  ranked_file
  
)


# ==========================================================
# 14. FINAL SUMMARY
# ==========================================================

overall_time <- as.numeric(
  
  difftime(
    Sys.time(),
    overall_start,
    units = "hours"
  )
  
)


cat("\n")
cat("========================================\n")
cat("ALL-CROSS SIMULATION COMPLETE\n")
cat("========================================\n")


cat(
  "\nNumber of crosses:",
  nrow(results),
  "\n"
)


cat(
  "Progeny per cross:",
  nProgeny,
  "\n"
)


cat(
  "Total simulated progeny:",
  nrow(results) * nProgeny,
  "\n"
)


cat(
  "Total time:",
  round(
    overall_time,
    2
  ),
  "hours\n"
)


cat("\nTop 20 crosses:\n")


print(
  head(
    results_ranked,
    20
  )
)


cat(
  "\nResults saved to:\n",
  ranked_file,
  "\n"
)