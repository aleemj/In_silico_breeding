############################################################
# Import external phased haplotypes into AlphaSimR
#
# Adrian Lee
# Adapted from AlphaSimR example by
# Chris Gaynor, Jon Bancic, Philip Greenspoon,
# Gregor Gorjanc
############################################################

rm(list = ls())

library(vcfR)
library(AlphaSimR)

############################################################
# Step 1. Read phased VCF
############################################################

vcf <- read.vcfR(
  "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/BEAGLE/arugula_100k_phased.vcf.gz"
)

############################################################
# Step 2. Extract marker information
############################################################

fix <- getFIX(vcf)

############################################################
# Step 3. Extract phased genotypes
############################################################

gt <- extract.gt(vcf, element = "GT")

############################################################
# Step 3.1 Keep only the 11 pseudochromosomes
############################################################

keepChr <- paste0("scaffold_", 1:11)
keep <- fix[, "CHROM"] %in% keepChr
fix <- fix[
  keep,
  ,
  drop = FALSE
]

gt <- gt[
  keep,
  ,
  drop = FALSE
]

genMap <- data.frame(
  
  marker = fix[, "ID"],
  chromosome = sub(
    "scaffold_",
    "",
    fix[, "CHROM"]
  ),
  position = as.numeric(fix[, "POS"])/1e8,
  
  stringsAsFactors = FALSE
  
)

genMap$chromosome <- as.integer(genMap$chromosome)

############################################################
# Step 4. Confirm the VCF is phased
############################################################

if(any(grepl("/", gt))){
  stop(
    "VCF contains unphased genotypes (/). ",
    "Please phase the VCF before importing."
  )
}

############################################################
# Step 5. Split phased haplotypes
############################################################

nMarkers <- nrow(gt)
nInd <- ncol(gt)

hap1 <- matrix(
  NA_integer_,
  nrow = nMarkers,
  ncol = nInd
)

hap2 <- matrix(
  NA_integer_,
  nrow = nMarkers,
  ncol = nInd
)

for(i in seq_len(nMarkers)){
  
  tmp <- strsplit(
    gt[i, ],
    "\\|"
  )
  
  hap1[i, ] <- as.integer(
    sapply(tmp, `[`, 1)
  )
  
  hap2[i, ] <- as.integer(
    sapply(tmp, `[`, 2)
  )
  
}

############################################################
# Step 6. Build AlphaSimR haplotype matrix
############################################################

haplo <- cbind(
  hap1,
  hap2
)

haplo <- t(haplo)

colnames(haplo) <- genMap$marker

rm(
  hap1,
  hap2
)

gc()

############################################################
# Step 7. Quality control
############################################################

cat("\n")

cat("Individuals :", nInd, "\n")
cat("Haplotypes  :", nrow(haplo), "\n")
cat("Markers     :", ncol(haplo), "\n\n")


############################################################
# Step 8. Recover original sample IDs
############################################################

sampleIDs <- colnames(gt)

cat("First sample IDs:\n")

print(
  head(sampleIDs)
)

cat("\nNumber of IDs:",
    length(sampleIDs),
    "\n")

cat("Unique IDs:",
    length(unique(sampleIDs)),
    "\n\n")


############################################################
# Step 9. Create founder pedigree
############################################################

ped <- data.frame(
  id = sampleIDs,
  mother = rep(0, nInd),
  father = rep(0, nInd),
  stringsAsFactors = FALSE
)


############################################################
# Step 10. Import founder population
############################################################
founderPop <- importHaplo(
  haplo = haplo,
  genMap = genMap,
  ploidy = 2,
  ped = ped
)


############################################################
# Step 11. Check founder population
############################################################

print(founderPop)
cat("\nPedigree:\n")
print(
  head(
    getPed(founderPop)
  )
)

############################################################
# Step 12. Save founder population
############################################################

saveRDS(
  founderPop,
  "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/arugula_founderPop_100k.rds"
)

saveRDS(
  genMap,
  "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/arugula_genmap_100k.rds"
)

saveRDS(
  haplo,
  "/Users/adrianlee/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/NUS/Breeding/ISB/arugula_haplo_100k.rds"
)

