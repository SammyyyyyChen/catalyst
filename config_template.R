# ---- Local configuration ------------------------------------------------
# Copy this file to config.R and fill in your paths.
# config.R is gitignored so your local paths stay private.

# Absolute or relative path to the folder that contains "raw files test/"
data_dir <- "path/to/Spheroid/Results/ExpSpheroid011_temporary/high dimensional analysis"

# Subdirectory inside data_dir that holds the .fcs files
fcs_subdir <- "raw files test"

# ---- Preserving your tSNE & clustering ----------------------------------
# If you already have an sce_IgE_clean.rds from a previous run (with the
# tSNE orientation and cluster assignments you like), copy it to:
#   output/rds/sce_IgE_clean.rds
# Then set RECOMPUTE <- FALSE in analysis/01_IgE_analysis.R (the default).
# The script will load that file and skip sections 1-5.
