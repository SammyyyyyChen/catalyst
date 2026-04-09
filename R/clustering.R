# ---- Merging tables -------------------------------------------------------
# Hand-curated from expression heatmaps. Each row maps an original FlowSOM
# meta-cluster to a merged cluster ID or phenotype label.

# Merging 1: for sce_IgE (meta12 -> 11 clusters)
merging_table_1 <- data.frame(
  cluster = as.character(1:12),
  merged_cluster = factor(
    c("2", "10", "7", "8", "9", "8", "3", "11", "4", "5", "6", "1"),
    levels = as.character(1:11)
  )
)

# Merging 2: for sce_IgE_clean (meta16 -> 14 clusters)
merging_table_2 <- data.frame(
  cluster = as.character(1:16),
  merged_cluster = factor(
    c("10", "11", "12", "8", "12", "2", "5", "9", "13",
      "1", "6", "14", "4", "7", "3", "3"),
    levels = as.character(1:14)
  )
)

# Merging 3: phenotype labels for sce_IgE_clean (cluster_ID -> Phenotype)
merging_table_3 <- data.frame(
  cluster = as.character(1:14),
  merged_cluster = c(
    "Classical, resting",
    "Classical, resting",
    "Non-classical, resting",
    "Classical, activated/antigen-presenting",
    "Intermediate, activated/antigen-presenting",
    "Non-classical, activated/antigen-presenting",
    "Non-classical, activated/antigen-presenting",
    "Classical, exhausted",
    "Classical, inhibited",
    "Classical, exhausted",
    "Classical, inhibited",
    "Classical, exhausted",
    "Classical, inhibited",
    "Classical, exhausted"
  )
)

# Global merging: collapse all clusters into one for global DS tests
merging_table_global <- data.frame(
  cluster        = as.character(1:16),
  merged_cluster = rep("1", 16)
)

# ---- Clustering helpers --------------------------------------------------

run_clustering <- function(sce, max_k = 16, xdim = 10, ydim = 10, seed = 123) {
  cluster(sce, features = "state",
          xdim = xdim, ydim = ydim, maxK = max_k, seed = seed)
}

apply_merges <- function(sce) {
  sce <- mergeClusters(sce, k = "meta16",
                       table = merging_table_2,
                       id = "cluster_ID", overwrite = TRUE)
  sce <- mergeClusters(sce, k = "cluster_ID",
                       table = merging_table_3,
                       id = "Phenotype", overwrite = TRUE)
  sce <- mergeClusters(sce, k = "meta16",
                       table = merging_table_global,
                       id = "global", overwrite = TRUE)
  sce
}
