# ---- Run both DA and DS for one contrast ---------------------------------

run_diffcyt_pair <- function(sce, design, contrast, clustering = "cluster_ID") {
  da <- diffcyt(sce,
                design         = design,
                contrast       = contrast,
                analysis_type  = "DA",
                method_DA      = "diffcyt-DA-edgeR",
                clustering_to_use = clustering)
  ds <- diffcyt(sce,
                design         = design,
                contrast       = contrast,
                analysis_type  = "DS",
                method_DS      = "diffcyt-DS-limma",
                clustering_to_use = clustering)
  list(da = da, ds = ds)
}

# ---- Format adjusted p-values for plot labels ----------------------------

prep_da_stats <- function(da_obj, label_name = "p") {
  res <- as.data.frame(rowData(da_obj$res))
  if (!"cluster_id" %in% colnames(res))
    res <- res %>% rownames_to_column("cluster_id")
  res %>%
    select(cluster_id, p_adj) %>%
    mutate(
      !!paste0("p_label_", label_name) :=
        ifelse(p_adj < 0.0001, "<0.0001",
               ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))
    ) %>%
    select(-p_adj)
}

prep_ds_stats <- function(ds_obj, label_name = "p") {
  res <- as.data.frame(rowData(ds_obj$res))
  res %>%
    select(cluster_id, marker_id, p_adj) %>%
    mutate(
      !!paste0("p_label_", label_name) :=
        ifelse(p_adj < 0.0001, "<0.0001",
               ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))
    ) %>%
    select(-p_adj)
}

# ---- Antibody comparisons (NIP vs IgE, PBS vs IgE) ----------------------
# Returns a named list with NIP_vs_IgE and PBS_vs_IgE, each containing da/ds.

run_antibody_comparisons <- function(sce, cols_design, clustering = "cluster_ID") {
  design <- createDesignMatrix(metadata(sce)$experiment_info,
                               cols_design = cols_design)
  n <- ncol(design)

  # NIP IgE vs CSPG4 IgE: coefficient for CSPG4 minus NIP
  # In the design with (Intercept, condNIP_IgE, condCSPG4_IgE, ...),
  # contrast = c(0, -1, 1, 0, ...) for NIP vs IgE
  contrast_nip  <- numeric(n)
  contrast_nip[2] <- -1
  contrast_nip[3] <-  1

  # PBS vs CSPG4 IgE: just the CSPG4 coefficient
  contrast_pbs  <- numeric(n)
  contrast_pbs[3] <- 1

  list(
    NIP_vs_IgE = run_diffcyt_pair(sce, design,
                                  createContrast(contrast_nip), clustering),
    PBS_vs_IgE = run_diffcyt_pair(sce, design,
                                  createContrast(contrast_pbs), clustering),
    design = design
  )
}

# ---- HV vs M comparison within a single condition subset -----------------

run_hv_vs_m <- function(sce, cond_value, cols_design, clustering = "cluster_ID") {
  sce_sub <- filterSCE(sce, condition == cond_value)
  design  <- createDesignMatrix(metadata(sce_sub)$experiment_info,
                                cols_design = cols_design)
  n <- ncol(design)
  contrast_vec <- numeric(n)
  contrast_vec[n] <- 1
  run_diffcyt_pair(sce_sub, design, createContrast(contrast_vec), clustering)
}

# ---- 2D vs 3D comparison within a single condition -----------------------

run_2d_vs_3d <- function(sce, cond_value, clustering = "cluster_ID") {
  sce_sub <- filterSCE(sce, condition == cond_value)
  design  <- createDesignMatrix(metadata(sce_sub)$experiment_info,
                                cols_design = c("setting", "sample_type"))
  contrast_vec <- c(0, 1, 0)
  run_diffcyt_pair(sce_sub, design, createContrast(contrast_vec), clustering)
}
