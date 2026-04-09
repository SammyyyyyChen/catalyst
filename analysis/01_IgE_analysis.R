# ============================================================================
# CATALYST IgE Flow Cytometry Analysis
# Monocyte phenotyping in melanoma spheroid co-culture with anti-CSPG4 IgE
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------
source("config.R")
source("R/packages.R")
source("R/theme.R")
source("R/data_loading.R")
source("R/clustering.R")
source("R/differential.R")
source("R/plotting.R")

fig_dir <- file.path("output", "figures")
rds_dir <- file.path("output", "rds")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# Set to TRUE to recompute clustering & tSNE from scratch.
# Set to FALSE (default) to load the saved checkpoint if it exists, which
# preserves your existing tSNE orientation and cluster assignments.
RECOMPUTE <- FALSE

checkpoint_file <- file.path(rds_dir, "sce_IgE_clean.rds")

if (!RECOMPUTE && file.exists(checkpoint_file)) {

  message("Loading cached checkpoint: ", checkpoint_file)
  sce_IgE_clean <- readRDS(checkpoint_file)

} else {

  message("Computing from scratch (RECOMPUTE=TRUE or no checkpoint found)...")

  # ---- 1. Load FCS data & build SCE ---------------------------------------
  fcs_path  <- file.path(data_dir, fcs_subdir)
  fileNames <- list_fcs_files(fcs_path)
  validate_fcs_files(fileNames)

  metadata <- build_metadata(fileNames)
  panel    <- build_panel(fileNames)
  write.csv(metadata, file.path(rds_dir, "metadata.csv"), row.names = FALSE)
  write.csv(panel,    file.path(rds_dir, "panel.csv"),    row.names = FALSE)

  sce <- load_and_prep_sce(fileNames, panel, metadata)
  sce_IgE <- prep_ige_sce(sce)

  table(colData(sce_IgE)$patient_id)
  table(colData(sce_IgE)$condition)
  table(colData(sce_IgE)$setting)

  # ---- 2. QC plots (pre-clean) --------------------------------------------
  save_figure("IgE_cell_counts_by_condition.png",
              plot_cell_counts(sce_IgE), width = 19, height = 9)

  save_figure("MDS_IgE_by_condition_setting.png",
              pbMDS(sce_IgE, color_by = "condition", shape_by = "setting",
                    label_by = NULL, size_by = TRUE) +
                scale_color_manual(values = condition_colors) + theme_pub())

  save_figure("MDS_IgE_by_batch.png",
              pbMDS(sce_IgE, color_by = "batch", label_by = NULL,
                    shape_by = "condition", size_by = TRUE) + theme_pub())

  save_figure("NRS_IgE_by_condition.png",
              plotNRS(sce_IgE, features = "state", color_by = "condition") +
                scale_color_manual(values = condition_colors) + theme_pub(),
              width = 12, height = 6)

  # ---- 3. Remove problematic batch & repeat QC ----------------------------
  batches_keep <- c("20251117", "20251205", "20251212", "20250131", "20250206")
  sce_IgE_clean <- filterSCE(sce_IgE, batch %in% batches_keep)

  save_figure("IgE_clean_cell_counts.png",
              plot_cell_counts(sce_IgE_clean), width = 19, height = 9)

  save_figure("MDS_IgE_clean_by_condition_setting.png",
              pbMDS(sce_IgE_clean, color_by = "condition", shape_by = "setting",
                    label_by = NULL, size_by = TRUE) +
                scale_color_manual(values = condition_colors) + theme_pub())

  save_figure("MDS_IgE_clean_by_batch.png",
              pbMDS(sce_IgE_clean, color_by = "batch", label_by = NULL,
                    shape_by = "condition", size_by = TRUE) + theme_pub())

  save_figure("NRS_IgE_clean_by_condition.png",
              plotNRS(sce_IgE_clean, features = "state", color_by = "condition") +
                scale_color_manual(values = condition_colors) + theme_pub(),
              width = 12, height = 6)

  # ---- 4. Clustering ------------------------------------------------------
  rowData(sce_IgE_clean)$marker_class <- as.character(rowData(sce_IgE_clean)$marker_class)
  rowData(sce_IgE_clean)$marker_class[rownames(sce_IgE_clean) == "CD32B"] <- "state"
  rowData(sce_IgE_clean)$marker_class <- factor(rowData(sce_IgE_clean)$marker_class)

  sce_IgE_clean <- run_clustering(sce_IgE_clean, max_k = 16)

  save_figure("IgE_clean_delta_plot.png",
              delta_area(sce_IgE_clean) + theme_pub(),
              width = 12, height = 9)

  sce_IgE_clean <- apply_merges(sce_IgE_clean)

  save_figure("IgE_clean_heatmap_cluster_ID.png",
              plotExprHeatmap(sce_IgE_clean, features = ordered_markers,
                              k = "cluster_ID", by = "cluster_id", fun = "median",
                              row_anno = TRUE, scale = "last", q = 0, bars = TRUE,
                              row_clust = FALSE, col_clust = FALSE, perc = TRUE,
                              k_pal = cluster_colors),
              width = 12, height = 7)

  # ---- 5. Dimensionality reduction ---------------------------------------
  sce_IgE_clean <- runDR(sce_IgE_clean, cell = 1000, dr = "TSNE",
                          features = "state", perplexity = 25, seed = 10)
  sce_IgE_clean <- runDR(sce_IgE_clean, dr = "UMAP", features = "state")

  # Save checkpoint with clustering + tSNE + UMAP baked in
  saveRDS(sce_IgE_clean, checkpoint_file)
  message("Saved checkpoint: ", checkpoint_file)

} # end if/else RECOMPUTE

# ---- 5b. tSNE visualisation (uses cached embedding) ----------------------
save_figure("tSNE_IgE_clean_all.png",
            plotDR(sce_IgE_clean, dr = "TSNE", color_by = "cluster_ID") +
              scale_color_manual(values = cluster_colors) + theme_pub())

# tSNE per condition
conditions <- unique(colData(sce_IgE_clean)$condition)
for (cond in conditions) {
  sce_sub <- filterSCE(sce_IgE_clean, k = "cluster_ID", condition == cond)
  save_figure(paste0("tSNE_", gsub(" ", "_", cond), ".png"),
              plotDR(sce_sub, dr = "TSNE", color_by = "cluster_ID") +
                scale_color_manual(values = cluster_colors) +
                ggtitle(cond) + theme_pub())
}

# tSNE by setting and sample type
save_figure("tSNE_IgE_clean_by_setting.png",
            plotDR(sce_IgE_clean, dr = "TSNE", color_by = "setting") +
              scale_colour_manual(values = condition_colors) + theme_pub())

save_figure("tSNE_IgE_clean_by_sample_type.png",
            plotDR(sce_IgE_clean, dr = "TSNE", color_by = "sample_type") +
              scale_colour_manual(values = condition_colors) + theme_pub())

# ---- 6. Abundance plots (descriptive) ------------------------------------
for (group_var in c("condition", "setting", "sample_type")) {
  save_figure(paste0("cluster_abundance_by_", group_var, ".png"),
              plotAbundances(sce_IgE_clean, k = "cluster_ID",
                             by = "sample_id", group_by = group_var) +
                scale_fill_manual(values = cluster_colors) +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                      text = element_text(size = 24),
                      plot.title = element_text(face = "bold", size = 28)),
              width = 10, height = 7)
}

# ---- 7. Differential analysis (combined, all settings) --------------------
# Reorder sample_id for heatmaps
new_order <- colData(sce_IgE_clean) %>% as.data.frame() %>%
  arrange(condition, sample_id) %>% pull(sample_id) %>% unique() %>% as.character()
sce_IgE_clean$sample_id <- factor(sce_IgE_clean$sample_id, levels = new_order)

diff_combined <- run_antibody_comparisons(
  sce_IgE_clean, cols_design = c("condition", "setting", "sample_type"))

message("DA: NIP vs IgE")
topTable(diff_combined$NIP_vs_IgE$da, format_vals = TRUE)
message("DA: PBS vs IgE")
topTable(diff_combined$PBS_vs_IgE$da, format_vals = TRUE)
message("DS: NIP vs IgE")
topTable(diff_combined$NIP_vs_IgE$ds, format_vals = TRUE)
message("DS: PBS vs IgE")
topTable(diff_combined$PBS_vs_IgE$ds, format_vals = TRUE)

# Abundance with p-values (exhausted clusters 11-14)
save_figure("abundance_exhausted_IgE_padj.png",
            plot_abundance_condition_pval(
              sce_IgE_clean,
              da_nip = diff_combined$NIP_vs_IgE$da,
              da_pbs = diff_combined$PBS_vs_IgE$da,
              clusters_to_show = c(11, 12, 13, 14)),
            width = 13, height = 4)

# Per-cluster marker expression loop
message("Generating per-cluster marker plots (combined)...")
plot_cluster_markers_loop(
  sce_IgE_clean,
  ds_nip  = diff_combined$NIP_vs_IgE$ds,
  ds_pbs  = diff_combined$PBS_vs_IgE$ds,
  markers = markers_ex_CD14,
  output_subdir = "pbExprs_cluster_IgE_clean")

# Global marker expression (using global clustering)
ds_global_pbs <- diffcyt(sce_IgE_clean, design = diff_combined$design,
                          contrast = createContrast(c(0, 0, 1, 0, 0)),
                          analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                          clustering_to_use = "global")
ds_global_nip <- diffcyt(sce_IgE_clean, design = diff_combined$design,
                          contrast = createContrast(c(0, -1, 1, 0, 0)),
                          analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                          clustering_to_use = "global")

save_figure("pbExprs_IgE_Inhibitory_padj.png",
            plot_global_markers_pval(sce_IgE_clean, ds_global_nip,
                                     ds_global_pbs, markers_inhibitory),
            width = 9, height = 4)

# Ridge / density plots (exhausted clusters)
save_figure("densities_IgE_exhausted.png",
            plot_ridge_composite(sce_IgE_clean, markers_ex_CD14,
                                 clusters_to_show = c("11", "12", "13", "14")),
            width = 9, height = 6)

# ---- 8. Pseudotime (whole population) ------------------------------------
cells_ok <- !is.na(reducedDim(sce_IgE_clean, "TSNE")[, 1])
sce_tsne <- sce_IgE_clean[, cells_ok]
sce_tsne$merged_ids <- cluster_ids(sce_tsne, "cluster_ID")

sce_tsne <- slingshot(sce_tsne, clusterLabels = "merged_ids",
                       reducedDim = "TSNE",
                       start.clus = c("1", "2", "3"), reweight = FALSE)

# Trajectory paths (from slingshot output, manually curated)
paths_all <- list(
  c("1", "4", "11", "13", "10", "14", "8", "9"),
  c("1", "4", "11", "13", "10", "6", "7", "3"),
  c("1", "4", "11", "13", "10", "5"),
  c("1", "4", "11", "12"),
  c("1", "4", "2")
)
save_figure("tSNE_IgE_clean_trajectory.png",
            plot_trajectory_mst(sce_tsne, paths_all))

# Pseudotime boxplot
df_pt <- as.data.frame(reducedDim(sce_tsne, "TSNE"))
colnames(df_pt) <- c("tSNE1", "tSNE2")
df_pt$Cluster    <- as.character(cluster_ids(sce_tsne, "cluster_ID"))
df_pt$Pseudotime <- rowMeans(slingPseudotime(sce_tsne), na.rm = TRUE)

save_figure("pseudotime_boxplot.png",
            ggplot(df_pt, aes(x = reorder(Cluster, Pseudotime, FUN = median),
                              y = Pseudotime, fill = Cluster)) +
              geom_boxplot(outlier.size = 0.1, alpha = 0.7) +
              scale_fill_manual(values = cluster_colors) +
              theme_minimal() +
              labs(title = "Pseudotime Progression by Cluster",
                   x = "Clusters (Ordered by Age)", y = "Pseudotime") +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Marker kinetics over pseudotime
pt_expr <- as.data.frame(t(assay(sce_tsne, "exprs")))
pt_expr$Pseudotime <- sce_tsne$pseudotime
pt_long <- pt_expr %>%
  pivot_longer(-Pseudotime, names_to = "variable", values_to = "value")

save_figure("marker_kinetics_pseudotime.png",
            ggplot(pt_long, aes(x = Pseudotime, y = value, color = variable)) +
              geom_smooth(method = "gam", se = FALSE) +
              facet_wrap(~variable, scales = "free_y") +
              theme_minimal() +
              labs(title = "Marker Expression Kinetics",
                   y = "Arcsinh Expression", x = "Pseudotime"),
            width = 12, height = 9)

# Classical monocytes only
keep_cl <- c("1", "2", "4", "8", "9", "11", "12", "13")
sce_classical <- sce_tsne[, sce_tsne$merged_ids %in% keep_cl]
sce_classical <- slingshot(sce_classical, clusterLabels = "merged_ids",
                            reducedDim = "TSNE",
                            start.clus = "1", reweight = FALSE)

paths_classical <- list(
  c("1", "4", "11", "12", "9", "8"),
  c("1", "4", "11", "13"),
  c("1", "4", "2")
)
save_figure("tSNE_IgE_clean_classical_trajectory.png",
            plot_trajectory_mst(sce_classical, paths_classical,
                                title = "Classical Monocyte Clusters on tSNE"))

# ---- 9. Per-setting analyses (2D / 3D) -----------------------------------
for (setting_val in c("2D", "3D")) {
  message("===== Setting: ", setting_val, " =====")
  sce_sub <- filterSCE(sce_IgE_clean, setting == setting_val)
  tag <- paste0("IgE_clean_", setting_val)

  # tSNE
  save_figure(paste0("tSNE_", tag, ".png"),
              plotDR(sce_sub, dr = "TSNE", color_by = "cluster_ID") +
                scale_color_manual(values = cluster_colors) + theme_pub())

  # Abundance
  for (gv in c("condition", "sample_type")) {
    save_figure(paste0("abundance_", tag, "_", gv, ".png"),
                plotAbundances(sce_sub, k = "cluster_ID",
                               by = "sample_id", group_by = gv) +
                  scale_fill_manual(values = cluster_colors) +
                  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        text = element_text(size = 24),
                        plot.title = element_text(face = "bold", size = 28)),
                width = 10, height = 5)
  }

  # Antibody comparisons (condition only, no setting covariate)
  diff_sub <- run_antibody_comparisons(
    sce_sub, cols_design = c("condition", "sample_type"))

  message("  DA NIP vs IgE (", setting_val, ")")
  print(topTable(diff_sub$NIP_vs_IgE$da, format_vals = TRUE))
  message("  DA PBS vs IgE (", setting_val, ")")
  print(topTable(diff_sub$PBS_vs_IgE$da, format_vals = TRUE))

  # HV vs M within each condition
  hv_m_results <- list()
  for (cond in c("PBS", "NIP", "CSPG4")) {
    hv_m_results[[cond]] <- tryCatch(
      run_hv_vs_m(sce_sub, cond, cols_design = c("condition", "sample_type")),
      error = function(e) { message("  Skipping HV vs M for ", cond, ": ", e$message); NULL }
    )
  }

  # Abundance with antibody p-values
  clusters_show <- if (setting_val == "2D") c(11, 12, 13, 14) else c(1, 2, 3)
  save_figure(paste0("abundance_", tag, "_padj.png"),
              plot_abundance_condition_pval(
                sce_sub,
                da_nip = diff_sub$NIP_vs_IgE$da,
                da_pbs = diff_sub$PBS_vs_IgE$da,
                clusters_to_show = clusters_show,
                shape_by = "sample_type"),
              width = 10, height = 4)

  # Per-cluster marker expression
  message("  Per-cluster marker plots (", setting_val, ")...")
  plot_cluster_markers_loop(
    sce_sub,
    ds_nip  = diff_sub$NIP_vs_IgE$ds,
    ds_pbs  = diff_sub$PBS_vs_IgE$ds,
    markers = markers_ex_CD14,
    output_subdir = paste0("pbExprs_cluster_", tag),
    shape_by = "sample_type")

  # Ridge / density plots
  ridge_clusters <- if (setting_val == "2D") c("1", "2", "3") else c("11", "12", "13", "14")
  save_figure(paste0("densities_", tag, ".png"),
              plot_ridge_composite(sce_sub, markers_ex_CD14,
                                   clusters_to_show = ridge_clusters),
              width = 9, height = 6)
}

# ---- 10. 2D vs 3D within each condition ----------------------------------
for (cond_val in c("PBS", "NIP IgE", "CSPG4 IgE")) {
  cond_tag <- gsub(" ", "_", cond_val)
  message("===== 2D vs 3D: ", cond_val, " =====")

  diff_setting <- run_2d_vs_3d(sce_IgE_clean, cond_val)

  message("  DA 2D vs 3D (", cond_val, ")")
  print(topTable(diff_setting$da, format_vals = TRUE))
  message("  DS 2D vs 3D (", cond_val, ")")
  print(topTable(diff_setting$ds, format_vals = TRUE))
}

# Save final SCE
saveRDS(sce_IgE_clean, file.path(rds_dir, "sce_IgE_clean_final.rds"))
message("Analysis complete.")
