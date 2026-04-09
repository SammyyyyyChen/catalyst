# ---- Save figure helper --------------------------------------------------

save_figure <- function(filename, plot_obj, width = 8, height = 7,
                        fig_dir = file.path("output", "figures")) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(fig_dir, filename)
  png(path, width = width, height = height, units = "in", res = 300)
  print(plot_obj)
  dev.off()
  message("Saved: ", path)
}

# ---- QC: cell counts bar plot -------------------------------------------

plot_cell_counts <- function(sce, collapse_ids = c("HV814" = "HV815")) {
  df <- as.data.frame(colData(sce))
  for (nm in names(collapse_ids))
    df <- df %>% mutate(patient_id = fct_collapse(patient_id, !!nm := c(nm, collapse_ids[[nm]])))

  ggplot(df, aes(x = patient_id, fill = condition)) +
    geom_bar(color = "white", linewidth = 0.1) +
    scale_fill_manual(values = condition_colors) +
    stat_count(geom = "text", aes(label = after_stat(count)),
               vjust = 1.2, size = 4, fontface = "bold", colour = "white") +
    theme_pub() +
    labs(title = "Cell Counts per Patient by Condition",
         x = "Patient ID", y = "Total Cell Count", fill = "Condition") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          legend.position = "bottom") +
    facet_wrap(~setting, scales = "free")
}

# ---- Abundance boxplot with antibody-comparison p-values -----------------
# Overlays NIP-vs-IgE and PBS-vs-IgE significance brackets on
# CATALYST plotAbundances output.

plot_abundance_condition_pval <- function(sce, da_nip, da_pbs,
                                         clusters_to_show = NULL,
                                         shape_by = "setting",
                                         ncol_facet = 5) {
  p <- plotAbundances(sce, k = "cluster_ID", by = "cluster_id", shape_by = shape_by) +
    facet_wrap(~cluster_id, scales = "free_y", ncol = 7) +
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = alpha(condition_colors, 0.6)) +
    theme(text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 24),
          axis.title = element_text(size = 20))

  if (!is.null(clusters_to_show))
    p$data <- p$data %>% dplyr::filter(cluster_id %in% clusters_to_show)

  s_nip <- prep_da_stats(da_nip, "nip")
  s_pbs <- prep_da_stats(da_pbs, "pbs")

  lbl <- p$data %>%
    group_by(cluster_id) %>%
    summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE),
              .groups = "drop") %>%
    left_join(s_nip, by = "cluster_id") %>%
    left_join(s_pbs, by = "cluster_id") %>%
    mutate(cluster_id = factor(cluster_id, levels = levels(p$data$cluster_id)),
           y_pbs = local_max * 1.15, y_pbs_t = local_max * 1.17,
           y_nip = local_max * 1.05, y_nip_t = local_max * 1.07)

  p +
    facet_wrap(~cluster_id, scales = "free_y", ncol = ncol_facet) +
    # PBS vs CSPG4 IgE bracket (wide)
    geom_segment(data = dplyr::filter(lbl, p_label_pbs != ""),
                 aes(x = 1, xend = 3, y = y_pbs, yend = y_pbs),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = dplyr::filter(lbl, p_label_pbs != ""),
              aes(x = 2, y = y_pbs_t, label = p_label_pbs),
              inherit.aes = FALSE, size = 4, vjust = 0) +
    # NIP IgE vs CSPG4 IgE bracket (narrow)
    geom_segment(data = dplyr::filter(lbl, p_label_nip != ""),
                 aes(x = 2, xend = 3, y = y_nip, yend = y_nip),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = dplyr::filter(lbl, p_label_nip != ""),
              aes(x = 2.5, y = y_nip_t, label = p_label_nip),
              inherit.aes = FALSE, size = 4, vjust = 0) +
    theme(panel.spacing.x = unit(1, "cm"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          text = element_text(size = 24),
          plot.title = element_text(face = "bold", size = 24))
}

# ---- Abundance boxplot split by a secondary factor (HV/M or 2D/3D) ------
# Shows boxplots per condition, split by `split_by`, with within-condition
# significance brackets from a list of DA results.

plot_abundance_split <- function(sce, da_list, clusters_to_show,
                                 split_by = "sample_type") {
  p <- plotAbundances(sce, k = "cluster_ID", by = "cluster_id",
                      shape_by = split_by)
  plot_data <- p$data %>% dplyr::filter(cluster_id %in% clusters_to_show)

  alpha_vals <- if (split_by == "sample_type") {
    c("HV" = 0.1, "M" = 0.3)
  } else {
    c("2D" = 0.1, "3D" = 0.5)
  }
  shape_vals <- if (split_by == "sample_type") {
    c("HV" = 21, "M" = 24)
  } else {
    c("2D" = 21, "3D" = 24)
  }

  # Build label data from the list of DA results
  n_comparisons <- length(da_list)
  lbl <- plot_data %>%
    group_by(cluster_id) %>%
    summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE),
              .groups = "drop")
  for (i in seq_along(da_list)) {
    si <- prep_da_stats(da_list[[i]], as.character(i))
    lbl <- left_join(lbl, si, by = "cluster_id")
  }
  lbl <- lbl %>%
    mutate(cluster_id = factor(cluster_id, levels = levels(plot_data$cluster_id)))
  for (i in seq_along(da_list))
    lbl <- lbl %>%
      mutate("y{i}"   := local_max * (1 + 0.05 * i),
             "y{i}_t" := local_max * (1 + 0.05 * i + 0.02),
             .names_repair = "minimal")

  p_out <- ggplot(plot_data, aes(x = condition, y = Freq, color = condition)) +
    geom_boxplot(aes(fill = condition, alpha = .data[[split_by]]),
                 width = 0.6, position = position_dodge(width = 0.8),
                 outlier.shape = NA, show.legend = FALSE) +
    geom_point(aes(shape = .data[[split_by]], group = .data[[split_by]]),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
               size = 1.5, stroke = 1, fill = "transparent", show.legend = TRUE) +
    facet_wrap(~cluster_id, scales = "free_y", ncol = 4) +
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    scale_alpha_manual(values = alpha_vals) +
    scale_shape_manual(values = shape_vals)

  # Add per-condition significance brackets
  for (i in seq_along(da_list)) {
    col_p <- paste0("p_label_", i)
    col_y <- paste0("y", i)
    col_yt <- paste0("y", i, "_t")
    # Compute y columns inline
    tier_lbl <- lbl %>%
      mutate(y_bar  = local_max * (1 + 0.05 * !!i),
             y_text = local_max * (1 + 0.05 * !!i + 0.02)) %>%
      dplyr::filter(.data[[col_p]] != "")
    x_center <- i - 0.2
    p_out <- p_out +
      geom_segment(data = tier_lbl,
                   aes(x = !!i - 0.2, xend = !!i + 0.2, y = y_bar, yend = y_bar),
                   inherit.aes = FALSE, color = "black") +
      geom_text(data = tier_lbl,
                aes(x = !!i, y = y_text, label = .data[[col_p]]),
                inherit.aes = FALSE, size = 5, color = "black")
  }

  p_out +
    theme_bw() +
    labs(y = "Proportion [%]", x = NULL) +
    theme(text = element_text(size = 20),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 20),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.spacing = unit(1, "cm"))
}

# ---- Per-cluster marker expression loop with p-value brackets ------------
# Generates one PNG per cluster with pseudobulk marker boxplots and
# NIP-vs-IgE / PBS-vs-IgE brackets.

plot_cluster_markers_loop <- function(sce, ds_nip, ds_pbs,
                                      markers, output_subdir,
                                      shape_by = "setting",
                                      clusters = as.character(1:14),
                                      fig_dir = file.path("output", "figures")) {
  out_dir <- file.path(fig_dir, output_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  s_nip <- prep_ds_stats(ds_nip, "nip")
  s_pbs <- prep_ds_stats(ds_pbs, "pbs")

  for (id in clusters) {
    message("  cluster ", id)
    sce_sub <- filterSCE(sce, cluster_ID == id)
    p <- plotPbExprs(sce_sub, k = "cluster_ID", features = markers,
                     facet_by = "cluster_id", shape_by = shape_by, ncol = 4) +
      scale_color_manual(values = condition_colors) +
      scale_fill_manual(values = condition_colors) +
      theme(text = element_text(size = 20),
            plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

    p$data <- p$data %>% dplyr::filter(cluster_id == id)
    vals <- as.numeric(as.character(p$data$value))
    y_max <- max(vals, na.rm = TRUE)
    y_min <- min(vals, na.rm = TRUE)
    jump  <- (y_max - y_min) * 0.12

    lbl <- p$data %>%
      mutate(value = as.numeric(as.character(value))) %>%
      group_by(antigen, cluster_id) %>%
      summarize(local_max = max(value, na.rm = TRUE), .groups = "drop") %>%
      left_join(s_nip, by = c("antigen" = "marker_id", "cluster_id")) %>%
      left_join(s_pbs, by = c("antigen" = "marker_id", "cluster_id")) %>%
      mutate(y_nip = local_max + jump * 0.5,  y_nip_t = local_max + jump * 0.7,
             y_pbs = local_max + jump * 1.5,  y_pbs_t = local_max + jump * 1.7)

    p$layers[[1]]$position <- position_dodge(width = 0.9)
    p$layers[[1]]$geom_params$width <- 0.6
    p$layers[[2]]$position <- position_dodge(width = 0.9)

    marker_levels <- levels(factor(p$data$antigen))
    lbl <- lbl %>%
      mutate(x_pos  = as.numeric(factor(antigen, levels = marker_levels)),
             x_s_nip = x_pos + offsets[2], x_e_nip = x_pos + offsets[3],
             x_s_pbs = x_pos + offsets[1], x_e_pbs = x_pos + offsets[3])

    p_final <- p +
      geom_segment(data = dplyr::filter(lbl, p_label_pbs != ""),
                   aes(x = x_s_pbs, xend = x_e_pbs, y = y_pbs, yend = y_pbs),
                   inherit.aes = FALSE) +
      geom_text(data = dplyr::filter(lbl, p_label_pbs != ""),
                aes(x = (x_s_pbs + x_e_pbs)/2, y = y_pbs_t, label = p_label_pbs),
                inherit.aes = FALSE, size = 3.5, vjust = 0) +
      geom_segment(data = dplyr::filter(lbl, p_label_nip != ""),
                   aes(x = x_s_nip, xend = x_e_nip, y = y_nip, yend = y_nip),
                   inherit.aes = FALSE) +
      geom_text(data = dplyr::filter(lbl, p_label_nip != ""),
                aes(x = (x_s_nip + x_e_nip)/2, y = y_nip_t, label = p_label_nip),
                inherit.aes = FALSE, size = 3.5, vjust = 0) +
      coord_cartesian(ylim = c(y_min, y_max + jump * 2.5))

    fname <- paste0("pbExprs_cluster_", id, "_padj.png")
    save_figure(fname, p_final, width = 15, height = 5,
                fig_dir = out_dir)
    rm(sce_sub, p, p_final); gc(verbose = FALSE)
  }
}

# ---- Global marker expression with p-values -----------------------------

plot_global_markers_pval <- function(sce, ds_nip, ds_pbs, markers) {
  p <- plotPbExprs(sce, features = markers, facet_by = "antigen",
                   shape_by = "setting", ncol = 7) +
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    theme(text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank())

  # Format stats (global has no cluster_id, just marker_id)
  fmt <- function(obj, nm) {
    res <- as.data.frame(rowData(obj$res))
    res %>%
      select(marker_id, p_adj) %>%
      mutate(!!paste0("p_label_", nm) :=
               ifelse(p_adj < 0.0001, "<0.0001",
                      ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
      select(-p_adj)
  }
  s_nip <- fmt(ds_nip, "nip")
  s_pbs <- fmt(ds_pbs, "pbs")

  lbl <- p$data %>%
    group_by(antigen) %>%
    summarize(local_max = max(as.numeric(as.character(value)), na.rm = TRUE),
              .groups = "drop") %>%
    left_join(s_nip, by = c("antigen" = "marker_id")) %>%
    left_join(s_pbs, by = c("antigen" = "marker_id")) %>%
    mutate(antigen = factor(antigen, levels = levels(p$data$antigen)),
           y_pbs = local_max * 1.10, y_pbs_t = local_max * 1.12,
           y_nip = local_max * 1.20, y_nip_t = local_max * 1.22)

  p +
    facet_wrap(~antigen, scales = "free_y", ncol = 7) +
    geom_segment(data = dplyr::filter(lbl, p_label_nip != ""),
                 aes(x = 1, xend = 3, y = y_nip, yend = y_nip),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = dplyr::filter(lbl, p_label_nip != ""),
              aes(x = 2, y = y_nip_t, label = p_label_nip),
              inherit.aes = FALSE, size = 3.5, vjust = 0) +
    geom_segment(data = dplyr::filter(lbl, p_label_pbs != ""),
                 aes(x = 2, xend = 3, y = y_pbs, yend = y_pbs),
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = dplyr::filter(lbl, p_label_pbs != ""),
              aes(x = 2.5, y = y_pbs_t, label = p_label_pbs),
              inherit.aes = FALSE, size = 3.5, vjust = 0)
}

# ---- Ridge / density plot row -------------------------------------------

make_ridge_row <- function(data, markers, title_text,
                           clusters_to_show, ref_medians,
                           show_x = FALSE) {
  df_sub <- data %>%
    filter(antigen %in% markers, cluster_id %in% clusters_to_show) %>%
    mutate(antigen = factor(antigen, levels = markers))
  meds_sub <- ref_medians %>%
    filter(antigen %in% markers) %>%
    mutate(antigen = factor(as.character(antigen), levels = markers))

  p <- ggplot(df_sub, aes(x = expression, y = cluster_id,
                           fill = cluster_id, color = cluster_id)) +
    geom_vline(data = meds_sub, aes(xintercept = m_val),
               color = "grey20", linetype = "dashed", linewidth = 0.5) +
    geom_density_ridges(alpha = 0.7, scale = 1.5, linewidth = 0.5) +
    scale_fill_manual(values = cluster_colors) +
    scale_color_manual(values = cluster_colors) +
    facet_wrap(~antigen, nrow = 1, ncol = 5, scales = "free_x", drop = FALSE) +
    facetted_pos_scales(
      x = list(antigen %in% c(" ", "  ") ~ scale_x_continuous(guide = "none"))
    ) +
    labs(subtitle = title_text, x = NULL, y = "Cluster ID") +
    theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 14),
          strip.background = element_blank(),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 14, margin = margin(b = 5)),
          panel.background = element_blank(),
          legend.position = "none",
          panel.spacing = unit(1, "lines"))

  if (!show_x)
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank())
  p
}

# Build the full 3-row ridge composite for a given SCE subset.
plot_ridge_composite <- function(sce, markers, clusters_to_show,
                                 ref_cluster = "1") {
  plot_df <- as.data.frame(t(assay(sce, "exprs")[markers, ]))
  cl_ids  <- cluster_ids(sce, "cluster_ID")
  plot_df$cluster_id <- factor(as.character(cl_ids),
                                levels = rev(as.character(1:14)))
  plot_df_long <- plot_df %>%
    pivot_longer(cols = -cluster_id, names_to = "antigen", values_to = "expression")
  plot_df_long$antigen <- factor(plot_df_long$antigen, levels = markers)

  medians_ref <- plot_df_long %>%
    filter(cluster_id == ref_cluster) %>%
    group_by(antigen) %>%
    summarize(m_val = median(expression, na.rm = TRUE), .groups = "drop")

  p1 <- make_ridge_row(plot_df_long, ridge_groups[[1]], names(ridge_groups)[1],
                        clusters_to_show, medians_ref)
  p2 <- make_ridge_row(plot_df_long, ridge_groups[[2]], names(ridge_groups)[2],
                        clusters_to_show, medians_ref)
  p3 <- make_ridge_row(plot_df_long, ridge_groups[[3]], names(ridge_groups)[3],
                        clusters_to_show, medians_ref, show_x = TRUE)

  p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
}

# ---- Trajectory / MST on tSNE -------------------------------------------

plot_trajectory_mst <- function(sce_tsne, paths, title = "Monocyte Clusters on tSNE",
                                subtitle = "Trajectory: Resting -> Activation -> Exhaustion -> Inhibition") {
  df_plot <- as.data.frame(reducedDim(sce_tsne, "TSNE"))
  colnames(df_plot) <- c("tSNE1", "tSNE2")
  df_plot$Cluster <- as.character(cluster_ids(sce_tsne, "cluster_ID"))

  centers <- df_plot %>%
    group_by(Cluster) %>%
    summarize(x = mean(tSNE1), y = mean(tSNE2), .groups = "drop") %>%
    as.data.frame()
  rownames(centers) <- centers$Cluster

  edge_coords <- data.frame()
  for (path in paths) {
    for (i in seq_len(length(path) - 1)) {
      s <- path[i]; e <- path[i + 1]
      if (s %in% centers$Cluster && e %in% centers$Cluster)
        edge_coords <- rbind(edge_coords, data.frame(
          x = centers[s, "x"], y = centers[s, "y"],
          xend = centers[e, "x"], yend = centers[e, "y"]))
    }
  }
  edge_coords <- unique(edge_coords)

  ggplot(df_plot, aes(x = tSNE1, y = tSNE2)) +
    geom_point(aes(color = Cluster), size = 0.5) +
    geom_segment(data = edge_coords,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "black", linewidth = 1.2, inherit.aes = FALSE) +
    geom_label(data = centers, aes(x = x, y = y, label = Cluster),
               size = 4, fontface = "bold", alpha = 0.8, inherit.aes = FALSE) +
    theme_minimal() +
    theme(text = element_text(size = 14),
          plot.title = element_text(face = "bold", size = 20),
          axis.title = element_text(size = 18),
          plot.subtitle = element_text(size = 14)) +
    scale_color_manual(values = cluster_colors) +
    labs(title = title, subtitle = subtitle)
}
