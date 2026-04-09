library(CATALYST)
library(flowCore)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(stringr)
library(dplyr)
library(diffcyt)
library(viridis)
library(pals)
library(RColorBrewer)
library(cyCombine)
library(ComplexHeatmap)
library(forcats)
library(ggridges)
library(tidyselect)

# Save RDS
saveRDS(sce_IgE, file = "sce_IgE.rds")
saveRDS(sce, file = "sce.rds")
saveRDS(sce_IgE_filter, file = "sce_IgE_filter.rds")
saveRDS(sce_IgE_clean, file = "sce_IgE_clean.rds")
sce_IgE <- readRDS("sce_IgE.rds")
sce_IgE_filter <- readRDS("sce_IgE_filter.rds")
#=========================================================================

#Load gated monocytic cells FCS files names
file_path <- "raw files test"
fileNames <- list.files(
  path = file_path,
  pattern = "*.fcs",
  full.names = TRUE,
  ignore.case = TRUE
)

# Loop through each file to check for the Header/Text offset error
bad_files <- c()

for (f in fileNames) {
  # Try reading the header and data
  check <- tryCatch({
    read.FCS(f, transformation = FALSE)
    TRUE
  }, error = function(e) {
    message(paste("!!! ERROR in file:", basename(f)))
    message(paste("    Message:", e$message))
    return(FALSE)
  })
  
  if (!check) {
    bad_files <- c(bad_files, f)
  }
}

# Final summary
if (length(bad_files) == 0) {
  message("All files read successfully!")
} else {
  message(paste("\nFound", length(bad_files), "corrupted files:"))
  print(basename(bad_files))
}



#Creating metadata file, sample_id has to be unique
file_base <- basename(fileNames)
metadata <- data.frame(file_name = file_base) %>%
  mutate(
    setting = str_extract(file_base, "2D|3D"),
    sample_type = str_extract(file_base, "HV|M"),
    patient_id = str_extract(file_base, "(HV|M)[0-9]+"),
    condition = str_extract(file_base, "PBS|NIP IgG1|CSPG4 IgG1|NIP IgE|CSPG4 IgE"),
    sample_id = paste(patient_id, condition, setting, sep = "_")
  ) %>%
  mutate(batch = case_when(
    patient_id %in% c("HV778", "HV807", "HV808") ~ "20251117",
    patient_id %in% c("HV810", "HV814", "HV815") ~ "20251205",
    patient_id %in% c("M764", "M770", "M781") ~ "20251212",
    patient_id %in% c("M817", "M818", "HV820") ~ "20251219",
    patient_id %in% c("M748", "M832", "HV830") ~ "20250131",
    patient_id %in% c("M834", "HV835", "HV836") ~ "20250206"))

head(metadata)
write.csv(metadata, "metadata.csv")

# Read one FCS file to extract parameters
ff <- read.FCS(fileNames[1], transformation = FALSE, truncate_max_range = FALSE)
params <- pData(parameters(ff))

panel <- tibble(
  fcs_colname = params$name,
  antigen = params$desc
) %>% mutate(
  antigen = gsub("[[:space:]]+[A-Za-z0-9-]+-A$", "", antigen)
)

# Clean up antigen names if needed
panel$antigen[is.na(panel$antigen)] <- panel$fcs_colname[is.na(panel$antigen)]

# Mark which channels are used for clustering
panel <- panel %>%
  mutate(
    marker_class = case_when(
      antigen %in% c("FSC-A", "SSC-A", "Time", "DAPI-A", "FSC-Width") ~ "none",
      TRUE ~ "state"
    )
  )
write.csv(panel, "panel.csv")

# Assign colour
p <- plotDR(sce_IgE_clean, dr = "TSNE", color_by = "meta16")
g <- ggplot_build(p)
cluster_colors <- unique(g$data[[1]]$colour)
print(cluster_colors)

condition_colors <- c(
  "PBS"        = "#BDBDBD",   # Neutral Grey
  "NIP IgE"    = "#e3c7b6", # Pale Red/Pink
  "CSPG4 IgE"  = "#9C2830", # Dark Red
  "2D" = "#89d6a4",
  "3D" = "#f5b871",
  "HV" = "lightskyblue2",
  "M" = "coral1"
)
cluster_colors <- c(
  "1"        = "#FB8072",   # Neutral Grey
  "2"   = "#E78AC3", # Pale Blue
  "3" = "#FDB462", # Dark Blue
  "4"    = "#DC050C", # Pale Red/Pink
  "5"  = "#33A02C", # Dark Red
  "6" = "#E7298A",
  "7" = "#B2DF8A",
  "8" = "#1D74CC",
  "9" = "#882E72",
  "10" = "#FF7F00",
  "11" = "#8DD3C7",
  "12" = "#7BAFDE",
  "13" = "#E6AB02",
  "14" = "#B17BA6"
)

[1] "#B2DF8A" "#FB8072" "#882E72" "#DC050C" "#55A1B1" "#E78AC3" "#7BAFDE" "#8DD3C7" "#1965B0" "#FDB462" "#B17BA6"
[12] "#A6761D" "#33A02C" "#E7298A" "#E6AB02" "#FF7F00"
#=========================================================================
# Import data into SCE
sce <- prepData(
  x = fileNames,
  panel = panel,
  md = metadata,
  features = panel$fcs_colname[panel$marker_class == "state"],
  cofactor = 400,
  md_cols = list(
    file = "file_name", 
    id = "sample_id",
    factors = c("sample_type", "patient_id", "condition", "setting", "batch")
  ), ignore.text.offset = TRUE, emptyValue = FALSE
)

sce_IgE <- filterSCE(sce, condition %in% c("PBS", "NIP IgE", "CSPG4 IgE"))

# Set PBS as the baseline for condition
metadata(sce_IgE)$experiment_info$condition <- factor(
  metadata(sce_IgE)$experiment_info$condition, 
  levels = c("PBS", "NIP IgE", "CSPG4 IgE")
)
sce_IgE$condition <- factor(sce_IgE$condition, levels = c("PBS", "NIP IgE", "CSPG4 IgE"))
# Set 2D as the baseline for setting
metadata(sce_IgE)$experiment_info$setting <- factor(
  metadata(sce_IgE)$experiment_info$setting, 
  levels = c("2D", "3D")
)
# Set HV as the baseline for sample type
metadata(sce_IgE)$experiment_info$sample_type <- factor(
  metadata(sce_IgE)$experiment_info$sample_type, 
  levels = c("HV", "M")
)


# Check sce project
table(colData(sce_IgE)$patient_id)
table(colData(sce_IgE)$condition)
table(colData(sce_IgE)$setting)

df_meta_IgE <- as.data.frame(colData(sce_IgE))
df_meta_IgE_combine <- as.data.frame(colData(sce_IgE)) %>%
  mutate(patient_id = fct_collapse(patient_id, "HV814" = c("HV814", "HV815")))

# QC Create the stacked bar plot
png("IgE cell counts per patient by condition.png", width = 19, height = 9, units = "in", res = 300)
ggplot(df_meta_IgE_combine, aes(x = patient_id, fill = condition)) +
  geom_bar(color = "white", linewidth = 0.1) + 
  scale_fill_manual(values = condition_colors) + 
  stat_count(geom = "text", 
             aes(label = after_stat(count)), 
             vjust = 1.2, 
             size = 4,
             fontface = "bold",
             colour = "white") +
  theme_minimal() +
  labs(title = "Cell Counts per Patient by Condition", 
       x = "Patient ID", 
       y = "Total Cell Count",
       fill = "Condition") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.position = "bottom"
  ) +
  facet_wrap(~setting, scales = "free")
dev.off()

# MDS_by_condition_setting
png("MDS_IgE_by_condition_setting.png", width = 8, height = 6, units = "in", res = 300)
pbMDS(sce_IgE, color_by = "condition", shape_by = "setting", label_by = NULL, size_by = TRUE) + 
  scale_color_manual(values = condition_colors) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

# MDS_by_batch
png("MDS_by_IgE_condition_batch.png", width = 8, height = 6, units = "in", res = 300)
pbMDS(sce_IgE, color_by = "batch", label_by = NULL, shape_by = "condition", size_by = TRUE) + 
  #scale_color_manual(values = condition_colors) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

#NRS_by_condition
png("NRS_IgE_by_condition.png", width = 12, height = 6, units = "in", res = 300)
plotNRS(sce_IgE, features = "state", color_by = "condition") + 
  scale_color_manual(values = condition_colors) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

#================================================================================================
# Remove sample with weird NIP
sce_IgE_clean <- filterSCE(sce_IgE, batch %in% c("20251117", "20251205", "20251212", "20250131", "20250206"))
saveRDS(sce_IgE_clean, file = "sce_IgE_clean.rds")

# Check sce project
table(colData(sce_IgE_clean)$patient_id)
table(colData(sce_IgE_clean)$condition)
table(colData(sce_IgE_clean)$setting)

# QC on clean data
df_meta_IgE_clean <- as.data.frame(colData(sce_IgE_clean)) %>%
  mutate(patient_id = fct_collapse(patient_id, "HV814" = c("HV814", "HV815")))
png("IgE clean cell counts per patient by condition.png", width = 19, height = 9, units = "in", res = 300)
ggplot(df_meta_IgE_clean, aes(x = patient_id, fill = condition)) +
  geom_bar(color = "white", linewidth = 0.1) + 
  scale_fill_manual(values = condition_colors) + 
  stat_count(geom = "text", 
             aes(label = after_stat(count)), 
             vjust = 1.2, 
             size = 4,
             fontface = "bold",
             colour = "white") +
  theme_minimal() +
  labs(title = "Cell Counts per Patient by Condition", 
       x = "Patient ID", 
       y = "Total Cell Count",
       fill = "Condition") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.position = "bottom"
  ) +
  facet_wrap(~setting, scales = "free")
dev.off()
# MDS_by_condition_setting
png("MDS_IgE_clean_by_condition_setting.png", width = 8, height = 6, units = "in", res = 300)
pbMDS(sce_IgE_clean, color_by = "condition", shape_by = "setting", label_by = NULL, size_by = TRUE) + 
  scale_color_manual(values = condition_colors) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

# MDS_by_batch
png("MDS_by_IgE_clean_condition_batch.png", width = 8, height = 6, units = "in", res = 300)
pbMDS(sce_IgE_clean, color_by = "batch", label_by = NULL, shape_by = "condition", size_by = TRUE) + 
  #scale_color_manual(values = condition_colors) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

#NRS_by_condition
png("NRS_IgE_clean_by_condition.png", width = 12, height = 6, units = "in", res = 300)
plotNRS(sce_IgE_clean, features = "state", color_by = "condition") + 
  scale_color_manual(values = condition_colors) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()


#==========================================================================================
# Run clustering on sce_IgE
sce_IgE <- cluster(sce_IgE, features = "state", 
                   xdim = 10, ydim = 10, maxK = 12, seed = 123)
# delta plot
png("IgE_delta_plot.png", width = 12, height = 9, units = "in", res = 300)
delta_area(sce_IgE) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

table(colData(sce_IgE)$cluster_id)

# Run clustering on sce_IgE_clean
rowData(sce_IgE_clean)$marker_class <- as.character(rowData(sce_IgE_clean)$marker_class)
rowData(sce_IgE_clean)$marker_class[rownames(sce_IgE_clean) == "CD32B"] <- "state"
rowData(sce_IgE_clean)$marker_class <- factor(rowData(sce_IgE_clean)$marker_class)
rowData(sce_IgE_clean)

sce_IgE_clean <- cluster(sce_IgE_clean, features = "state", 
                   xdim = 10, ydim = 10, maxK = 16, seed = 123)
# delta plot
png("IgE_clean_delta_plot.png", width = 12, height = 9, units = "in", res = 300)
delta_area(sce_IgE_clean) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

table(colData(sce_IgE_clean)$cluster_id)
#==========================================================================================

#assign markers group
ordered_markers <- c("CD14", "CD16", "CD64", "CD32B", "FCER1", "CD23", "CD80", "CD86", "CD40", "CCR2", "HLADR", "CD163", "CD206", "PDL1")
FcRs <- c("FCER1", "CD23")
Activation <- c("CD80", "CD86", "CD40", "CCR2", "HLADR")
Inhibitory <- c("CD163", "CD206", "PDL1")


p1 <- plotExprHeatmap(sce_IgE, features = ordered_markers, k = "meta12",by = "cluster_id", fun = "median", row_anno = TRUE,
                scale = "last", q = 0, bars = TRUE, row_clust = FALSE, col_clust = FALSE, perc = TRUE)
p2 <- plotExprHeatmap(sce_IgE_clean, features = ordered_markers, k = "meta16",by = "cluster_id", fun = "median", row_anno = TRUE,
                scale = "last", q = 0, bars = TRUE, row_clust = FALSE, col_clust = FALSE, perc = TRUE)
p1 + p2
png("IgE_clean_exprs_heatmap_by_meta16.png", width = 12, height = 7, units = "in", res = 300)
p2
dev.off()

# Merging 1 (for sce_IgE)
merging_table_1<- data.frame(
  cluster = as.character(1:12),
  merged_cluster = c(
    "2", # 1 HLADR+ PDL1+
    "10",  # 2 HLADR+
    "7",   # 3 CD23+ PDL1+
    "8", # 4 CD40+ CD80+ CD86+
    "9",   # 5 PDL1+
    "8", # 6 CD163+ CD206+ PDL1+
    "3",  # 7 HLADR+
    "11",  # 8 HLADR+
    "4",   # 9 CD163+ CD206+ 
    "5", # 10 CD163+ CD206+ PDL1+
    "6", # 11 nothing
    "1" # 12
  ))
numeric_order <- as.character(1:11)
merging_table_1$merged_cluster <- factor(
  merging_table_1$merged_cluster, 
  levels = numeric_order
)
sce_IgE <- mergeClusters(sce_IgE, k = "meta12", table = merging_table_1, id = "merging1", overwrite = TRUE)
sce_IgE_filter <- filterSCE(sce_IgE, k = "merging1", cluster_id %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

# Merging 2 (for sce_IgE_clean)
merging_table_2<- data.frame(
  cluster = as.character(1:16),
  merged_cluster = c(
    "10", # 1 
    "11",  # 2 
    "12",   # 3 
    "8", # 4 
    "12",   # 5 
    "2", # 6 
    "5",  # 7 
    "9",  # 8 
    "13",   # 9 
    "1", # 10
    "6", # 11 
    "14", # 12
    "4",# 13
    "7", # 14
    "3", # 15
    "3" # 16
  ))
numeric_order <- as.character(1:14)
merging_table_2$merged_cluster <- factor(
  merging_table_2$merged_cluster, 
  levels = numeric_order
)
sce_IgE_clean <- mergeClusters(sce_IgE_clean, k = "meta16", table = merging_table_2, id = "cluster_ID", overwrite = TRUE)

# Merging 3 (for sce_IgE_clean)
merging_table_3<- data.frame(
  cluster = as.character(1:14),
  merged_cluster = c(
    "Classical, resting", # 1 HLADR+ PDL1+
    "Classical, resting",  # 2 HLADR+
    "Non-classical, resting", # 3 CD40+ CD80+ CD86+
    "Classical, activated/antigen-presenting",   # 4 PDL1+
    "Intermediate, activated/antigen-presenting",   # 5 CD23+ PDL1+
    "Non-classical, activated/antigen-presenting", # 6 CD163+ CD206+ PDL1+
    "Non-classical, activated/antigen-presenting",  # 7 HLADR+
    "Classical, exhausted",  # 8 HLADR+
    "Classical, inhibited",   # 9 CD163+ CD206+ 
    "Classical, exhausted", # 10 CD163+ CD206+ PDL1+
    "Classical, inhibited", # 11 nothing
    "Classical, exhausted", # 12
    "Classical, inhibited",# 13
    "Classical, exhausted" # 14
  ))
numeric_order <- as.character(1:16)
merging_table_2$merged_cluster <- factor(
  merging_table_2$merged_cluster, 
  levels = numeric_order
)
sce_IgE_clean <- mergeClusters(sce_IgE_clean, k = "cluster_ID", table = merging_table_3, id = "Phenotype", overwrite = TRUE)

png("IgE_clean_exprs_heatmap_by_clusters.png", width = 12, height = 7, units = "in", res = 300)
plotExprHeatmap(sce_IgE_clean, features = ordered_markers, k = "cluster_ID",by = "cluster_id", fun = "median", row_anno = TRUE,
                scale = "last", q = 0, bars = TRUE, row_clust = FALSE, col_clust = FALSE, perc = TRUE, k_pal = cluster_colors)
dev.off()

#UMAP
sce_IgE_clean <- runDR(sce_IgE_clean, dr = "UMAP", features = "state")
plotDR(sce_IgE_clean, dr = "UMAP", color_by = "cluster_ID")

#tSNE
sce_IgE <- runDR(sce_IgE, dr = "TSNE", features = "state")
sce_IgE_clean <- runDR(sce_IgE_clean, cell = 1000, dr = "TSNE", features = "state", perplexity= 25, seed = 10)
png("IgE_clean_tSNE_all.png", width = 8, height = 7, units = "in", res = 300)
plotDR(sce_IgE_clean,dr = "TSNE", color_by = "cluster_ID")  +
  scale_color_manual(values = cluster_colors) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

# tSNE_each_condition
conditions <- unique(colData(sce_IgE_clean)$condition)
dir.create("tSNE_IgE_clean_by_each_condition", showWarnings = FALSE)
tSNE_plots <- lapply(conditions, function(m) {
  m_safe <- gsub(" ", "_", m)
  sce_sub <- filterSCE(sce_IgE_clean,k = "cluster_ID", condition == m)
  file_path <- file.path("tSNE_IgE_clean_by_each_condition", paste0("tSNE_", m_safe, ".png"))
  png(file_path, width = 8, height = 7, units = "in", res = 300)
  p <- plotDR(sce_sub, dr = "TSNE", color_by = "cluster_ID") + 
    scale_color_manual(values = cluster_colors) +
    ggtitle(paste(m)) + 
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
          axis.title = element_text(size = 20))
  print(p)
  dev.off()
  rm(p, sce_sub)
  gc() 
  return(NULL)
})

# tSNE_setting
png("tSNE_IgE_clean_3D.png", width = 8, height = 7, units = "in", res = 300)
plotDR(filterSCE(sce_IgE_clean, k = "cluster_ID", setting %in% c("3D")), dr = "TSNE", color_by = "setting") + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20)) + 
  scale_colour_manual(values = condition_colors)
dev.off()

# tSNE_sample type
png("tSNE_IgE_clean_sample.png", width = 8, height = 7, units = "in", res = 300)
plotDR(filterSCE(sce_IgE_clean, sample_type %in% c("HV", "M")), dr = "TSNE", color_by = "sample_type") + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20)) + 
  scale_colour_manual(values = condition_colors)
dev.off()

# Abundance by condition
png("cluster_abundance_sample_id_IgE_clean.png", width = 10, height = 5, units = "in", res = 300)
plotAbundances(sce_IgE_clean, k = "cluster_ID", by = "sample_id", group_by = "condition") + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    plot.title = element_text(face = "bold", size = 28)
  )
dev.off()

# Abundance by setting
png("cluster_abundance_setting_IgE_clean.png", width = 10, height = 7, units = "in", res = 300)
plotAbundances(sce_IgE_clean, k = "cluster_ID", by = "sample_id", group_by = "setting") + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    plot.title = element_text(face = "bold", size = 28)
  )
dev.off()

# Abundance by sample type
png("cluster_abundance_sample_IgE_clean.png", width = 10, height = 7, units = "in", res = 300)
plotAbundances(sce_IgE_clean, k = "cluster_ID", by = "sample_id", group_by = "sample_type") + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    plot.title = element_text(face = "bold", size = 28)
  )
dev.off()

# Pseudobulk expression boxplot
plotPbExprs(sce_IgE_clean, k = "cluster_ID", features = markers_ex_CD14, 
            group_by = "cluster_id", color_by = "condition", 
            size_by = FALSE, ncol = 5)
#=========================================================================================

# Differential analysis NIP IgE vs CSPG4 IgE
# Create a design matrix. This tells the model to look at 'condition' while accounting for 'patient_id' (random effect)
design <- createDesignMatrix(metadata(sce_IgE_clean)$experiment_info, 
                             cols_design = c("condition", "setting", "sample_type"))
# Define the contrast. This tells R exactly what to compare (e.g., CSPG4 IgE vs PBS)
contrast_NIP_vs_IgE <- createContrast(c(0, -1, 1, 0, 0))
# Run DA analysis using GLMM (Generalized Linear Mixed Models)
da_res_NIP_vs_IgE <- diffcyt(sce_IgE_clean, 
                             design = design, 
                             contrast = contrast_NIP_vs_IgE,
                             analysis_type = "DA", 
                             method_DA = "diffcyt-DA-edgeR", 
                             clustering_to_use = "cluster_ID")

topTable(da_res_NIP_vs_IgE, format_vals = TRUE)
#Differential Abundance Heatmap
# reorder SCE
new_sample_order <- colData(sce_IgE_clean) %>%
  as.data.frame() %>%
  arrange(condition, sample_id) %>% 
  pull(sample_id) %>%
  unique() %>%
  as.character()
sce_IgE_clean$sample_id <- factor(sce_IgE_clean$sample_id, levels = new_sample_order)
columns_to_show <- c("condition", "sample_type", "setting")

ds_res_NIP_vs_IgE <- diffcyt(sce_IgE_clean, 
                             design = design, 
                             contrast = contrast_NIP_vs_IgE,
                             analysis_type = "DS", 
                             method_DS = "diffcyt-DS-limma",
                             clustering_to_use = "cluster_ID")

topTable(ds_res_NIP_vs_IgE, format_vals = TRUE)

# Differential analysis PBS vs CSPG4 IgE
# Create a design matrix. This tells the model to look at 'condition' while accounting for 'patient_id' (random effect)
design <- createDesignMatrix(metadata(sce_IgE_clean)$experiment_info, 
                             cols_design = c("condition", "setting", "sample_type"))
# Define the contrast. This tells R exactly what to compare (e.g., CSPG4 IgE vs PBS)
contrast_PBS_vs_IgE <- createContrast(c(0, 0, 1, 0, 0)) # Use colnames(design) to see the exact names to use here

# Run DA analysis using GLMM (Generalized Linear Mixed Models)
da_res_PBS_vs_IgE <- diffcyt(sce_IgE_clean, 
                             design = design, 
                             contrast = contrast_PBS_vs_IgE,
                             analysis_type = "DA", 
                             method_DA = "diffcyt-DA-edgeR", 
                             clustering_to_use = "cluster_ID")

topTable(da_res_PBS_vs_IgE, format_vals = TRUE)
#Differential Abundance Heatmap
# reorder SCE
new_sample_order <- colData(sce_IgE_clean) %>%
  as.data.frame() %>%
  arrange(condition, sample_id) %>% 
  pull(sample_id) %>%
  unique() %>%
  as.character()
sce_IgE_clean$sample_id <- factor(sce_IgE_clean$sample_id, levels = new_sample_order)
columns_to_show <- c("condition", "sample_type", "setting")

ds_res_PBS_vs_IgE <- diffcyt(sce_IgE_clean, 
                             design = design, 
                             contrast = contrast_PBS_vs_IgE,
                             analysis_type = "DS", 
                             method_DS = "diffcyt-DS-limma",
                             clustering_to_use = "cluster_ID")

topTable(ds_res_PBS_vs_IgE, format_vals = TRUE)

#=================================================================================
# ADD P value
# cluster_abundance_condition_IgE

p <- plotAbundances(sce_IgE_clean, k = "cluster_ID", by = "cluster_id", shape_by = "setting") + 
  facet_wrap(~cluster_id, scales = "free_y", ncol = 7) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = alpha(condition_colors, 0.6)) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))

p$data <- p$data %>% 
  dplyr::filter(cluster_id %in% c(11, 12, 13, 14))
plot_data <- p$data

# Function to clean and format stats
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  if (!"cluster_id" %in% colnames(res_table)) {
    res_table <- res_table %>% rownames_to_column("cluster_id")
  }
  res_table %>%
    select(cluster_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}

# Process all 4 comparisons
s2 <- prep_stats(da_res_NIP_vs_IgE,  "2") # NIP_IgE(3) vs CSPG4_IgE(5)
s4 <- prep_stats(da_res_PBS_vs_IgE,  "4") # PBS(1) vs CSPG4_IgE(5)

label_data_final <- p$data %>%
  group_by(cluster_id) %>%
  summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE)) %>% 
  #left_join(s1, by = "cluster_id") %>%
  left_join(s2, by = "cluster_id") %>%
  #left_join(s3, by = "cluster_id") %>%
  left_join(s4, by = "cluster_id") %>%
  mutate(cluster_id = factor(cluster_id, levels = levels(p$data$cluster_id))) %>%
  mutate(
    y4 = local_max * 1.15, y4_t = local_max * 1.17, # Tier 2
    y2 = local_max * 1.05, y2_t = local_max * 1.07  # Tier 4
  )

p_final <- p + 
  facet_wrap(~cluster_id, scales = "free_y", ncol = 5) +
  
  # Tier 4: PBS (3) vs CSPG4 IgE (5)
  geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
               aes(x = 1, xend = 3, y = y4, yend = y4), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
            aes(x = 2, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 4, vjust = 0) +
  
  # Tier 2: NIP IgE (1) vs CSPG4_IgE (5)
  geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
               aes(x = 2, xend = 3, y = y2, yend = y2), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
            aes(x = 2.5, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 4, vjust = 0) +
  theme(panel.spacing.x = unit(1, "cm"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        text = element_text(size = 24),
        plot.title = element_text(face = "bold", size = 24))

p_final
png("cluster_abundance_exhausted_IgE_padj.png", width = 13, height = 4, units = "in", res = 300)
p_final
dev.off()

#=================================================================================
# p value of marker expression in IgE in specific cluster
# Function to clean and format stats
markers_ex_FcyR <- c("FCER1", "CD23", "CD80", "CD86", "CD40", "CCR2", "HLADR", "CD163", "CD206", "PDL1")
markers_ex_CD14 <- c("CD64", "CD32B","CD16", "FCER1", "CD23", "CD80", "CD86", "CD40", "CCR2", "HLADR", "CD163", "CD206", "PDL1")
n_cond <- 3
d_width <- 0.9
step <- d_width / n_cond
offsets <- (1:n_cond - (n_cond + 1) / 2) * step

p <- plotPbExprs(sce_IgE_clean, k = "cluster_ID", 
                 features = markers_ex_FcyR, facet_by = "cluster_id", shape_by = "setting", ncol = 4) + 
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

p$data <- p$data %>% 
  dplyr::filter(cluster_id %in% c(1))
plot_data <- p$data


prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  res_table %>%
    select(cluster_id, marker_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}


s2 <- prep_stats(ds_res_NIP_vs_IgE,  "2") # NIP_IgE(3) vs CSPG4_IgE(5)
s4 <- prep_stats(ds_res_PBS_vs_IgE,  "4") # PBS(1) vs CSPG4_IgE(5)

label_data_final <- p$data %>%
  group_by(antigen, cluster_id) %>%
  summarize(local_max = max(as.numeric(as.character(value)), na.rm = TRUE)) %>% 
  left_join(s2, by = c("antigen" = "marker_id", "cluster_id")) %>%
  left_join(s4, by = c("antigen" = "marker_id", "cluster_id")) %>%
  mutate(antigen = factor(antigen, levels = levels(p$data$antigen))) %>%
  mutate(
    y2 = local_max * 1.05, y2_t = local_max * 1.07, # Tier 1
    y4 = local_max * 1.15, y4_t = local_max * 1.17, # Tier 3
  )

#label_data_final <- label_data_final %>% dplyr::filter(p_label_2 != "" | p_label_4 != "")
#sig_antigens <- label_data_final$antigen
#p$data <- p$data %>% dplyr::filter(antigen %in% sig_antigens)

marker_levels <- levels(factor(p$data$antigen))
label_data_final <- label_data_final %>%
  mutate(x_pos = as.numeric(factor(antigen, levels = marker_levels)),
         # s2: NIP_IgE (2) vs CSPG4_IgE (3)
         x_s2 = x_pos + offsets[2], x_e2 = x_pos + offsets[3],
         # s4: PBS (1) vs CSPG4_IgE (3)
         x_s4 = x_pos + offsets[1], x_e4 = x_pos + offsets[3])

p$layers[[1]]$position <- position_dodge(width = 0.9)
p$layers[[1]]$geom_params$width <- 0.6 
p$layers[[2]]$position <- position_dodge(width = 0.9)

p_final <- p + 
  facet_wrap(~cluster_id, scales = "free_y", ncol = 4) +
  
  # Tier 4: PBS (1) vs CSPG4_IgE (4)
  geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
               aes(x = x_s4, xend = x_e4, y = y4, yend = y4), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
            aes(x = (x_s4 + x_e4)/2, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 4, vjust = 0) +
  # Tier 2: NIP_IgG1 (2) vs CSPG4_IgG1 (4)
  geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
               aes(x = x_s2, xend = x_e2, y = y2, yend = y2), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
            aes(x = (x_s2 + x_e2)/2, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 4, vjust = 0)

p_final
png("pbExprs_cluster_IgE_sig_padj.png", width = 5, height = 4, units = "in", res = 300)
p_final
dev.off()


#=================================================================================
# loop for marker expression in IgE in multiple clusters
# 1. Define the clusters you want to loop through
# You can use unique(sce_IgE_clean$cluster_ID) for all clusters
target_clusters <- as.character(1:14)
dir.create("pbExprs_cluster_IgE_clean", showWarnings = FALSE)
# 2. Pre-prepare the stats tables once (outside the loop)
s2 <- prep_stats(ds_res_NIP_vs_IgE, "2")
s4 <- prep_stats(ds_res_PBS_vs_IgE, "4")

# 3. Start the loop
lapply(target_clusters, function(id) {
  
  message(paste("Generating plot for cluster:", id))
  
  # Filter SCE for the specific cluster to keep plot data clean
  sce_sub <- filterSCE(sce_IgE_clean, cluster_ID == id)

  
  # Generate base plot for this cluster
  p <- plotPbExprs(sce_sub, k = "cluster_ID", 
                   features = markers_ex_CD14, 
                   facet_by = "cluster_id", 
                   shape_by = "setting", ncol = 4) + 
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    theme(text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0,    
                                     hjust = 0.5,  
                                     vjust = 0.5)
          )
  p$data <- p$data %>% 
    dplyr::filter(cluster_id == id)
  # Calculate dynamic Y offsets based on THIS cluster's data range
  vals <- as.numeric(as.character(p$data$value))
  y_max <- max(vals, na.rm = TRUE)
  y_min <- min(vals, na.rm = TRUE)
  y_range <- y_max - y_min
  jump <- y_range * 0.12
  
  label_data_final <- p$data %>%
    mutate(value = as.numeric(as.character(value))) %>%
    group_by(antigen, cluster_id) %>%
    summarize(local_max = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(s2, by = c("antigen" = "marker_id", "cluster_id")) %>%
    left_join(s4, by = c("antigen" = "marker_id", "cluster_id")) %>%
    mutate(
      y2 = local_max + (jump * 0.5),   # Bar 1
      y2_t = local_max + (jump * 0.7), # Text 1
      y4 = local_max + (jump * 1.5),   # Bar 2
      y4_t = local_max + (jump * 1.7)  # Text 2
    )
  
  # Adjust geom positions for dodging
  p$layers[[1]]$position <- position_dodge(width = 0.9)
  p$layers[[1]]$geom_params$width <- 0.6 
  p$layers[[2]]$position <- position_dodge(width = 0.9)

  
  # Calculate x-positions based on the antigen factor levels in the plot
  marker_levels <- levels(factor(p$data$antigen))
  label_data_final <- label_data_final %>%
    mutate(
      x_pos = as.numeric(factor(antigen, levels = marker_levels)),
      # s2: NIP (2) vs IgE (3)
      x_s2 = x_pos + offsets[2], 
      x_e2 = x_pos + offsets[3],
      # s4: PBS (1) vs IgE (3)
      x_s4 = x_pos + offsets[1], 
      x_e4 = x_pos + offsets[3])
  
  # Add layers for Tier 2 and Tier 4 stats
  p_final <- p + 
    geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
                 aes(x = x_s4, xend = x_e4, y = y4, yend = y4), inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
              aes(x = (x_s4 + x_e4)/2, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 3.5, vjust = 0) +
    geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
                 aes(x = x_s2, xend = x_e2, y = y2, yend = y2), inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
              aes(x = (x_s2 + x_e2)/2, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 3.5, vjust = 0) +
    coord_cartesian(ylim = c(y_min, y_max + (jump * 2.5))) # Add headroom at the top
  
  # Save the plot with a unique name for each cluster
  file_name <- paste0("pbExprs_cluster_IgE_clean/pbExprs_cluster_", id, "_sig_padj.png")
  png(file_name, width = 15, height = 5, units = "in", res = 300)
  print(p_final)
  dev.off()

})

#=================================================================================
#global p value of marker expression by condition
merging_table_global<- data.frame(
  cluster = as.character(1:16),
  merged_cluster = c(
    "1",
    "1",
    "1",
    "1", 
    "1",
    "1", 
    "1", 
    "1",  
    "1", 
    "1",
    "1",
    "1",
    "1", 
    "1",
    "1",
    "1"
  ))

sce_IgE_clean <- mergeClusters(sce_IgE_clean, k = "meta16", table = merging_table_global, id = "global", overwrite = TRUE)


ds_res_global_PBS_vs_IgE <- diffcyt(sce_IgE_clean, 
                                    design = design, 
                                    contrast = contrast_PBS_vs_IgE,
                                    analysis_type = "DS", 
                                    method_DS = "diffcyt-DS-limma",
                                    clustering_to_use = "global")
topTable(ds_res_global_PBS_vs_IgE)

ds_res_global_NIP_vs_IgE <- diffcyt(sce_IgE_clean, 
                                    design = design, 
                                    contrast = contrast_NIP_vs_IgE,
                                    analysis_type = "DS", 
                                    method_DS = "diffcyt-DS-limma",
                                    clustering_to_use = "global")
topTable(ds_res_global_NIP_vs_IgE)


p <- plotPbExprs(sce_IgE_clean, features = Inhibitory, facet_by = "antigen", shape_by = "setting", ncol = 7) + 
  scale_color_manual(values = condition_colors) + 
  scale_fill_manual(values = condition_colors) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())


plot_data <- p$data

# Function to clean and format stats
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  res_table %>%
    select(marker_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}

# Process all 4 comparisons
s2 <- prep_stats(ds_res_global_NIP_vs_IgE,  "2") # NIP_IgE(3) vs CSPG4_IgE(5)
s4 <- prep_stats(ds_res_global_PBS_vs_IgE,  "4") # PBS(1) vs CSPG4_IgE(5)

label_data_final <- p$data %>%
  group_by(antigen) %>%
  summarize(local_max = max(as.numeric(as.character(value)), na.rm = TRUE)) %>% 
  #left_join(s1, by = c("antigen" = "marker_id")) %>%
  left_join(s2, by = c("antigen" = "marker_id")) %>%
  #left_join(s3, by = c("antigen" = "marker_id")) %>%
  left_join(s4, by = c("antigen" = "marker_id")) %>%
  mutate(antigen = factor(antigen, levels = levels(p$data$antigen))) %>%
  mutate(
    #y3 = local_max * 1.10, y3_t = local_max * 1.12, # Tier 1
    y4 = local_max * 1.10, y4_t = local_max * 1.12, # Tier 2
    #y1 = local_max * 1.30, y1_t = local_max * 1.32, # Tier 3
    y2 = local_max * 1.20, y2_t = local_max * 1.22  # Tier 4
  )

p_final <- p + 
  facet_wrap(~antigen, scales = "free_y", ncol = 7) +
  
  # Tier 2: PBS (1) vs CSPG4_IgE (5)
  geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
               aes(x = 1, xend = 3, y = y2, yend = y2), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
            aes(x = 2, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 3.5, vjust = 0) +
  
  # Tier 4: NIP IgE (3) vs CSPG4 IgE (5)
  geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
               aes(x = 2, xend = 3, y = y4, yend = y4), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
            aes(x = 2.5, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 3.5, vjust = 0)
p_final
png("pbExprs_IgE_Inhibitory_padj.png", width = 9, height = 4, units = "in", res = 300)
p_final
dev.off()
#=================================================================
# Density plot
library(ggh4x)

plotClusterExprs(sce_IgE_clean, k = "cluster_ID", features = markers_ex_CD14)

plot_df <- as.data.frame(t(assay(sce_IgE_clean, "exprs")[markers_ex_CD14, ]))
cl_ids <- cluster_ids(sce_IgE_clean, "cluster_ID")
plot_df$cluster_id <- factor(as.character(cl_ids), levels = rev(as.character(1:14)))
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(
    cols = -cluster_id, 
    names_to = "antigen", 
    values_to = "expression")
plot_df_long$antigen <- factor(plot_df_long$antigen, levels = as.character(markers_ex_CD14))

groups <- list(
  "FcRs" = c("CD64", "CD32B", "CD16", "FCER1", "CD23"),
  "Activation/Antigen-presenting" = c("CD80", "CD86", "CD40", "CCR2", "HLADR"),
  "Inhibition/Exhaustion" = c("CD163", "CD206", "PDL1", " ", "  ")
)
medians_cl1 <- plot_df_long %>%
  filter(cluster_id == "1") %>%
  group_by(antigen) %>%
  summarize(m_val = median(expression, na.rm = TRUE))


facet_themes <- list(
  " " = theme(panel.grid.major = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_blank()),
  "  " = theme(panel.grid.major = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_blank())
)

make_ridge_row <- function(data, markers, title_text, show_x = FALSE) {
  df_sub <- data %>% 
    filter(antigen %in% markers, cluster_id %in% c("11", "12", "13", "14")) %>%
    mutate(antigen = factor(antigen, levels = markers))  
  meds_sub <- medians_cl1 %>% 
    filter(antigen %in% markers) %>%
    mutate(antigen = factor(as.character(antigen), levels = markers))
p <- ggplot(df_sub, aes(x = expression, y = cluster_id, fill = cluster_id, color = cluster_id)) +
  # Add vertical lines first so they sit BEHIND the ridges
  geom_vline(data = meds_sub, aes(xintercept = m_val), 
             color = "grey20", linetype = "dashed", linewidth = 0.5) +
  # The ridges
  geom_density_ridges(alpha = 0.7, scale = 1.5, linewidth = 0.5) +
  # Apply your cluster_colour vector
  scale_fill_manual(values = cluster_colors) + 
  scale_color_manual(values = cluster_colors) +
  # Separate by antigen group
  facet_wrap(~antigen, nrow = 1, ncol = 5, scales = "free_x", drop = FALSE) +
  facetted_pos_scales(x = list(antigen %in% c(" ", "  ") ~ scale_x_continuous(guide = "none"))) +
  labs(subtitle = title_text, x = NULL, y = "Cluster ID") +
  # Styling to match the original
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
  theme(
    plot.title = element_text(size = 20),
    plot.subtitle = element_text(size = 14),
    strip.background = element_blank(), # Removes grey box behind antigen names
    panel.grid.major = element_line(color = "grey90"), # Keeps the background grid
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 14, margin = margin(b = 5)), # Makes antigen names bold
    panel.background = element_blank(),      # Ensures clean background
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  )
p <- p + theme(panel.grid.major.y = element_blank())
if (!show_x) { 
  p <- p + theme(axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank()) 
}

return(p)
}

p1 <- make_ridge_row(plot_df_long, groups[[1]], names(groups)[1])
p2 <- make_ridge_row(plot_df_long, groups[[2]], names(groups)[2])
p3 <- make_ridge_row(plot_df_long, groups[[3]], names(groups)[3], show_x = TRUE)

final_plot <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
final_plot
png("densities_IgE_exhausted_padj.png", width = 9, height = 6, units = "in", res = 300)
final_plot
dev.off()

#==========================================================================
# Psedotime analysis with whole population
library(slingshot)
library(SingleCellExperiment)

cells_with_tsne <- !is.na(reducedDim(sce_IgE_clean, "TSNE")[,1])
sce_IgE_tsne <- sce_IgE_clean[, cells_with_tsne]
sce_IgE_tsne$merged_ids <- cluster_ids(sce_IgE_tsne, "cluster_ID")

sce_IgE_tsne <- slingshot(sce_IgE_tsne, 
                     clusterLabels = "merged_ids", 
                     reducedDim = 'TSNE', 
                     start.clus = c("1", "2", "3"),
                     reweight = FALSE)
SlingshotDataSet(sce_IgE_tsne)

plot(reducedDim(sce_IgE_tsne, "TSNE"), 
     pch = 16, cex = 0.5,
     main = "Slingshot Trajectory on tSNE")
lines(SlingshotDataSet(sce_IgE_tsne), lwd = 3, col = 'black')

df_plot <- as.data.frame(reducedDim(sce_IgE_tsne, "TSNE"))
colnames(df_plot) <- c("tSNE1", "tSNE2")
merged_labels <- cluster_ids(sce_IgE_tsne, "cluster_ID")
df_plot$Cluster <- merged_labels[!is.na(reducedDim(sce_IgE_tsne, "TSNE")[,1])]

# plotting curve on tSNE plot
sds <- SlingshotDataSet(sce_IgE_tsne)
curve_list <- lapply(seq_along(slingCurves(sds)), function(i) {
  curr_curve <- slingCurves(sds)[[i]]
  df <- as.data.frame(curr_curve$s[curr_curve$ord, ])
  colnames(df) <- c("tSNE1", "tSNE2")
  df$Lineage <- as.character(i)
  return(df)
})
all_curves <- do.call(rbind, curve_list)

ggplot(df_plot, aes(x = tSNE1, y = tSNE2, color = Cluster)) +
  geom_point(size = 0.5) +
  geom_path(data = all_curves, aes(group = Lineage), 
            color = "black", size = 1, arrow = arrow(length = unit(0.2, "cm"))) +
  theme_minimal() +
  scale_color_manual(values = cluster_colors) +
  labs(title = "Monocyte Clusters on tSNE",
       subtitle = "Trajectory Path: Resting -> Activation -> Exhaustion -> Inhibition") +
  guides(color = guide_legend(override.aes = list(size = 3)))

# Plotting mean spinning tree (MST) on tSNE plot
#Calculate the mean tSNE position for every cluster
centers_manual <- df_plot %>%
  group_by(Cluster) %>%
  summarize(x = mean(tSNE1), y = mean(tSNE2)) %>%
  as.data.frame()

# Set the row names to the Cluster ID so we can look them up
rownames(centers_manual) <- centers_manual$Cluster

# Define your exact lineages again (from your previous output)
paths <- list(
  c("1", "4", "11", "13", "10", "14", "8", "9"),
  c("1", "4", "11", "13", "10", "6", "7", "3"),
  c("1", "4", "11", "13", "10", "5"),
  c("1", "4", "11", "12"),
  c("1", "4", "2")
)

# Build the edge_coords by connecting the manual centers
edge_coords <- data.frame()

for (p in paths) {
  for (i in 1:(length(p) - 1)) {
    start_node <- p[i]
    end_node   <- p[i+1]
    
    # Only add if both clusters exist in our manual centers
    if (start_node %in% centers_manual$Cluster & end_node %in% centers_manual$Cluster) {
      edge_coords <- rbind(edge_coords, data.frame(
        x = centers_manual[start_node, "x"],
        y = centers_manual[start_node, "y"],
        xend = centers_manual[end_node, "x"],
        yend = centers_manual[end_node, "y"]
      ))
    }
  }
}

edge_coords <- unique(edge_coords)
print(nrow(edge_coords))


p <- ggplot(df_plot, aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = Cluster), size = 0.5) +
  geom_segment(data = edge_coords, 
               aes(x = x, y = y, xend = xend, yend = yend),
               #arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  geom_label(data = centers_manual, 
             aes(x = x, y = y, label = Cluster),
             size = 4, fontface = "bold", alpha = 0.8, inherit.aes = FALSE) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14)) +
  scale_color_manual(values = cluster_colors) +
  labs(title = "Monocyte Clusters on tSNE",
       subtitle = "Trajectory Path: Resting -> Activation -> Exhaustion -> Inhibition")
p
png("tSNE_IgE_clean_trajectory.png", width = 8, height = 7, units = "in", res = 300)
p
dev.off()

df_plot$Pseudotime <- rowMeans(slingPseudotime(sce_IgE_tsne), na.rm = TRUE)

# Create the boxplot to check the "age" of each cluster
ggplot(df_plot, aes(x = reorder(Cluster, Pseudotime, FUN = median), y = Pseudotime, fill = Cluster)) +
  geom_boxplot(outlier.size = 0.1, alpha = 0.7) +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal() +
  labs(title = "Pseudotime Progression by Cluster",
       x = "Clusters (Ordered by 'Age')",
       y = "Developmental Pseudotime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Marker expression kinetics
plot_data <- as.data.frame(t(assay(sce_IgE_tsne, "exprs")))
plot_data$Pseudotime <- sce_IgE_tsne$pseudotime

# Melt for ggplot
plot_melted <- melt(plot_data, id.vars = "Pseudotime")

# Plot marker expression over pseudotime
ggplot(plot_melted, aes(x = Pseudotime, y = value, color = variable)) +
  geom_smooth(method = "gam", se = FALSE) + 
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Marker Expression Kinetics: Resting to Exhausted",
       y = "Arcsinh Expression", x = "Developmental Progress (Pseudotime)")

#Pseudotime analysis with calssical monocytes only
cells_with_tsne <- !is.na(reducedDim(sce_IgE_clean, "TSNE")[,1])
sce_IgE_tsne <- sce_IgE_clean[, cells_with_tsne]
sce_IgE_tsne$merged_ids <- cluster_ids(sce_IgE_tsne, "cluster_ID")
keep_clusters <- c("1", "2", "4", "8", "9", "11", "12", "13")
sce_classical_tsne <- sce_IgE_tsne[, sce_IgE_tsne$merged_ids %in% keep_clusters]

sce_classical_tsne <- slingshot(sce_classical_tsne, 
                          clusterLabels = "merged_ids", 
                          reducedDim = 'TSNE', 
                          start.clus = c("1"),
                          reweight = FALSE)
SlingshotDataSet(sce_classical_tsne)

df_plot <- as.data.frame(reducedDim(sce_classical_tsne, "TSNE"))
colnames(df_plot) <- c("tSNE1", "tSNE2")
merged_labels <- cluster_ids(sce_classical_tsne, "cluster_ID")
df_plot$Cluster <- merged_labels[!is.na(reducedDim(sce_classical_tsne, "TSNE")[,1])]
# Plotting mean spinning tree (MST) on tSNE plot
#Calculate the mean tSNE position for every cluster
centers_manual <- df_plot %>%
  group_by(Cluster) %>%
  summarize(x = mean(tSNE1), y = mean(tSNE2)) %>%
  as.data.frame()

# Set the row names to the Cluster ID so we can look them up
rownames(centers_manual) <- centers_manual$Cluster

# Define your exact lineages again (from your previous output)
paths <- list(
  c("1", "4", "11", "12", "9", "8"),
  c("1", "4", "11", "13"),
  c("1", "4", "2")
)

# Build the edge_coords by connecting the manual centers
edge_coords <- data.frame()

for (p in paths) {
  for (i in 1:(length(p) - 1)) {
    start_node <- p[i]
    end_node   <- p[i+1]
    
    # Only add if both clusters exist in our manual centers
    if (start_node %in% centers_manual$Cluster & end_node %in% centers_manual$Cluster) {
      edge_coords <- rbind(edge_coords, data.frame(
        x = centers_manual[start_node, "x"],
        y = centers_manual[start_node, "y"],
        xend = centers_manual[end_node, "x"],
        yend = centers_manual[end_node, "y"]
      ))
    }
  }
}

edge_coords <- unique(edge_coords)
print(nrow(edge_coords))


p <- ggplot(df_plot, aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = Cluster), size = 0.5) +
  geom_segment(data = edge_coords, 
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  geom_label(data = centers_manual, 
             aes(x = x, y = y, label = Cluster),
             nudge_x = 3, nudge_y = 3,
             size = 4, fontface = "bold", alpha = 0.8, inherit.aes = FALSE) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = cluster_colors) +
  labs(title = "Classical Monocyte Clusters on tSNE",
       subtitle = "Trajectory Path: Resting -> Activation -> Exhaustion -> Inhibition")
p
png("tSNE_IgE_clean_classical_trajectory.png", width = 8, height = 7, units = "in", res = 300)
p
dev.off()

#======================================================================================
#======================================================================================
# Analysis with 2D data
sce_IgE_clean_2D <- filterSCE(sce_IgE_clean, setting == "2D")
png("IgE_clean_2D_tSNE_all.png", width = 8, height = 7, units = "in", res = 300)
plotDR(sce_IgE_clean_2D,dr = "TSNE", color_by = "cluster_ID")  +
  scale_color_manual(values = cluster_colors) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

# Abundance by condition
png("cluster_abundance_sample_id_IgE_clean_2D.png", width = 10, height = 5, units = "in", res = 300)
plotAbundances(sce_IgE_clean_2D, k = "cluster_ID", by = "sample_id", group_by = "condition") + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    plot.title = element_text(face = "bold", size = 28)
  )
dev.off()

# Abundance by sample type
png("cluster_abundance_sample_IgE_clean_2D.png", width = 8, height = 5, units = "in", res = 300)
plotAbundances(sce_IgE_clean_2D, k = "cluster_ID", by = "sample_id", group_by = "sample_type") + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    plot.title = element_text(face = "bold", size = 28)
  )
dev.off()
#======================================================================================
# Differential analysis NIP IgE vs CSPG4 IgE
design <- createDesignMatrix(metadata(sce_IgE_clean_2D)$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_NIP_vs_IgE_2D <- createContrast(c(0, -1, 1, 0))
da_res_NIP_vs_IgE_2D <- diffcyt(sce_IgE_clean_2D, 
                             design = design, 
                             contrast = contrast_NIP_vs_IgE_2D,
                             analysis_type = "DA", 
                             method_DA = "diffcyt-DA-edgeR", 
                             clustering_to_use = "cluster_ID")

topTable(da_res_NIP_vs_IgE_2D, format_vals = TRUE)

ds_res_NIP_vs_IgE_2D <- diffcyt(sce_IgE_clean_2D, 
                             design = design, 
                             contrast = contrast_NIP_vs_IgE_2D,
                             analysis_type = "DS", 
                             method_DS = "diffcyt-DS-limma",
                             clustering_to_use = "cluster_ID")

topTable(ds_res_NIP_vs_IgE_2D, format_vals = TRUE)

# Differential analysis PBS vs CSPG4 IgE
design <- createDesignMatrix(metadata(sce_IgE_clean_2D)$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_PBS_vs_IgE_2D <- createContrast(c(0, 0, 1, 0)) 
da_res_PBS_vs_IgE_2D <- diffcyt(sce_IgE_clean_2D, 
                             design = design, 
                             contrast = contrast_PBS_vs_IgE_2D,
                             analysis_type = "DA", 
                             method_DA = "diffcyt-DA-edgeR", 
                             clustering_to_use = "cluster_ID")

topTable(da_res_PBS_vs_IgE_2D, format_vals = TRUE)

ds_res_PBS_vs_IgE_2D <- diffcyt(sce_IgE_clean_2D, 
                             design = design, 
                             contrast = contrast_PBS_vs_IgE_2D,
                             analysis_type = "DS", 
                             method_DS = "diffcyt-DS-limma",
                             clustering_to_use = "cluster_ID")

topTable(ds_res_PBS_vs_IgE_2D, format_vals = TRUE)

# Differential analysis PBS HV vs M
design <- createDesignMatrix(metadata(filterSCE(sce_IgE_clean_2D, condition = "PBS"))$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_PBS_HV_vs_M_2D <- createContrast(c(0, 0, 0, 1)) 
da_res_PBS_HV_vs_M_2D <- diffcyt(filterSCE(sce_IgE_clean_2D, condition = "PBS"), 
                                design = design, 
                                contrast = contrast_PBS_HV_vs_M_2D,
                                analysis_type = "DA", 
                                method_DA = "diffcyt-DA-edgeR", 
                                clustering_to_use = "cluster_ID")

topTable(da_res_PBS_HV_vs_M_2D, format_vals = TRUE)

ds_res_PBS_HV_vs_M_2D <- diffcyt(filterSCE(sce_IgE_clean_2D, condition = "PBS"), 
                                design = design, 
                                contrast = contrast_PBS_HV_vs_M_2D,
                                analysis_type = "DS", 
                                method_DS = "diffcyt-DS-limma",
                                clustering_to_use = "cluster_ID")

topTable(ds_res_PBS_HV_vs_M_2D, format_vals = TRUE)

# Differential analysis NIP HV vs M
design <- createDesignMatrix(metadata(filterSCE(sce_IgE_clean_2D, condition = "NIP"))$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_NIP_HV_vs_M_2D <- createContrast(c(0, 0, 0, 1)) 
da_res_NIP_HV_vs_M_2D <- diffcyt(filterSCE(sce_IgE_clean_2D, condition = "NIP"), 
                                 design = design, 
                                 contrast = contrast_NIP_HV_vs_M_2D,
                                 analysis_type = "DA", 
                                 method_DA = "diffcyt-DA-edgeR", 
                                 clustering_to_use = "cluster_ID")

topTable(da_res_NIP_HV_vs_M_2D, format_vals = TRUE)

ds_res_NIP_HV_vs_M_2D <- diffcyt(filterSCE(sce_IgE_clean_2D, condition = "NIP"), 
                                 design = design, 
                                 contrast = contrast_NIP_HV_vs_M_2D,
                                 analysis_type = "DS", 
                                 method_DS = "diffcyt-DS-limma",
                                 clustering_to_use = "cluster_ID")

topTable(ds_res_NIP_HV_vs_M_2D, format_vals = TRUE)

# Differential analysis CSPG4 HV vs M
design <- createDesignMatrix(metadata(filterSCE(sce_IgE_clean_2D, condition = "CSPG4"))$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_CSPG4_HV_vs_M_2D <- createContrast(c(0, 0, 0, 1)) 
da_res_CSPG4_HV_vs_M_2D <- diffcyt(filterSCE(sce_IgE_clean_2D, condition = "CSPG4"), 
                                 design = design, 
                                 contrast = contrast_CSPG4_HV_vs_M_2D,
                                 analysis_type = "DA", 
                                 method_DA = "diffcyt-DA-edgeR", 
                                 clustering_to_use = "cluster_ID")

topTable(da_res_CSPG4_HV_vs_M_2D, format_vals = TRUE)

ds_res_CSPG4_HV_vs_M_2D <- diffcyt(filterSCE(sce_IgE_clean_2D, condition = "CSPG4"), 
                                 design = design, 
                                 contrast = contrast_CSPG4_HV_vs_M_2D,
                                 analysis_type = "DS", 
                                 method_DS = "diffcyt-DS-limma",
                                 clustering_to_use = "cluster_ID")

topTable(ds_res_CSPG4_HV_vs_M_2D, format_vals = TRUE)
#=================================================================================
# ADD P value
# cluster_abundance_condition_IgE

p <- plotAbundances(sce_IgE_clean_2D, k = "cluster_ID", by = "cluster_id", shape_by = "sample_type") + 
  facet_wrap(~cluster_id, scales = "free_y", ncol = 7) +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = alpha(condition_colors, 0.6)) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))

p$data <- p$data %>% 
  dplyr::filter(cluster_id %in% c(11, 12, 13, 14))
plot_data <- p$data

# Function to clean and format stats
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  if (!"cluster_id" %in% colnames(res_table)) {
    res_table <- res_table %>% rownames_to_column("cluster_id")
  }
  res_table %>%
    select(cluster_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}

# Process all 4 comparisons
s2 <- prep_stats(da_res_NIP_vs_IgE_2D,  "2") # NIP_IgE(3) vs CSPG4_IgE(5)
s4 <- prep_stats(da_res_PBS_vs_IgE_2D,  "4") # PBS(1) vs CSPG4_IgE(5)

label_data_final <- p$data %>%
  group_by(cluster_id) %>%
  summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE)) %>% 
  left_join(s2, by = "cluster_id") %>%
  left_join(s4, by = "cluster_id") %>%
  mutate(cluster_id = factor(cluster_id, levels = levels(p$data$cluster_id))) %>%
  mutate(
    y4 = local_max * 1.15, y4_t = local_max * 1.17, # Tier 2
    y2 = local_max * 1.05, y2_t = local_max * 1.07  # Tier 4
  )

p_final <- p + 
  facet_wrap(~cluster_id, scales = "free_y", ncol = 5) + 
  geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
               aes(x = 1, xend = 3, y = y4, yend = y4), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
            aes(x = 2, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 4, vjust = 0) +
  geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
               aes(x = 2, xend = 3, y = y2, yend = y2), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
            aes(x = 2.5, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 4, vjust = 0) +
  theme(panel.spacing.x = unit(1, "cm"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        text = element_text(size = 24),
        plot.title = element_text(face = "bold", size = 24))

p_final
png("cluster_abundance_resting_IgE_2D_padj.png", width = 10, height = 4, units = "in", res = 300)
p_final
dev.off()

#=============================================================================================
#Differential abundance analysis with sub boxplot with P value
# 1. Generate the base object to get the data structure
p <- plotAbundances(sce_IgE_clean_2D, k = "cluster_ID", by = "cluster_id", shape_by = "sample_type")

plot_data <- p$data %>% 
  dplyr::filter(cluster_id %in% c(3, 7, 9))

# Function to clean and format stats
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  if (!"cluster_id" %in% colnames(res_table)) {
    res_table <- res_table %>% rownames_to_column("cluster_id")
  }
  res_table %>%
    select(cluster_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}

s1 <- prep_stats(da_res_PBS_HV_vs_M_2D,  "1")
s2 <- prep_stats(da_res_NIP_HV_vs_M_2D,  "2")
s3 <- prep_stats(da_res_CSPG4_HV_vs_M_2D,  "3")

label_data_final <- plot_data %>%
  group_by(cluster_id) %>%
  summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE)) %>% 
  left_join(s1, by = "cluster_id") %>%
  left_join(s2, by = "cluster_id") %>%
  left_join(s3, by = "cluster_id") %>%
  mutate(cluster_id = factor(cluster_id, levels = levels(plot_data$cluster_id))) %>%
  mutate(
    y1 = local_max * 1.05, y1_t = local_max * 1.07,
    y2 = local_max * 1.10, y2_t = local_max * 1.12,
    y3 = local_max * 1.15, y3_t = local_max * 1.17
  )

p_final <- ggplot(plot_data, aes(x = condition, y = Freq, color = condition)) +

  geom_boxplot(aes(fill = condition, alpha = sample_type),
               width = 0.6,
               position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               show.legend = FALSE) + 
  
  geom_point(aes(shape = sample_type, group = sample_type), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5, 
             stroke = 1, 
             fill = "transparent",
             show.legend = TRUE) + 

  facet_wrap(~cluster_id, scales = "free_y", ncol = 4) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  scale_alpha_manual(values = c("HV" = 0.1, "M" = 0.3)) + 
  scale_shape_manual(values = c("HV" = 21, "M" = 24)) +
  
  geom_segment(data = dplyr::filter(label_data_final, p_label_1 != ""),
               aes(x =0.8, xend = 1.2, y = y1, yend = y1), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_1 != ""),
            aes(x = 1, y = y1_t, label = p_label_1), 
            inherit.aes = FALSE, size = 5, color = "black") +
  
  geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
               aes(x = 1.8, xend = 2.2, y = y2, yend = y2), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
            aes(x = 2, y = y2_t, label = p_label_2), 
            inherit.aes = FALSE, size = 5, color = "black") +
  
  geom_segment(data = dplyr::filter(label_data_final, p_label_3 != ""),
               aes(x = 2.8, xend = 3.2, y = y3, yend = y3), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_3 != ""),
            aes(x = 3, y = y3_t, label = p_label_3), 
            inherit.aes = FALSE, size = 5, color = "black") +
  
  theme_bw() +
  labs(y = "Proportion [%]", x = NULL) +
  theme(
    text = element_text(size = 20),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold", size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "cm")
  )

p_final

#=================================================================================
# loop for marker expression in IgE in multiple clusters
target_clusters <- as.character(1:14)
dir.create("pbExprs_cluster_IgE_clean_2D", showWarnings = FALSE)
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  res_table %>%
    select(cluster_id, marker_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}
s2 <- prep_stats(ds_res_NIP_vs_IgE_2D, "2")
s4 <- prep_stats(ds_res_PBS_vs_IgE_2D, "4")

lapply(target_clusters, function(id) {
  message(paste("Generating plot for cluster:", id))
  sce_sub <- filterSCE(sce_IgE_clean_2D, cluster_ID == id)
  p <- plotPbExprs(sce_sub, k = "cluster_ID", 
                   features = markers_ex_CD14, 
                   facet_by = "cluster_id", 
                   shape_by = "sample_type", ncol = 4) + 
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    scale_shape_manual(values = c(1, 2)) +
    theme(text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0,    
                                     hjust = 0.5,  
                                     vjust = 0.5))
  p$data <- p$data %>% 
    dplyr::filter(cluster_id == id)
  vals <- as.numeric(as.character(p$data$value))
  y_max <- max(vals, na.rm = TRUE)
  y_min <- min(vals, na.rm = TRUE)
  y_range <- y_max - y_min
  jump <- y_range * 0.12
  
  label_data_final <- p$data %>%
    mutate(value = as.numeric(as.character(value))) %>%
    group_by(antigen, cluster_id) %>%
    summarize(local_max = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(s2, by = c("antigen" = "marker_id", "cluster_id")) %>%
    left_join(s4, by = c("antigen" = "marker_id", "cluster_id")) %>%
    mutate(
      y2 = local_max + (jump * 0.5),   # Bar 1
      y2_t = local_max + (jump * 0.7), # Text 1
      y4 = local_max + (jump * 1.5),   # Bar 2
      y4_t = local_max + (jump * 1.7)  # Text 2
    )
  
  p$layers[[1]]$position <- position_dodge(width = 0.9)
  p$layers[[1]]$geom_params$width <- 0.6 
  p$layers[[2]]$position <- position_dodge(width = 0.9)
  
  marker_levels <- levels(factor(p$data$antigen))
  label_data_final <- label_data_final %>%
    mutate(
      x_pos = as.numeric(factor(antigen, levels = marker_levels)),
      x_s2 = x_pos + offsets[2], 
      x_e2 = x_pos + offsets[3],
      x_s4 = x_pos + offsets[1], 
      x_e4 = x_pos + offsets[3])
  
  p_final <- p + 
    geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
                 aes(x = x_s4, xend = x_e4, y = y4, yend = y4), inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
              aes(x = (x_s4 + x_e4)/2, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 3.5, vjust = 0) +
    geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
                 aes(x = x_s2, xend = x_e2, y = y2, yend = y2), inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
              aes(x = (x_s2 + x_e2)/2, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 3.5, vjust = 0) +
    coord_cartesian(ylim = c(y_min, y_max + (jump * 2.5))) 
  
  file_name <- paste0("pbExprs_cluster_IgE_clean_2D/pbExprs_cluster_", id, "_sig_padj.png")
  png(file_name, width = 15, height = 5, units = "in", res = 300)
  print(p_final)
  dev.off()
  })
#=================================================================================
# Density plot
library(ggh4x)

plotClusterExprs(sce_IgE_clean_2D, k = "cluster_ID", features = markers_ex_CD14)

plot_df <- as.data.frame(t(assay(sce_IgE_clean_2D, "exprs")[markers_ex_CD14, ]))
cl_ids <- cluster_ids(sce_IgE_clean_2D, "cluster_ID")
plot_df$cluster_id <- factor(as.character(cl_ids), levels = rev(as.character(1:14)))
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(
    cols = -cluster_id, 
    names_to = "antigen", 
    values_to = "expression")
plot_df_long$antigen <- factor(plot_df_long$antigen, levels = as.character(markers_ex_CD14))

groups <- list(
  "FcRs" = c("CD64", "CD32B", "CD16", "FCER1", "CD23"),
  "Activation/Antigen-presenting" = c("CD80", "CD86", "CD40", "CCR2", "HLADR"),
  "Inhibition/Exhaustion" = c("CD163", "CD206", "PDL1", " ", "  ")
)
medians_cl1 <- plot_df_long %>%
  filter(cluster_id == "1") %>%
  group_by(antigen) %>%
  summarize(m_val = median(expression, na.rm = TRUE))


facet_themes <- list(
  " " = theme(panel.grid.major = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_blank()),
  "  " = theme(panel.grid.major = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_blank())
)

make_ridge_row <- function(data, markers, title_text, show_x = FALSE) {
  df_sub <- data %>% 
    filter(antigen %in% markers, cluster_id %in% c("1", "2", "3")) %>%
    mutate(antigen = factor(antigen, levels = markers))  
  meds_sub <- medians_cl1 %>% 
    filter(antigen %in% markers) %>%
    mutate(antigen = factor(as.character(antigen), levels = markers))
  p <- ggplot(df_sub, aes(x = expression, y = cluster_id, fill = cluster_id, color = cluster_id)) +
    # Add vertical lines first so they sit BEHIND the ridges
    geom_vline(data = meds_sub, aes(xintercept = m_val), 
               color = "grey20", linetype = "dashed", linewidth = 0.5) +
    # The ridges
    geom_density_ridges(alpha = 0.7, scale = 1.5, linewidth = 0.5) +
    # Apply your cluster_colour vector
    scale_fill_manual(values = cluster_colors) + 
    scale_color_manual(values = cluster_colors) +
    # Separate by antigen group
    facet_wrap(~antigen, nrow = 1, ncol = 5, scales = "free_x", drop = FALSE) +
    facetted_pos_scales(x = list(antigen %in% c(" ", "  ") ~ scale_x_continuous(guide = "none"))) +
    labs(subtitle = title_text, x = NULL, y = "Cluster ID") +
    # Styling to match the original
    theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
    theme(
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 14),
      strip.background = element_blank(), # Removes grey box behind antigen names
      panel.grid.major = element_line(color = "grey90"), # Keeps the background grid
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 14, margin = margin(b = 5)), # Makes antigen names bold
      panel.background = element_blank(),      # Ensures clean background
      legend.position = "none",
      panel.spacing = unit(1, "lines")
    )
  p <- p + theme(panel.grid.major.y = element_blank())
  if (!show_x) { 
    p <- p + theme(axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank()) 
  }
  return(p)
}
p1 <- make_ridge_row(plot_df_long, groups[[1]], names(groups)[1])
p2 <- make_ridge_row(plot_df_long, groups[[2]], names(groups)[2])
p3 <- make_ridge_row(plot_df_long, groups[[3]], names(groups)[3], show_x = TRUE)

final_plot <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
final_plot
png("densities_IgE_2D_resting_padj.png", width = 9, height = 6, units = "in", res = 300)
final_plot
dev.off()
#==========================================================================
#======================================================================================
#======================================================================================
# Analysis with 3D data
sce_IgE_clean_3D <- filterSCE(sce_IgE_clean, setting == "3D")
png("IgE_clean_3D_tSNE_all.png", width = 8, height = 7, units = "in", res = 300)
plotDR(sce_IgE_clean_3D,dr = "TSNE", color_by = "cluster_ID")  +
  scale_color_manual(values = cluster_colors) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()

# Abundance by condition
png("cluster_abundance_sample_id_IgE_clean_3D.png", width = 10, height = 5, units = "in", res = 300)
plotAbundances(sce_IgE_clean_3D, k = "cluster_ID", by = "sample_id", group_by = "condition") + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    plot.title = element_text(face = "bold", size = 28)
  )
dev.off()

# Abundance by sample type
png("cluster_abundance_sample_IgE_clean_3D.png", width = 8, height = 5, units = "in", res = 300)
plotAbundances(sce_IgE_clean_3D, k = "cluster_ID", by = "sample_id", group_by = "sample_type") + 
  scale_fill_manual(values = cluster_colors) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    text = element_text(size = 24),
    plot.title = element_text(face = "bold", size = 28)
  )
dev.off()
#======================================================================================
# Differential analysis NIP IgE vs CSPG4 IgE
design <- createDesignMatrix(metadata(sce_IgE_clean_3D)$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_NIP_vs_IgE_3D <- createContrast(c(0, -1, 1, 0))
da_res_NIP_vs_IgE_3D <- diffcyt(sce_IgE_clean_3D, 
                                design = design, 
                                contrast = contrast_NIP_vs_IgE_3D,
                                analysis_type = "DA", 
                                method_DA = "diffcyt-DA-edgeR", 
                                clustering_to_use = "cluster_ID")

topTable(da_res_NIP_vs_IgE_3D, format_vals = TRUE)

ds_res_NIP_vs_IgE_3D <- diffcyt(sce_IgE_clean_3D, 
                                design = design, 
                                contrast = contrast_NIP_vs_IgE_3D,
                                analysis_type = "DS", 
                                method_DS = "diffcyt-DS-limma",
                                clustering_to_use = "cluster_ID")

topTable(ds_res_NIP_vs_IgE_3D, format_vals = TRUE)

# Differential analysis PBS vs CSPG4 IgE
design <- createDesignMatrix(metadata(sce_IgE_clean_3D)$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_PBS_vs_IgE_3D <- createContrast(c(0, 0, 1, 0)) 
da_res_PBS_vs_IgE_3D <- diffcyt(sce_IgE_clean_3D, 
                                design = design, 
                                contrast = contrast_PBS_vs_IgE_3D,
                                analysis_type = "DA", 
                                method_DA = "diffcyt-DA-edgeR", 
                                clustering_to_use = "cluster_ID")

topTable(da_res_PBS_vs_IgE_3D, format_vals = TRUE)

ds_res_PBS_vs_IgE_3D <- diffcyt(sce_IgE_clean_3D, 
                                design = design, 
                                contrast = contrast_PBS_vs_IgE_3D,
                                analysis_type = "DS", 
                                method_DS = "diffcyt-DS-limma",
                                clustering_to_use = "cluster_ID")

topTable(ds_res_PBS_vs_IgE_3D, format_vals = TRUE)

# Differential analysis PBS HV vs M
design <- createDesignMatrix(metadata(filterSCE(sce_IgE_clean_3D, condition = "PBS"))$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_PBS_HV_vs_M_3D <- createContrast(c(0, 0, 0, 1)) 
da_res_PBS_HV_vs_M_3D <- diffcyt(filterSCE(sce_IgE_clean_3D, condition = "PBS"), 
                                 design = design, 
                                 contrast = contrast_PBS_HV_vs_M_3D,
                                 analysis_type = "DA", 
                                 method_DA = "diffcyt-DA-edgeR", 
                                 clustering_to_use = "cluster_ID")

topTable(da_res_PBS_HV_vs_M_3D, format_vals = TRUE)

ds_res_PBS_HV_vs_M_3D <- diffcyt(filterSCE(sce_IgE_clean_3D, condition = "PBS"), 
                                 design = design, 
                                 contrast = contrast_PBS_HV_vs_M_3D,
                                 analysis_type = "DS", 
                                 method_DS = "diffcyt-DS-limma",
                                 clustering_to_use = "cluster_ID")

topTable(ds_res_PBS_HV_vs_M_3D, format_vals = TRUE)

# Differential analysis NIP HV vs M
design <- createDesignMatrix(metadata(filterSCE(sce_IgE_clean_3D, condition = "NIP"))$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_NIP_HV_vs_M_3D <- createContrast(c(0, 0, 0, 1)) 
da_res_NIP_HV_vs_M_3D <- diffcyt(filterSCE(sce_IgE_clean_3D, condition = "NIP"), 
                                 design = design, 
                                 contrast = contrast_NIP_HV_vs_M_3D,
                                 analysis_type = "DA", 
                                 method_DA = "diffcyt-DA-edgeR", 
                                 clustering_to_use = "cluster_ID")

topTable(da_res_NIP_HV_vs_M_3D, format_vals = TRUE)

ds_res_NIP_HV_vs_M_3D <- diffcyt(filterSCE(sce_IgE_clean_3D, condition = "NIP"), 
                                 design = design, 
                                 contrast = contrast_NIP_HV_vs_M_3D,
                                 analysis_type = "DS", 
                                 method_DS = "diffcyt-DS-limma",
                                 clustering_to_use = "cluster_ID")

topTable(ds_res_NIP_HV_vs_M_3D, format_vals = TRUE)

# Differential analysis CSPG4 HV vs M
design <- createDesignMatrix(metadata(filterSCE(sce_IgE_clean_3D, condition = "CSPG4"))$experiment_info, 
                             cols_design = c("condition", "sample_type"))
contrast_CSPG4_HV_vs_M_3D <- createContrast(c(0, 0, 0, 1)) 
da_res_CSPG4_HV_vs_M_3D <- diffcyt(filterSCE(sce_IgE_clean_3D, condition = "CSPG4"), 
                                   design = design, 
                                   contrast = contrast_CSPG4_HV_vs_M_3D,
                                   analysis_type = "DA", 
                                   method_DA = "diffcyt-DA-edgeR", 
                                   clustering_to_use = "cluster_ID")

topTable(da_res_CSPG4_HV_vs_M_3D, format_vals = TRUE)

ds_res_CSPG4_HV_vs_M_3D <- diffcyt(filterSCE(sce_IgE_clean_3D, condition = "CSPG4"), 
                                   design = design, 
                                   contrast = contrast_CSPG4_HV_vs_M_3D,
                                   analysis_type = "DS", 
                                   method_DS = "diffcyt-DS-limma",
                                   clustering_to_use = "cluster_ID")

topTable(ds_res_CSPG4_HV_vs_M_3D, format_vals = TRUE)
#=================================================================================
# ADD P value
# cluster_abundance_condition_IgE

p <- plotAbundances(sce_IgE_clean_3D, k = "cluster_ID", by = "cluster_id", shape_by = "sample_type") + 
  facet_wrap(~cluster_id, scales = "free_y", ncol = 7) +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = alpha(condition_colors, 0.6)) + 
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))

p$data <- p$data %>% 
  dplyr::filter(cluster_id %in% c(1, 2, 3))
plot_data <- p$data

# Function to clean and format stats
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  if (!"cluster_id" %in% colnames(res_table)) {
    res_table <- res_table %>% rownames_to_column("cluster_id")
  }
  res_table %>%
    select(cluster_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}

# Process all 4 comparisons
s2 <- prep_stats(da_res_NIP_vs_IgE_3D,  "2") # NIP_IgE(3) vs CSPG4_IgE(5)
s4 <- prep_stats(da_res_PBS_vs_IgE_3D,  "4") # PBS(1) vs CSPG4_IgE(5)

label_data_final <- p$data %>%
  group_by(cluster_id) %>%
  summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE)) %>% 
  left_join(s2, by = "cluster_id") %>%
  left_join(s4, by = "cluster_id") %>%
  mutate(cluster_id = factor(cluster_id, levels = levels(p$data$cluster_id))) %>%
  mutate(
    y4 = local_max * 1.15, y4_t = local_max * 1.17, # Tier 2
    y2 = local_max * 1.05, y2_t = local_max * 1.07  # Tier 4
  )

p_final <- p + 
  facet_wrap(~cluster_id, scales = "free_y", ncol = 5) + 
  geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
               aes(x = 1, xend = 3, y = y4, yend = y4), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
            aes(x = 2, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 4, vjust = 0) +
  geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
               aes(x = 2, xend = 3, y = y2, yend = y2), inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
            aes(x = 2.5, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 4, vjust = 0) +
  theme(panel.spacing.x = unit(1, "cm"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        text = element_text(size = 24),
        plot.title = element_text(face = "bold", size = 24))

p_final
png("cluster_abundance_resting_IgE_3D_padj.png", width = 10, height = 4, units = "in", res = 300)
p_final
dev.off()

#=============================================================================================
#Differential abundance analysis with sub boxplot with P value
# 1. Generate the base object to get the data structure
p <- plotAbundances(sce_IgE_clean_3D, k = "cluster_ID", by = "cluster_id", shape_by = "sample_type")

plot_data <- p$data %>% 
  dplyr::filter(cluster_id %in% c(3, 7, 9))

# Function to clean and format stats
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  if (!"cluster_id" %in% colnames(res_table)) {
    res_table <- res_table %>% rownames_to_column("cluster_id")
  }
  res_table %>%
    select(cluster_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}

s1 <- prep_stats(da_res_PBS_HV_vs_M_3D,  "1")
s2 <- prep_stats(da_res_NIP_HV_vs_M_3D,  "2")
s3 <- prep_stats(da_res_CSPG4_HV_vs_M_3D,  "3")

label_data_final <- plot_data %>%
  group_by(cluster_id) %>%
  summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE)) %>% 
  left_join(s1, by = "cluster_id") %>%
  left_join(s2, by = "cluster_id") %>%
  left_join(s3, by = "cluster_id") %>%
  mutate(cluster_id = factor(cluster_id, levels = levels(plot_data$cluster_id))) %>%
  mutate(
    y1 = local_max * 1.05, y1_t = local_max * 1.07,
    y2 = local_max * 1.10, y2_t = local_max * 1.12,
    y3 = local_max * 1.15, y3_t = local_max * 1.17
  )

p_final <- ggplot(plot_data, aes(x = condition, y = Freq, color = condition)) +
  
  geom_boxplot(aes(fill = condition, alpha = sample_type),
               width = 0.6,
               position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               show.legend = FALSE) + 
  
  geom_point(aes(shape = sample_type, group = sample_type), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5, 
             stroke = 1, 
             fill = "transparent",
             show.legend = TRUE) + 
  
  facet_wrap(~cluster_id, scales = "free_y", ncol = 4) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  scale_alpha_manual(values = c("HV" = 0.1, "M" = 0.3)) + 
  scale_shape_manual(values = c("HV" = 21, "M" = 24)) +
  
  geom_segment(data = dplyr::filter(label_data_final, p_label_1 != ""),
               aes(x =0.8, xend = 1.2, y = y1, yend = y1), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_1 != ""),
            aes(x = 1, y = y1_t, label = p_label_1), 
            inherit.aes = FALSE, size = 5, color = "black") +
  
  geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
               aes(x = 1.8, xend = 2.2, y = y2, yend = y2), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
            aes(x = 2, y = y2_t, label = p_label_2), 
            inherit.aes = FALSE, size = 5, color = "black") +
  
  geom_segment(data = dplyr::filter(label_data_final, p_label_3 != ""),
               aes(x = 2.8, xend = 3.2, y = y3, yend = y3), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_3 != ""),
            aes(x = 3, y = y3_t, label = p_label_3), 
            inherit.aes = FALSE, size = 5, color = "black") +
  
  theme_bw() +
  labs(y = "Proportion [%]", x = NULL) +
  theme(
    text = element_text(size = 20),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold", size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "cm")
  )

p_final

#=================================================================================
# loop for marker expression in IgE in multiple clusters
target_clusters <- as.character(1:14)
dir.create("pbExprs_cluster_IgE_clean_3D", showWarnings = FALSE)
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  res_table %>%
    select(cluster_id, marker_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}
s2 <- prep_stats(ds_res_NIP_vs_IgE_3D, "2")
s4 <- prep_stats(ds_res_PBS_vs_IgE_3D, "4")

lapply(target_clusters, function(id) {
  message(paste("Generating plot for cluster:", id))
  sce_sub <- filterSCE(sce_IgE_clean_3D, cluster_ID == id)
  p <- plotPbExprs(sce_sub, k = "cluster_ID", 
                   features = markers_ex_CD14, 
                   facet_by = "cluster_id", 
                   shape_by = "sample_type", ncol = 4) + 
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    scale_shape_manual(values = c(1, 2)) +
    theme(text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0,    
                                     hjust = 0.5,  
                                     vjust = 0.5))
  p$data <- p$data %>% 
    dplyr::filter(cluster_id == id)
  vals <- as.numeric(as.character(p$data$value))
  y_max <- max(vals, na.rm = TRUE)
  y_min <- min(vals, na.rm = TRUE)
  y_range <- y_max - y_min
  jump <- y_range * 0.12
  
  label_data_final <- p$data %>%
    mutate(value = as.numeric(as.character(value))) %>%
    group_by(antigen, cluster_id) %>%
    summarize(local_max = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(s2, by = c("antigen" = "marker_id", "cluster_id")) %>%
    left_join(s4, by = c("antigen" = "marker_id", "cluster_id")) %>%
    mutate(
      y2 = local_max + (jump * 0.5),   # Bar 1
      y2_t = local_max + (jump * 0.7), # Text 1
      y4 = local_max + (jump * 1.5),   # Bar 2
      y4_t = local_max + (jump * 1.7)  # Text 2
    )
  
  p$layers[[1]]$position <- position_dodge(width = 0.9)
  p$layers[[1]]$geom_params$width <- 0.6 
  p$layers[[2]]$position <- position_dodge(width = 0.9)
  
  marker_levels <- levels(factor(p$data$antigen))
  label_data_final <- label_data_final %>%
    mutate(
      x_pos = as.numeric(factor(antigen, levels = marker_levels)),
      x_s2 = x_pos + offsets[2], 
      x_e2 = x_pos + offsets[3],
      x_s4 = x_pos + offsets[1], 
      x_e4 = x_pos + offsets[3])
  
  p_final <- p + 
    geom_segment(data = dplyr::filter(label_data_final, p_label_4 != ""),
                 aes(x = x_s4, xend = x_e4, y = y4, yend = y4), inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(label_data_final, p_label_4 != ""),
              aes(x = (x_s4 + x_e4)/2, y = y4_t, label = p_label_4), inherit.aes = FALSE, size = 3.5, vjust = 0) +
    geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
                 aes(x = x_s2, xend = x_e2, y = y2, yend = y2), inherit.aes = FALSE) +
    geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
              aes(x = (x_s2 + x_e2)/2, y = y2_t, label = p_label_2), inherit.aes = FALSE, size = 3.5, vjust = 0) +
    coord_cartesian(ylim = c(y_min, y_max + (jump * 2.5))) 
  
  file_name <- paste0("pbExprs_cluster_IgE_clean_3D/pbExprs_cluster_", id, "_sig_padj.png")
  png(file_name, width = 15, height = 5, units = "in", res = 300)
  print(p_final)
  dev.off()
})
#=================================================================================
# Density plot
library(ggh4x)

plotClusterExprs(sce_IgE_clean_3D, k = "cluster_ID", features = markers_ex_CD14)

plot_df <- as.data.frame(t(assay(sce_IgE_clean_3D, "exprs")[markers_ex_CD14, ]))
cl_ids <- cluster_ids(sce_IgE_clean_3D, "cluster_ID")
plot_df$cluster_id <- factor(as.character(cl_ids), levels = rev(as.character(1:14)))
plot_df_long <- plot_df %>%
  tidyr::pivot_longer(
    cols = -cluster_id, 
    names_to = "antigen", 
    values_to = "expression")
plot_df_long$antigen <- factor(plot_df_long$antigen, levels = as.character(markers_ex_CD14))

groups <- list(
  "FcRs" = c("CD64", "CD32B", "CD16", "FCER1", "CD23"),
  "Activation/Antigen-presenting" = c("CD80", "CD86", "CD40", "CCR2", "HLADR"),
  "Inhibition/Exhaustion" = c("CD163", "CD206", "PDL1", " ", "  ")
)
medians_cl1 <- plot_df_long %>%
  filter(cluster_id == "1") %>%
  group_by(antigen) %>%
  summarize(m_val = median(expression, na.rm = TRUE))


facet_themes <- list(
  " " = theme(panel.grid.major = element_blank(), 
              panel.border = element_blank(),
              axis.line = element_blank()),
  "  " = theme(panel.grid.major = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_blank())
)

make_ridge_row <- function(data, markers, title_text, show_x = FALSE) {
  df_sub <- data %>% 
    filter(antigen %in% markers, cluster_id %in% c("11", "12", "13", "14")) %>%
    mutate(antigen = factor(antigen, levels = markers))  
  meds_sub <- medians_cl1 %>% 
    filter(antigen %in% markers) %>%
    mutate(antigen = factor(as.character(antigen), levels = markers))
  p <- ggplot(df_sub, aes(x = expression, y = cluster_id, fill = cluster_id, color = cluster_id)) +
    # Add vertical lines first so they sit BEHIND the ridges
    geom_vline(data = meds_sub, aes(xintercept = m_val), 
               color = "grey20", linetype = "dashed", linewidth = 0.5) +
    # The ridges
    geom_density_ridges(alpha = 0.7, scale = 1.5, linewidth = 0.5) +
    # Apply your cluster_colour vector
    scale_fill_manual(values = cluster_colors) + 
    scale_color_manual(values = cluster_colors) +
    # Separate by antigen group
    facet_wrap(~antigen, nrow = 1, ncol = 5, scales = "free_x", drop = FALSE) +
    facetted_pos_scales(x = list(antigen %in% c(" ", "  ") ~ scale_x_continuous(guide = "none"))) +
    labs(subtitle = title_text, x = NULL, y = "Cluster ID") +
    # Styling to match the original
    theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
    theme(
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 14),
      strip.background = element_blank(), # Removes grey box behind antigen names
      panel.grid.major = element_line(color = "grey90"), # Keeps the background grid
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 14, margin = margin(b = 5)), # Makes antigen names bold
      panel.background = element_blank(),      # Ensures clean background
      legend.position = "none",
      panel.spacing = unit(1, "lines")
    )
  p <- p + theme(panel.grid.major.y = element_blank())
  if (!show_x) { 
    p <- p + theme(axis.text.x = element_blank(), 
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank()) 
  }
  return(p)
}
p1 <- make_ridge_row(plot_df_long, groups[[1]], names(groups)[1])
p2 <- make_ridge_row(plot_df_long, groups[[2]], names(groups)[2])
p3 <- make_ridge_row(plot_df_long, groups[[3]], names(groups)[3], show_x = TRUE)

final_plot <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))
final_plot
png("densities_IgE_3D_exhausted_padj.png", width = 9, height = 6, units = "in", res = 300)
final_plot
dev.off()
#==================================================================================================
#==================================================================================================
#==================================================================================================
# Analysis with different conditions
# Differential analysis 2D vs 3D in PBS
sce_IgE_clean_PBS <- filterSCE(sce_IgE_clean, condition == "PBS")

design <- createDesignMatrix(metadata(sce_IgE_clean_PBS)$experiment_info, 
                             cols_design = c("setting", "sample_type"))
contrast_2D_vs_3D_PBS  <- createContrast(c(0, 1, 0))
da_res_2D_vs_3D_PBS <- diffcyt(sce_IgE_clean_PBS, 
                                design = design, 
                                contrast = contrast_2D_vs_3D_PBS,
                                analysis_type = "DA", 
                                method_DA = "diffcyt-DA-edgeR", 
                                clustering_to_use = "cluster_ID")

topTable(da_res_2D_vs_3D_PBS, format_vals = TRUE)

ds_res_2D_vs_3D_PBS <- diffcyt(sce_IgE_clean_PBS, 
                                design = design, 
                                contrast = contrast_2D_vs_3D_PBS,
                                analysis_type = "DS", 
                                method_DS = "diffcyt-DS-limma",
                                clustering_to_use = "cluster_ID")

topTable(ds_res_2D_vs_3D_PBS, format_vals = TRUE)

# Differential analysis 2D vs 3D in NIP IgE
sce_IgE_clean_NIP <- filterSCE(sce_IgE_clean, condition == "NIP IgE")

design <- createDesignMatrix(metadata(sce_IgE_clean_NIP)$experiment_info, 
                             cols_design = c("setting", "sample_type"))
contrast_2D_vs_3D_NIP  <- createContrast(c(0, 1, 0))
da_res_2D_vs_3D_NIP <- diffcyt(sce_IgE_clean_NIP, 
                               design = design, 
                               contrast = contrast_2D_vs_3D_NIP,
                               analysis_type = "DA", 
                               method_DA = "diffcyt-DA-edgeR", 
                               clustering_to_use = "cluster_ID")

topTable(da_res_2D_vs_3D_NIP, format_vals = TRUE)

ds_res_2D_vs_3D_NIP <- diffcyt(sce_IgE_clean_NIP, 
                               design = design, 
                               contrast = contrast_2D_vs_3D_NIP,
                               analysis_type = "DS", 
                               method_DS = "diffcyt-DS-limma",
                               clustering_to_use = "cluster_ID")

topTable(ds_res_2D_vs_3D_NIP, format_vals = TRUE)

# Differential analysis 2D vs 3D in CSPG4 IgE
sce_IgE_clean_CSPG4 <- filterSCE(sce_IgE_clean, condition == "CSPG4 IgE")

design <- createDesignMatrix(metadata(sce_IgE_clean_CSPG4)$experiment_info, 
                             cols_design = c("setting", "sample_type"))
contrast_2D_vs_3D_CSPG4  <- createContrast(c(0, 1, 0))
da_res_2D_vs_3D_CSPG4 <- diffcyt(sce_IgE_clean_CSPG4, 
                               design = design, 
                               contrast = contrast_2D_vs_3D_CSPG4,
                               analysis_type = "DA", 
                               method_DA = "diffcyt-DA-edgeR", 
                               clustering_to_use = "cluster_ID")

topTable(da_res_2D_vs_3D_CSPG4, format_vals = TRUE)

ds_res_2D_vs_3D_CSPG4 <- diffcyt(sce_IgE_clean_CSPG4, 
                               design = design, 
                               contrast = contrast_2D_vs_3D_CSPG4,
                               analysis_type = "DS", 
                               method_DS = "diffcyt-DS-limma",
                               clustering_to_use = "cluster_ID")

topTable(ds_res_2D_vs_3D_CSPG4, format_vals = TRUE)
#=================================================================================
# ADD P value
# cluster_abundance_condition_IgE
p <- plotAbundances(sce_IgE_clean_PBS, k = "cluster_ID", by = "cluster_id", shape_by = "setting")

plot_data <- p$data #%>% 
  #dplyr::filter(cluster_id %in% c(1, 2, 3))

# Function to clean and format stats
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  if (!"cluster_id" %in% colnames(res_table)) {
    res_table <- res_table %>% rownames_to_column("cluster_id")
  }
  res_table %>%
    select(cluster_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}

s1 <- prep_stats(da_res_2D_vs_3D_PBS,  "1")
#s2 <- prep_stats(da_res_2D_vs_3D_NIP,  "2")
#s3 <- prep_stats(da_res_2D_vs_3D_CSPG4,  "3")

label_data_final <- plot_data %>%
  group_by(cluster_id) %>%
  summarize(local_max = max(as.numeric(as.character(Freq)), na.rm = TRUE)) %>% 
  left_join(s1, by = "cluster_id") %>%
  #left_join(s2, by = "cluster_id") %>%
  #left_join(s3, by = "cluster_id") %>%
  mutate(cluster_id = factor(cluster_id, levels = levels(plot_data$cluster_id))) %>%
  mutate(
    y1 = local_max * 1.03, y1_t = local_max * 1.07,
   # y2 = local_max * 1.07, y2_t = local_max * 1.12,
   # y3 = local_max * 1.13, y3_t = local_max * 1.17
  )

label_data_final <- label_data_final %>% dplyr::filter(p_label_1 != "")
sig_cluster_id <- label_data_final$cluster_id
plot_data <- plot_data %>% dplyr::filter(cluster_id %in% sig_cluster_id)

p_final <- ggplot(plot_data, aes(x = condition, y = Freq, color = condition)) +
  
  geom_boxplot(aes(fill = condition, alpha = setting),
               width = 0.6,
               position = position_dodge(width = 0.8), 
               outlier.shape = NA, 
               show.legend = FALSE) + 
  
  geom_point(aes(shape = setting, group = setting), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 1.5, 
             stroke = 1, 
             fill = "transparent",
             show.legend = TRUE) + 
  
  facet_wrap(~cluster_id, scales = "free_y", ncol = 4) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  scale_alpha_manual(values = c("2D" = 0.1, "3D" = 0.5)) + 
  scale_shape_manual(values = c("2D" = 21, "3D" = 24)) +
  
  geom_segment(data = dplyr::filter(label_data_final, p_label_1 != ""),
               aes(x =0.8, xend = 1.2, y = y1, yend = y1), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = dplyr::filter(label_data_final, p_label_1 != ""),
            aes(x = 1, y = y1_t, label = p_label_1), 
            inherit.aes = FALSE, size = 4, color = "black") +
  
 # geom_segment(data = dplyr::filter(label_data_final, p_label_2 != ""),
              # aes(x = 1.8, xend = 2.2, y = y2, yend = y2), 
              # inherit.aes = FALSE, color = "black") +
 # geom_text(data = dplyr::filter(label_data_final, p_label_2 != ""),
           # aes(x = 2, y = y2_t, label = p_label_2), 
           #inherit.aes = FALSE, size = 4, color = "black") +
  
 # geom_segment(data = dplyr::filter(label_data_final, p_label_3 != ""),
               #aes(x = 2.8, xend = 3.2, y = y3, yend = y3), 
              # inherit.aes = FALSE, color = "black") +
 # geom_text(data = dplyr::filter(label_data_final, p_label_3 != ""),
           # aes(x = 3, y = y3_t, label = p_label_3), 
           # inherit.aes = FALSE, size = 4, color = "black") +
  
  theme_bw() +
  labs(y = "Proportion [%]", x = NULL) +
  theme(
    text = element_text(size = 20),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold", size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "cm")
  )

p_final
png("cluster_abundance_resting_IgE_PBS_2Dvs3D_padj.png", width = 13, height = 4, units = "in", res = 300)
p_final
dev.off()
#=====================================================================================================
# loop for marker expression in IgE in multiple clusters
target_clusters <- as.character(1:14)
dir.create("pbExprs_cluster_IgE_clean_PBS", showWarnings = FALSE)
prep_stats <- function(da_obj, name) {
  res_table <- as.data.frame(rowData(da_obj$res))
  res_table %>%
    select(cluster_id, marker_id, p_adj) %>%
    mutate(!!paste0("p_label_", name) := ifelse(p_adj < 0.0001, "<0.0001", 
                                                ifelse(p_adj < 0.05, sprintf("%.4f", p_adj), ""))) %>%
    select(-p_adj)
}
s1 <- prep_stats(ds_res_2D_vs_3D_PBS, "1")

lapply(target_clusters, function(id) {
  message(paste("Generating plot for cluster:", id))
  sce_sub <- filterSCE(sce_IgE_clean_CSPG4, cluster_ID == id)
  p <- plotPbExprs(sce_sub, k = "cluster_ID", 
                   features = markers_ex_CD14, 
                   facet_by = "cluster_id", 
                   shape_by = "setting", ncol = 4) + 
    scale_color_manual(values = condition_colors) +
    scale_fill_manual(values = condition_colors) +
    scale_shape_manual(values = c(1, 2)) +
    theme(text = element_text(size = 20),
          plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0,    
                                     hjust = 0.5,  
                                     vjust = 0.5))
  plot_data <- p$data %>% 
    dplyr::filter(cluster_id == id)
  vals <- as.numeric(as.character(p$data$value))
  y_max <- max(vals, na.rm = TRUE)
  y_min <- min(vals, na.rm = TRUE)
  y_range <- y_max - y_min
  jump <- y_range * 0.12
  
  label_data_final <- p$data %>%
    mutate(value = as.numeric(as.character(value))) %>%
    group_by(antigen, cluster_id) %>%
    summarize(local_max = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(s1, by = c("antigen" = "marker_id", "cluster_id")) %>%
    mutate(
      y1 = local_max + (jump * 0.5),   # Bar 1
      y1_t = local_max + (jump * 0.57), # Text 1
    ) %>% dplyr::filter(cluster_id == id)
  
  label_data_final <- label_data_final %>% dplyr::filter(p_label_1 != "")
  sig_antigen <- label_data_final$antigen
  plot_data <- plot_data %>% dplyr::filter(antigen %in% sig_antigen)  
  
  p$layers[[1]]$position <- position_dodge(width = 0.9)
  p$layers[[1]]$geom_params$width <- 0.6 
  p$layers[[2]]$position <- position_dodge(width = 0.9)
  
  marker_levels <- levels(factor(plot_data$antigen))
  label_data_final <- label_data_final %>%
    mutate(
      x_pos = as.numeric(factor(antigen, levels = marker_levels)),
      x_s1 = x_pos + offsets[1], 
      x_e1 = x_pos + offsets[2])
  
  p_final <- ggplot(plot_data, aes(x = antigen, y = value, color = condition)) +
    
    geom_boxplot(aes(fill = condition, alpha = setting),
                 width = 0.7,
                 position = position_dodge(width = 0.8), 
                 outlier.shape = NA, 
                 show.legend = FALSE) + 
    
    geom_point(aes(shape = setting, group = setting), 
               position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
               size = 1.5, 
               stroke = 0.8, 
               fill = "transparent",
               show.legend = TRUE) + 
    
    # Facet by antigen to keep markers in separate panels (similar to your wrap logic)
    scale_color_manual(values = condition_colors, limits = force) +
    scale_fill_manual(values = condition_colors, limits = force) +
    scale_alpha_manual(values = c("2D" = 0.2, "3D" = 0.6)) + 
    scale_shape_manual(values = c("2D" = 21, "3D" = 24)) +
    
    # Significance Brackets (Dodged over the center of the 2D and 3D boxes)
    geom_segment(data = dplyr::filter(label_data_final, p_label_1 != ""),
                 aes(x = x_s1, xend = x_e1, y = y1, yend = y1), 
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = dplyr::filter(label_data_final, p_label_1 != ""),
              aes(x = x_pos, y = y1_t, label = p_label_1), 
              inherit.aes = FALSE, size = 4, color = "black") +
    
    theme_bw() +
    labs(title = paste(id), y = "median expression", x = NULL) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      text = element_text(size = 18),
      strip.background = element_blank(), 
      axis.text.x = element_text(size = 16)
    )
  
  file_name <- paste0("pbExprs_cluster_IgE_clean_CSPG4/pbExprs_cluster_", id, "_sig_padj.png")
  png(file_name, width = 15, height = 5, units = "in", res = 300)
  print(p_final)
  dev.off()
})
#=================================================================================

sce_IgE_clean_PBS_test <-sce_IgE_clean_PBS
sce_IgE_clean_PBS_test <- runDR(sce_IgE_clean_PBS_test, cell = 1000, dr = "TSNE", features = "state", perplexity= 25, seed = 10)
png("IgE_clean_tSNE_PBS.png", width = 8, height = 7, units = "in", res = 300)
plotDR(sce_IgE_clean_PBS_test,dr = "TSNE", color_by = "cluster_ID")  +
  scale_color_manual(values = cluster_colors) +
  theme(text = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 24),
        axis.title = element_text(size = 20))
dev.off()









