# ---- Colour palettes -----------------------------------------------------

condition_colors <- c(
  "PBS"       = "#BDBDBD",
  "NIP IgE"   = "#e3c7b6",
  "CSPG4 IgE" = "#9C2830",
  "2D" = "#89d6a4",
  "3D" = "#f5b871",
  "HV" = "lightskyblue2",
  "M"  = "coral1"
)

cluster_colors <- c(
  "1"  = "#FB8072",
  "2"  = "#E78AC3",
  "3"  = "#FDB462",
  "4"  = "#DC050C",
  "5"  = "#33A02C",
  "6"  = "#E7298A",
  "7"  = "#B2DF8A",
  "8"  = "#1D74CC",
  "9"  = "#882E72",
  "10" = "#FF7F00",
  "11" = "#8DD3C7",
  "12" = "#7BAFDE",
  "13" = "#E6AB02",
  "14" = "#B17BA6"
)

# ---- Marker groups -------------------------------------------------------

ordered_markers <- c(
  "CD14", "CD16", "CD64", "CD32B", "FCER1", "CD23",
  "CD80", "CD86", "CD40", "CCR2", "HLADR",
  "CD163", "CD206", "PDL1"
)

markers_ex_CD14 <- c(
  "CD64", "CD32B", "CD16", "FCER1", "CD23",
  "CD80", "CD86", "CD40", "CCR2", "HLADR",
  "CD163", "CD206", "PDL1"
)

markers_FcRs      <- c("FCER1", "CD23")
markers_activation <- c("CD80", "CD86", "CD40", "CCR2", "HLADR")
markers_inhibitory <- c("CD163", "CD206", "PDL1")

ridge_groups <- list(
  "FcRs"                       = c("CD64", "CD32B", "CD16", "FCER1", "CD23"),
  "Activation/Antigen-presenting" = c("CD80", "CD86", "CD40", "CCR2", "HLADR"),
  "Inhibition/Exhaustion"      = c("CD163", "CD206", "PDL1", " ", "  ")
)

# ---- Shared ggplot theme helpers -----------------------------------------

theme_pub <- function(base_size = 20) {

  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title  = element_text(face = "bold", size = base_size + 4, hjust = 0.5),
      axis.title  = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size),
      strip.text   = element_text(size = base_size)
    )
}

# Dodge parameters used for boxplot significance brackets
n_cond  <- 3
d_width <- 0.9
step    <- d_width / n_cond
offsets <- (seq_len(n_cond) - (n_cond + 1) / 2) * step
