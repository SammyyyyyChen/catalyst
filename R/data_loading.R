# ---- FCS file helpers -----------------------------------------------------

list_fcs_files <- function(fcs_dir) {
  list.files(fcs_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
}

validate_fcs_files <- function(file_names) {
  bad <- c()
  for (f in file_names) {
    ok <- tryCatch({
      read.FCS(f, transformation = FALSE)
      TRUE
    }, error = function(e) {
      message("!!! ERROR in file: ", basename(f), "\n    ", e$message)
      FALSE
    })
    if (!ok) bad <- c(bad, f)
  }
  if (length(bad) == 0) {
    message("All ", length(file_names), " files read successfully.")
  } else {
    message("Found ", length(bad), " corrupted files:")
    print(basename(bad))
  }
  invisible(bad)
}

# ---- Metadata from filenames ---------------------------------------------

build_metadata <- function(file_names) {
  fb <- basename(file_names)
  data.frame(file_name = fb) %>%
    mutate(
      setting     = str_extract(file_name, "2D|3D"),
      sample_type = str_extract(file_name, "HV|M"),
      patient_id  = str_extract(file_name, "(HV|M)[0-9]+"),
      condition   = str_extract(file_name, "PBS|NIP IgG1|CSPG4 IgG1|NIP IgE|CSPG4 IgE"),
      sample_id   = paste(patient_id, condition, setting, sep = "_"),
      batch = case_when(
        patient_id %in% c("HV778", "HV807", "HV808") ~ "20251117",
        patient_id %in% c("HV810", "HV814", "HV815") ~ "20251205",
        patient_id %in% c("M764", "M770", "M781")    ~ "20251212",
        patient_id %in% c("M817", "M818", "HV820")   ~ "20251219",
        patient_id %in% c("M748", "M832", "HV830")   ~ "20250131",
        patient_id %in% c("M834", "HV835", "HV836")  ~ "20250206"
      )
    )
}

# ---- Panel from first FCS ------------------------------------------------

build_panel <- function(file_names) {
  ff <- read.FCS(file_names[1], transformation = FALSE, truncate_max_range = FALSE)
  params <- pData(parameters(ff))

  panel <- tibble(
    fcs_colname = params$name,
    antigen     = params$desc
  ) %>%
    mutate(antigen = gsub("[[:space:]]+[A-Za-z0-9-]+-A$", "", antigen))

  panel$antigen[is.na(panel$antigen)] <- panel$fcs_colname[is.na(panel$antigen)]

  panel <- panel %>%
    mutate(
      marker_class = case_when(
        antigen %in% c("FSC-A", "SSC-A", "Time", "DAPI-A", "FSC-Width") ~ "none",
        TRUE ~ "state"
      )
    )
  panel
}

# ---- Build SingleCellExperiment ------------------------------------------

load_and_prep_sce <- function(file_names, panel, metadata) {
  sce <- prepData(
    x        = file_names,
    panel    = panel,
    md       = metadata,
    features = panel$fcs_colname[panel$marker_class == "state"],
    cofactor = 400,
    md_cols  = list(
      file    = "file_name",
      id      = "sample_id",
      factors = c("sample_type", "patient_id", "condition", "setting", "batch")
    ),
    ignore.text.offset = TRUE,
    emptyValue = FALSE
  )
  sce
}

# ---- Filter to IgE conditions and set factor baselines -------------------

prep_ige_sce <- function(sce) {
  sce_ige <- filterSCE(sce, condition %in% c("PBS", "NIP IgE", "CSPG4 IgE"))

  ei <- metadata(sce_ige)$experiment_info
  ei$condition   <- factor(ei$condition,   levels = c("PBS", "NIP IgE", "CSPG4 IgE"))
  ei$setting     <- factor(ei$setting,     levels = c("2D", "3D"))
  ei$sample_type <- factor(ei$sample_type, levels = c("HV", "M"))
  metadata(sce_ige)$experiment_info <- ei

  sce_ige$condition <- factor(sce_ige$condition, levels = c("PBS", "NIP IgE", "CSPG4 IgE"))

  sce_ige
}
