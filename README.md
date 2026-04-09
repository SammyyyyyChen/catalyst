# CATALYST IgE Flow Cytometry Analysis

Analysis of monocyte phenotyping by high-dimensional flow cytometry in a
melanoma spheroid co-culture model treated with anti-CSPG4 IgE antibody.

## Experiment overview

- **Cell line**: A2058 melanoma spheroids
- **Culture settings**: 2D monolayer vs 3D spheroid
- **Treatments**: PBS, NIP IgE (isotype control), CSPG4 IgE (therapeutic)
- **Monocyte sources**: Healthy volunteers (HV) and melanoma patients (M)

## Setup

1. Clone this repo and open `catalyst.Rproj` in RStudio.
2. Copy `config_template.R` to `config.R` and fill in the path to your data
   directory (the folder containing `raw files test/` with `.fcs` files).
3. Install required R packages (see `R/packages.R`).
4. Run `analysis/01_IgE_analysis.R` section by section in RStudio.

## Project structure

```
R/                 Helper functions (sourced by analysis scripts)
  packages.R       Package loading
  theme.R          Colour palettes, marker lists, shared themes
  data_loading.R   FCS import, metadata & panel construction
  clustering.R     FlowSOM clustering and merge-table definitions
  differential.R   diffcyt wrappers and statistics helpers
  plotting.R       Reusable figure-generation functions
analysis/          Main analysis scripts
  01_IgE_analysis.R
output/            Generated figures and saved objects (gitignored)
  figures/
  rds/
config_template.R  Template for local path configuration
```
