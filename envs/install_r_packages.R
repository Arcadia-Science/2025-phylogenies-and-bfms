# Install R packages not available (or not version-pinned) via conda.
# Run this once after creating the conda environment from envs/dev.yml:
#   conda env create -f envs/dev.yml
#   conda activate phylogenies-bfms
#   Rscript envs/install_r_packages.R

# remotes provides install_version() for installing specific CRAN versions.
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_version("here", version = "1.0.1", repos = "https://cloud.r-project.org")
remotes::install_version("BiocManager", version = "1.30.25", repos = "https://cloud.r-project.org")

# arcadiathemeR is only available from GitHub; pin to a tagged release with `ref`.
remotes::install_github("Arcadia-Science/arcadiathemeR", ref = "v0.1.0")
