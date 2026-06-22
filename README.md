# Phylogenies and Biological Foundation Models

## Purpose
This repo contains code associated with the pub "Phylogenies and Biological Foundation Models".

## Installation and Setup
This repository uses R and conda to manage software environments and installations. All analyses were developed on an Apple MacBook Pro running macOS Monterey (see [Compute Specifications](#compute-specifications)). You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).

1. **Create and activate the conda environment.** This installs R, Python, and most dependencies.
   ```
   conda env create -f envs/dev.yml
   conda activate phylogenies-bfms
   ```

2. **Install the remaining R packages.** A few R packages (`here`, `BiocManager`, and the GitHub-only `arcadiathemeR`) are not available through conda and must be installed with the helper script. Run this once after activating the environment:
   ```
   Rscript envs/install_r_packages.R
   ```

3. **Set the NVIDIA API key (Figure 3 only).** Figure 3 queries the hosted Evo 2 model through NVIDIA's API. Obtain a key from [build.nvidia.com/arc/evo2-40b](https://build.nvidia.com/arc/evo2-40b) and make it available to R, for example by adding the following line to `~/.Renviron`:
   ```
   NVIDIA_API_KEY=your-key-here
   ```

## Running the analyses
Run all scripts from the repository root so that `here()` and the data-download commands resolve paths correctly. The figure scripts source `code/utils/utils.R` automatically and download or extract the data they need.

```
Rscript code/analysis/fig1_analysis.R
Rscript code/analysis/fig2_analysis.R
Rscript code/analysis/fig3_analysis.R   # requires NVIDIA_API_KEY (see step 3 above)
Rscript code/analysis/fig4_analysis.R   # requires Zenodo data (see Data section below)
```

The scripts are independent of one another and can be run in any order, but Figure 4 requires the Pfam data to be downloaded first (see below).

## Data

### Overview
The repository contains code to collect most of the necessary data to reproduce analyses in the pub. Figures 2 and 3 download their inputs automatically when the corresponding script is run. Figure 4 requires the example Pfam protein families, which are distributed via [Zenodo](https://zenodo.org/records/15644457).

### Figure 4 data (Zenodo)
Download the example Pfam families archive from the [Zenodo deposit](https://zenodo.org/records/15644457) and extract its contents so that the FASTA alignments live in `data/example_pfam_families/`. From the repository root:
```
mkdir -p data/example_pfam_families
# Download the archive from https://zenodo.org/records/15644457, then extract it, e.g.:
unzip example_pfam_families.zip -d data/example_pfam_families/
```
After extraction, `data/example_pfam_families/` should contain the `.fasta` multiple sequence alignment files that `code/analysis/fig4_analysis.R` reads. The clustering step writes its results to `data/example_pfam_families/clustering_results/`.

### Description of the folder structure
The repository is organized into the following top-level directories:

**code**: R scripts and Jupyter notebooks used for data processing, phylogenetic analysis, and figure generation.
**data**: Phylogenetic trees (.tre files), protein family sequences (FASTA and Stockholm formats), and clustering results.

```
├── code
│   ├── analysis
│   │   ├── fig1_analysis.R
│   │   ├── fig2_analysis.R
│   │   ├── fig3_analysis.R
│   │   └── fig4_analysis.R
│   └── utils
│       ├── fig2_extract_trees.py
│       ├── fig3_download_human_genes.py
│       ├── fig4_cluster_fasta_directory.sh
│       └── utils.R
├── data
│   ├── cox1.tre
```

## Methods

1. **fig1_analysis.R**: Creates side-by-side fan plots of cox1 gene trees for plant (Streptophyta) and animal (Chordata) phylogenies with outlier removal.
2. **fig2_analysis.R**: Downloads Ensembl Compara trees, calculates Hill's diversity across >14,000 protein families, and visualizes diversity patterns with branch length variance analysis.
3. **fig3_analysis.R**: Integrates gene age data with evolutionary likelihood estimates using rolling window analysis to examine temporal trends across evolutionary time.
4. **fig4_analysis.R**: Compares Hill's diversity of multiple sequence alignments with MMseqs clustering patterns across 218 Pfam families.

## Compute Specifications
All analyses were done on an Apple MacBook Pro running macOS Montery with 32GB RAM, 10 cores, and 1TB of storage.

## Contributing
This is a research repository. For questions or collaboration inquiries, please contact the repository maintainers.
