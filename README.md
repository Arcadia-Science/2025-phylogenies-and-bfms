# Phylogenies and Biological Foundation Models

## Purpose
This repo contains code associated with the pub "Phylogenies and Biological Foundation Models".

## Installation and Setup
This repository uses R and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html). 

### R Dependencies
The project requires the following R packages:
- ape (phylogenetic analysis)
- ggplot2 (visualization)
- gridExtra (plot layouts)
- Biostrings (sequence analysis)
- arcadiathemeR (custom themes)
- pracma (mathematical functions)
- httr (HTTP requests)
- jsonlite (JSON parsing)
- parallel (parallel processing)
- pbapply (progress bars)
- dplyr (data manipulation)
- tidyr (data tidying)
- scales (plot scaling)

## Data

### Overview
The repository contains phylogenetic trees, protein family alignments, and clustering results for analyzing relationships between phylogenetic diversity and protein evolution.

### Description of the folder structure
The repository is organized into the following top-level directories:

**code**: R scripts and Jupyter notebooks used for data processing, phylogenetic analysis, and figure generation.
**data**: Phylogenetic trees (.tre files), protein family sequences (FASTA and Stockholm formats), and clustering results.

```
├── code
│   ├── notebooks
│   │   ├── fig1_analysis.ipynb
│   │   ├── fig2_analysis.ipynb
│   │   ├── fig3_analysis.ipynb
│   │   ├── fig4_analysis.ipynb
│   │   └── fig5_analysis.ipynb
│   └── utils
│       ├── fig2_extract_trees.py
│       ├── fig3_download_human_genes.py
│       ├── fig4_cluster_fasta_directory.sh
│       └── utils.R
├── data
│   ├── cox1.tre
```

## Methods

1. **fig1_analysis.ipynb**: Creates side-by-side fan plots of cox1 gene trees for plant (Streptophyta) and animal (Chordata) phylogenies with outlier removal.
2. **fig2_analysis.ipynb**: Downloads Ensembl Compara trees, calculates Hill's diversity across >14,000 protein families, and visualizes diversity patterns with branch length variance analysis.
3. **fig3_analysis.ipynb**: Integrates gene age data with evolutionary likelihood estimates using rolling window analysis to examine temporal trends across evolutionary time.
4. **fig4_analysis.ipynb**: Compares Hill's diversity of multiple sequence alignments with MMseqs clustering patterns across 218 Pfam families.
5. **fig5_analysis.ipynb**: Applies distance-based clustering to phylogenetic trees and calculates effective sequence ratios to quantify redundancy across >54,000 protein families.

## Compute Specifications
All analyses were done on an Apple MacBook Pro running macOS Montery with 32GB RAM, 10 cores, and 1TB of storage.

## Contributing
This is a research repository. For questions or collaboration inquiries, please contact the repository maintainers.
