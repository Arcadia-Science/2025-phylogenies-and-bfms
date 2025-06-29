# Phylogenies and Biological Foundation Models

## Purpose
This repo contains code associated with the pub "Phylogenies and Biological Foundation Models".

## Installation and Setup
This repository uses R and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html). 

## Data

### Overview
The repository contains code to collect the necessary data to reproduce analyses in the pub. Example Pfam protein families used in Figure 4 are available on [Zenodo](https://zenodo.org/records/15644457).

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
