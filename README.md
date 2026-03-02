# Morphodynamics of Human Early Brain Organoid Development

Reproducing Figure 1 from [Jain et al. 2025](https://www.nature.com/articles/s41586-025-09151-3) (Nature), which maps the transcriptomic landscape of human brain organoid development from Day 5 to Day 30. This work was done as part of my rotation in the Treutlein lab (ETH Zurich / QuadBio, Basel).

## What I'm reproducing

**Figure 1b-d**: Timecourse scRNA-seq integration across 6 timepoints (~41,000 cells), including:
- QC filtering, log-normalization, and cell cycle regression
- Integration with Cluster Similarity Spectrum (CSS)
- UMAP visualization with cluster annotation
- Dot plot analysis for cell type identification
- Stacked bar plot of cell type proportions over time

## Cell types identified

Neurectoderm, Late Neurectoderm, Prosencephalic progenitors, Telencephalic progenitors, Late Prosencephalic progenitors, Diencephalic progenitors, Tel/Die neurons

## Scripts

| File | Description |
|------|-------------|
| `fig1_morphodynamics_prescale.R` | Full pipeline: QC → normalization → CSS integration → clustering → annotation → UMAP |
| `fig1_morphodynamics.R` | Initial exploratory analysis |

## Tech Stack

R · Seurat · simspec (CSS) · ggplot2 · MetBrewer
