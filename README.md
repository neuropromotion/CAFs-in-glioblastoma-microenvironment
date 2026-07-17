# CAFs in the Glioblastoma Microenvironment

Overview

This repository contains code for an integrated analysis of single-cell and bulk RNA-seq data aimed at characterizing cancer-associated fibroblasts (CAFs) in the glioblastoma (GBM) tumor microenvironment.

Glioblastoma is the most aggressive primary brain tumor in adults, characterized by profound remodeling of its tumor microenvironment. However, the existence, origin, and functional relevance of CAFs in glioblastoma remain incompletely understood.

In this study, we analyze 54 single-cell RNA-seq datasets and 88 bulk RNA-seq samples, integrating independent discovery and validation cohorts to systematically characterize CAFs in GBM. Our analysis identifies a reproducible transcriptional continuum linking endothelial cells, pericytes, and CAFs, and defines a robust cross-cohort CAF gene signature. We further investigate CAF-associated cell–cell communication programs and quantify CAF abundance in bulk tumors using both signature-based scoring and deconvolution approaches.

---
## Original scRNA-seq dataset: discovery - [GSE173278]; validation - GSM3828672 + GSE103224 + GSE135045 + 10x Genomics GBM 3' v3 
## https://www.kaggle.com/datasets/ismailovaly/caf-in-glioblastoma
---
## Repository Structure

| File | Description |
|------|-------------|
| `load_data.R` | Data loading and processing workflow for the discovery and validation datasets |
| `functions.R` | auxillary functions |
| `cellchat.R` | **CellChat** pipeline for ligand–receptor interaction inference |
| `GAM.R` | Processing workflow for generative additive models implementation |
| `ESCAPE.R` | Workflow for functional enrichment analysis |
| `Diffusion_maps.R` | Workflow for diffusion maps and time trajectory analysis |
| `CNA.R` | InferCNA - copy number alteration inferring framework |
| `for_deconvolution.Rmd` | Script for data preparation for bulk-deconvolution analysis |
| `gse190504_CAF_rank_score.py` | Python script for CAF signature rank expression evaluation across bulk gse190504 data conditions |
| `gse190504_cellanneal_deconvolution.py` | Python script for cellanneal bulk-deconvolution analysis within bulk gse190504 data |

---


---

## Contact
Work email: amismailov@hse.ru
Personal email: ismailov.aly@gmail.com

---
