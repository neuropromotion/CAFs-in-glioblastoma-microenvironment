# CAFs in the Glioblastoma Microenvironment

This repository contains code for an end-to-end **single-cell RNA-seq analysis** of **cancer-associated fibroblasts (CAFs)** and other stromal cell populations within the **glioblastoma (GBM) tumor microenvironment**.  
The workflow highlights the **multifaceted role of stromal cells** in GBM, including:

- Support of **glioblastoma stem cell (GSC)** proliferation  
- **Pro-angiogenic** remodeling and vascular niche maintenance  
- **Regional immunosuppression** and immune evasion  
- Crosstalk between CAFs, tumor cells, endothelial cells, and infiltrating immune cells

---

## Repository Structure

| File | Description |
|------|-------------|
| `main_branch.R` | Processing workflow for the main dataset |
| `Validation_1.R` | Processing workflow for validation dataset #1 |
| `Validation_2.R` | Processing workflow for validation dataset #2 |
| `GSVA.R` | **GSVA-based subtype annotation** of glioblastoma tumor cells. Requires `GBM_subtypes_DEG.xlsx` |
| `GBM_subtypes_DEG.xlsx` | DEG signatures for glioblastoma transcriptional subtypes (PMIDs: 31327527, 40191211) |
| `cellchat.R` | **CellChat** pipeline for ligand–receptor interaction inference. Uses configuration in `cellchat_configue/` |
| `DEG.R` | Differential expression analysis (cluster-wise or group-wise) |
| `save_plots.R` | Script for figure post-processing and export |

---

## Analysis Pipeline Overview

### 1. **Ambient RNA Correction**
- **DecontX** — applied to the main dataset and validation dataset #1  
- **SoupX** — applied to validation dataset #2  

### 2. **Doublet Detection**
- Performed using **DoubletFinder**

### 3. **Batch Integration**
- **Harmony** used for removing batch effects (main dataset)

### 4. **Core scRNA-seq Workflow**
- Normalization  
- Scaling  
- PCA  
- UMAP projection  
- Clustering using the **Louvain algorithm**

### 5. **Tumor Cell Subtype Annotation**
- **GSVA** used to classify GBM tumor cells into six transcriptional states:  
  - **Glioblastoma Stem Cells**  
  - **MES1-like**, **MES2-like**  
  - **NPC1-like**, **NPC2-like**  
  - **AC-like**  
  - **OPC-like**

### 6. **Cell–Cell Communication Analysis**
- Conducted using **CellChat**, focusing on:
  - CAF–tumor signaling
  - CAF–endothelium interactions
  - CAF-immune interactions 
  - stromal-tumor interactions
  - stromal-immune interactions

---

## Contact
Work email: amismailov@hse.ru
Personal email: neuro.promotion@gmail.com 

---
