# CAFs-in-glioblastoma-microenvironment
Code for end‑to‑end scRNA‑seq analysis of cancer‑associated fibroblasts (CAFs) and other stromal compartments in the glioblastoma microenvironment. Analysis reveals the multifaceted impact of stromal cells on GBM biology, including support of glioblastoam stem cells proliferation, pro‑angiogenic signaling, and regional immunosuppression

---

## Files

| File | Description |
|------|-------------|
| `main_branch.R` | Processing of main dataset |
| `Validation_1.R` | Processing of validation dataset №1 |
| `Validation_2.R` | Processing of validation dataset №2 |
| `GSVA.R` | GSVA-based annotation of glioblastoma cells (subtypes revealing). Requires `GBM_subtypes_DEG.xlsx` |
| `GBM_subtypes_DEG.xlsx` | Differentially expressed genes of glioblastoma subtypes (PMIDs: 31327527, 40191211) |
| `cellchat.R` | CellChat pipeline. Requires files stored in `cellchat_configue` |
| `DEG.R` | Differentially expressed genes evaluation |
| `save_plots.R` | Figures precessing |

---
