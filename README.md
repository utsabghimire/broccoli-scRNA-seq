# 🥦 Broccoli Single-Cell RNA-seq Analysis Pipeline

This repository contains a complete pipeline for analyzing single-cell RNA-seq (scRNA-seq) data in **broccoli (*Brassica oleracea*)** using the **Seurat** package for clustering and **Monocle 3** for pseudotime/trajectory analysis.

---

## 📌 Project Goal

The aim is to investigate gene expression heterogeneity and senescence progression in broccoli at the **single-cell level**, particularly focusing on preharvest and postharvest conditions. This pipeline enables:

- Quality control and filtering  
- Normalization and dimensionality reduction  
- Cell clustering and UMAP visualization  
- Marker gene identification  
- Pseudotime trajectory analysis  

---

## 🧪 Dataset

The input dataset is generated using **10X Genomics** single-nucleus RNA-seq (snRNA-seq) on broccoli inflorescence tissue. It includes the following files typically found in the `filtered_feature_bc_matrix` directory:

```
matrix.mtx  
barcodes.tsv.gz  
features.tsv.gz  
```

---

## 🔧 Requirements

Install the required R packages:

```r
install.packages("Matrix")
install.packages("Seurat")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("monocle3")
```

---

## 🚀 How to Run

Clone this repository:

```bash
git clone https://github.com/your-username/broccoli-scRNA-seq.git
cd broccoli-scRNA-seq
```

Open and run the pipeline in R:

```r
source("broccoli_seurat_pipeline.R")
```

Make sure to update the `setwd()` path in the script to where your data files are located.

---

## 📊 Output Files

- `cluster_markers.csv`: Top marker genes per cluster  
- `pseudotime_values.csv`: Pseudotime values for each cell  

---

## 📈 Visualizations Generated

- Violin plots of QC metrics (`nFeature_RNA`, `nCount_RNA`)  
- PCA and Elbow plots  
- UMAP with clusters  
- Top variable features plot  
- Pseudotime trajectory plot  

---

## 🧬 Biological Application

This pipeline is designed to identify **senescence-associated genes (SAGs)** and understand **cell type-specific responses** during postharvest senescence. By integrating pseudotime analysis, we aim to reconstruct the **senescence trajectory** and reveal regulatory mechanisms underlying **broccoli quality decline after harvest**.
