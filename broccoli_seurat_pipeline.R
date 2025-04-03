# ====================================================
# 🥦 Broccoli Single-Cell RNA-seq Analysis Pipeline
# Author: Utsab Ghimire
# ====================================================

# Install required packages (run only once)
# install.packages("Matrix")
# install.packages("Seurat")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("monocle3")

# Load libraries
library(Matrix)
library(Seurat)
library(monocle3)

# Set working directory (edit this to your local path)
setwd("/Users/utsabghimire/Downloads/Pre1/outs/filtered_feature_bc_matrix")

# Load 10X data files
expression_matrix <- readMM("matrix.mtx")
barcodes <- readLines("barcodes.tsv.gz")
features <- read.delim("features.tsv.gz", header = FALSE)

# Assign row and column names
rownames(expression_matrix) <- features[, 2]
colnames(expression_matrix) <- barcodes

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = expression_matrix, project = "Broccoli_Senescence")

# Quality control visualization
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Filter cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)
VariableFeaturePlot(seurat_obj) + LabelPoints(points = top10, repel = TRUE)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
ElbowPlot(seurat_obj)

# Clustering and UMAP
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Marker genes
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, "cluster_markers.csv", row.names = FALSE)

# ----------------------------------------------------
# 📈 Pseudotime Analysis with Monocle 3
# ----------------------------------------------------

# Convert to Monocle object
cds <- as.cell_data_set(seurat_obj)

# Preprocessing & dimensionality reduction
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")

# Clustering and trajectory
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)

# Plot pseudotime
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE, label_leaves = FALSE)

# Save pseudotime values
pseudotime <- pseudotime(cds)
write.csv(as.data.frame(pseudotime), "pseudotime_values.csv")
