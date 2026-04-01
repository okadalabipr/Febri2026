# 01_preprocessing.R
# Description:
# Preprocessing of scRNA-seq data including normalization, dimensionality reduction, clustering, and GO analysis.

# Load libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

set.seed(1)

# Define paths
data_path <- "data/allsample.RDS"
output_dir <- "results/preprocessing/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
allsample <- readRDS(data_path)

# Metadata decoding
allsample$orig.ident <- recode(allsample$orig.ident,
                               "EGF_4H"   = "EGF_4h",
                               "EGF_8H"   = "EGF_8h",
                               "EGF_16H"  = "EGF_16h",
                               "HRG_4H"   = "HRG_4h",
                               "HRG_8H"   = "HRG_8h",
                               "HRG_16H"  = "HRG_16h")


# Normalization
DefaultAssay(allsample) <- "RNA"
allsample <- JoinLayers(allsample)

allsample <- NormalizeData(allsample)
allsample <- FindVariableFeatures(allsample, nfeatures = 2000)

# Scaling & PCA
allsample <- ScaleData(allsample)
allsample <- RunPCA(allsample)

# Dimensionality reduction
allsample <- RunUMAP(allsample, dims = 1:20)
allsample <- RunTSNE(allsample, dims = 1:20)

# Visualization
color_palette <- c(
  "EGF_4h"  = "#DEEBF7",
  "EGF_8h"  = "#9ECAE1",
  "EGF_16h" = "blue",
  "HRG_4h"  = "#FEE5D9",
  "HRG_8h"  = "#FCAE91",
  "HRG_16h" = "red"
)

p_umap <- DimPlot(allsample, group.by = "orig.ident") +
  scale_color_manual(values = color_palette)

ggsave(file.path(output_dir, "UMAP.png"), p_umap)

# Clustering
allsample <- FindNeighbors(allsample, dims = 1:20)
allsample <- FindClusters(allsample, resolution = 0.05)

# GO analysis
markers <- FindAllMarkers(allsample, only.pos = TRUE)

marker_list <- split(markers$gene, markers$cluster)

marker_list_entrez <- lapply(marker_list, function(genes) {
  bitr(genes, fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)$ENTREZID
})

go_results <- lapply(marker_list_entrez, function(gene_ids) {
  enrichGO(gene = gene_ids,
           OrgDb = org.Hs.eg.db,
           ont = "BP")
})

# Save objects
saveRDS(allsample, file = file.path(output_dir, "processed_seurat.rds"))