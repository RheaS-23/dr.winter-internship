library(dplyr)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(hdf5r)
library(cowplot)

# Load data
HCD <- readRDS("~/Desktop/RStudio/SingleCellProject/SingleCells/Week3_PostQC/HCD_postQC.rds")

# Normalize Data
HCD <- NormalizeData(HCD)

# Different # of Variable Genes
HCD1 <- FindVariableFeatures(HCD, selection.method = "vst", nfeatures = 500)
top10_1 <- head(VariableFeatures(HCD1), 10)
HCD1_plot <-LabelPoints(VariableFeaturePlot(HCD1), points = top10_1, repel = TRUE)

HCD2 <- FindVariableFeatures(HCD, selection.method = "vst", nfeatures = 1000)
top10_2 <- head(VariableFeatures(HCD2), 10)
HCD2_plot <-LabelPoints(VariableFeaturePlot(HCD2), points = top10_2, repel = TRUE)

HCD3 <- FindVariableFeatures(HCD, selection.method = "vst", nfeatures = 2000)
top10_3 <- head(VariableFeatures(HCD3), 10)
HCD3_plot <-LabelPoints(VariableFeaturePlot(HCD3), points = top10_3, repel = TRUE)

HCD1_plot + HCD2_plot + HCD3_plot

# Scale Data
HCD1 <- ScaleData(HCD1)
HCD2 <- ScaleData(HCD2)
HCD3 <- ScaleData(HCD3)

# Performing PCA & Elbow Plot
HCD1 <- RunPCA(HCD1, features = VariableFeatures(object = HCD1))
HCD2 <- RunPCA(HCD2, features = VariableFeatures(object = HCD2))
HCD3 <- RunPCA(HCD3, features = VariableFeatures(object = HCD3))
ElbowPlot(HCD1) + ElbowPlot(HCD2) + ElbowPlot(HCD3)


dim.list <- list(1:17, 1:18, 1:19, 1:20)
top10_list <- list(top10_1, top10_2, top10_3)


umap.and.plots <- function(object, dims, top10) {
  object <- RunUMAP(object, dims = dims)
  dim_plot <- DimPlot(object, reduction = "umap", ncol = 2)
  feature_plot <- FeaturePlot(object, features = top10)
  print(dim_plot)
  print(feature_plot)
  return(list(dim_plot, feature_plot))
}

mapply(FUN = umap.and.plots, object = list(HCD1), dims = dim.list, top10 = top10_list, SIMPLIFY = FALSE)



