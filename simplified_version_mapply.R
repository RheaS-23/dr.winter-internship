#himani's ml 
library(Seurat)
library(dplyr)
library('umap')
umap.method="umap-learn"
library(ggplot2)
library(patchwork)

#reading the obj (my sample: IWD)
IWD <- readRDS("~/Downloads/IWD_postQC.rds")

#choosing diff number of variable genes
IWD_1 <-FindVariableFeatures(IWD, selection.method = "vst", nfeatures = 500)
IWD_2 <-FindVariableFeatures(IWD, selection.method = "vst", nfeatures = 1000)
IWD_3<-FindVariableFeatures(IWD, selection.method = "vst", nfeatures = 2000)

#creating the variable feature plot for each one
IWD_1_10 = head(VariableFeatures(IWD_1), 10)
IWD_2_10 = head(VariableFeatures(IWD_2), 10)
IWD_3_10 = head(VariableFeatures(IWD_3), 10)
#performing PCA & creating elbow plot for each one

IWD_1 <-NormalizeData(IWD_1)
IWD_1 <- ScaleData(IWD_1)
IWD_1 <- RunPCA(IWD_1)

IWD_2 <-NormalizeData(IWD_2)
IWD_2 <- ScaleData(IWD_2)
IWD_2 <- RunPCA(IWD_2)

IWD_3 <-NormalizeData(IWD_3)
IWD_3 <- ScaleData(IWD_3)
IWD_3 <- RunPCA(IWD_3)

#the functions

dim.list <- list(1:17, 1:18, 1:19, 1:20)
nnlist <- list(5L, 10L, 20L, 50L)
lrlist <- list(0.1, 0.5, 1, 1.5)
top10_list <- list(IWD_1_10, IWD_2_10, IWD_3_10)

umap.and.plots <- function(object, dims, nn, lr, top10) {
  object <- RunUMAP(object, dims = dims)
  dim_plot <- DimPlot(object, reduction = "umap", ncol = 2)+
    ggtitle("500 Variable Plot dimensions")
  feature_plot <- FeaturePlot(object, features = top10)
  print(dim_plot + feature_plot)
  object <- RunUMAP(object, dims = 1:18, n.neighbor = nn)
  dim_plot <- DimPlot(object, reduction = "umap", ncol = 2)+
    ggtitle("500 Variable Plot n.neighbor")
  feature_plot <- FeaturePlot(object, features = top10)
  print(dim_plot + feature_plot)
  object <- RunUMAP(object, dims = 1:18, learning.rate = lr)
  dim_plot <- DimPlot(object, reduction = "umap", ncol = 2)+
    ggtitle("500 Variable Plot learning rate")
  feature_plot <- FeaturePlot(object, features = top10)
  print(dim_plot + feature_plot)
}

mapply(FUN = umap.and.plots, object = list(IWD_1), dims = dim.list, nn = nnlist, lr = lrlist, top10 = list(IWD_1_10), SIMPLIFY = FALSE)

