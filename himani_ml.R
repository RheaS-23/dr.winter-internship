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

LabelPoints(VariableFeaturePlot(IWD_1), points = head(VariableFeatures(IWD_1), 10), repel = TRUE, xnudge = 0, ynudge = 0)
LabelPoints(VariableFeaturePlot(IWD_2), points = head(VariableFeatures(IWD_2), 10), repel = TRUE, xnudge = 0, ynudge = 0)
LabelPoints(VariableFeaturePlot(IWD_3), points = head(VariableFeatures(IWD_3), 10), repel = TRUE, xnudge = 0, ynudge = 0)

#performing PCA & creating elbow plot for each one

IWD_1 <-NormalizeData(IWD_1)
IWD_1 <- ScaleData(IWD_1)
IWD_1 <- RunPCA(IWD_1)
ElbowPlot(IWD_1) + labs(title = "500 Variable Features") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

IWD_2 <-NormalizeData(IWD_2)
IWD_2 <- ScaleData(IWD_2)
IWD_2 <- RunPCA(IWD_2)
ElbowPlot(IWD_2)+ labs(title = "1000 Variable Features") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

IWD_3 <-NormalizeData(IWD_3)
IWD_3 <- ScaleData(IWD_3)
IWD_3 <- RunPCA(IWD_3)
ElbowPlot(IWD_3)+ labs(title = "2000 Variable Features") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#generating UMAPS ----------------------------------------------

#changing #PCS

IWD_1 <- RunUMAP(IWD_1, dims = 1:17)
DimPlot(IWD_1, reduction = "umap") + labs(title = "Dims = 1:17") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_1 <- RunUMAP(IWD_1, dims = 1:18)
DimPlot(IWD_1, reduction = "umap")+ labs(title = "Dims = 1:18") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_1 <- RunUMAP(IWD_1, dims = 1:19)
DimPlot(IWD_1, reduction = "umap")+ labs(title = "Dims = 1:19") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_1 <- RunUMAP(IWD_1, dims = 1:20) 
DimPlot(IWD_1, reduction = "umap")+ labs(title = "Dims = 1:20") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#IWD-2 

IWD_2 <- RunUMAP(IWD_2, dims = 1:17)
DimPlot(IWD_2, reduction = "umap") + labs(title = "Dims = 1:17") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_2 <- RunUMAP(IWD_2, dims = 1:18)
DimPlot(IWD_2, reduction = "umap")+ labs(title = "Dims = 1:18") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_2 <- RunUMAP(IWD_2, dims = 1:19)
DimPlot(IWD_2, reduction = "umap")+ labs(title = "Dims = 1:19") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_2 <- RunUMAP(IWD_2, dims = 1:20) 
DimPlot(IWD_2, reduction = "umap")+ labs(title = "Dims = 1:20") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#IWD-3

IWD_3 <- RunUMAP(IWD_3, dims = 1:17)
DimPlot(IWD_3, reduction = "umap") + labs(title = "Dims = 1:17") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_3 <- RunUMAP(IWD_3, dims = 1:18)
DimPlot(IWD_3, reduction = "umap")+ labs(title = "Dims = 1:18") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_3 <- RunUMAP(IWD_3, dims = 1:19)
DimPlot(IWD_3, reduction = "umap")+ labs(title = "Dims = 1:19") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_3 <- RunUMAP(IWD_3, dims = 1:20) 
DimPlot(IWD_3, reduction = "umap")+ labs(title = "Dims = 1:20") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#changing the amt of n_neighbors

umap1 <- RunUMAP(IWD_1, dims = 1:18, n.neighbors = 5L)
umap2 <- RunUMAP(IWD_1, dims = 1:18, n.neighbors = 10L)
umap3 <- RunUMAP(IWD_1, dims = 1:18, n.neighbors = 20L)
umap4 <- RunUMAP(IWD_1, dims = 1:18, n.neighbors = 50L)

DimPlot(umap1, reduction = "umap") +
DimPlot(umap2, reduction = "umap")+
DimPlot(umap3, reduction = "umap")+
DimPlot(umap4, reduction = "umap")

#2

umap1 <- RunUMAP(IWD_2, dims = 1:18, n.neighbors = 5L)
umap2 <- RunUMAP(IWD_2, dims = 1:18, n.neighbors = 10L)
umap3 <- RunUMAP(IWD_2, dims = 1:18, n.neighbors = 20L)
umap4 <- RunUMAP(IWD_2, dims = 1:18, n.neighbors = 50L)

DimPlot(umap1, reduction = "umap") +
  DimPlot(umap2, reduction = "umap")+
  DimPlot(umap3, reduction = "umap")+
  DimPlot(umap4, reduction = "umap")

#3

umap1 <- RunUMAP(IWD_3, dims = 1:18, n.neighbors = 5L)
umap2 <- RunUMAP(IWD_3, dims = 1:18, n.neighbors = 10L)
umap3 <- RunUMAP(IWD_3, dims = 1:18, n.neighbors = 20L)
umap4 <- RunUMAP(IWD_3, dims = 1:18, n.neighbors = 50L)

DimPlot(umap1, reduction = "umap") +
  DimPlot(umap2, reduction = "umap")+
  DimPlot(umap3, reduction = "umap")+
  DimPlot(umap4, reduction = "umap")

#UMAPS by learning rate

umap1 <- RunUMAP(IWD_1, dims = 1:18, learning.rate = 0.1)
umap2 <- RunUMAP(IWD_1, dims = 1:18, learning.rate = 0.5)
umap3 <- RunUMAP(IWD_1, dims = 1:18, learning.rate = 1)
umap4 <- RunUMAP(IWD_1, dims = 1:18, learning.rate = 1.5)

DimPlot(umap1, reduction = "umap") +
  DimPlot(umap2, reduction = "umap")+
  DimPlot(umap3, reduction = "umap")+
  DimPlot(umap4, reduction = "umap")

umap1 <- RunUMAP(IWD_2, dims = 1:18, learning.rate = 0.1)
umap2 <- RunUMAP(IWD_2, dims = 1:18, learning.rate = 0.5)
umap3 <- RunUMAP(IWD_2, dims = 1:18, learning.rate = 1)
umap4 <- RunUMAP(IWD_2, dims = 1:18, learning.rate = 1.5)

DimPlot(umap1, reduction = "umap") +
  DimPlot(umap2, reduction = "umap")+
  DimPlot(umap3, reduction = "umap")+
  DimPlot(umap4, reduction = "umap")

umap1 <- RunUMAP(IWD_3, dims = 1:18, learning.rate = 0.1)
umap2 <- RunUMAP(IWD_3, dims = 1:18, learning.rate = 0.5)
umap3 <- RunUMAP(IWD_3, dims = 1:18, learning.rate = 1)
umap4 <- RunUMAP(IWD_3, dims = 1:18, learning.rate = 1.5)

DimPlot(umap1, reduction = "umap") +
  DimPlot(umap2, reduction = "umap")+
  DimPlot(umap3, reduction = "umap")+
  DimPlot(umap4, reduction = "umap")


#pca plots

VizDimLoadings(IWD_1, dims = 1:2, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
VizDimLoadings(IWD_2, dims = 1:2, reduction = "pca")+
VizDimLoadings(IWD_3, dims = 1:2, reduction = "pca")

VizDimLoadings(IWD_1, dims = 3:4, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 3:4, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 3:4, reduction = "pca")

VizDimLoadings(IWD_1, dims = 5:6, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 5:6, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 5:6, reduction = "pca")

VizDimLoadings(IWD_1, dims = 7:8, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 7:8, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 7:8, reduction = "pca")

VizDimLoadings(IWD_1, dims = 9:10, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 9:10, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 9:10, reduction = "pca")

VizDimLoadings(IWD_1, dims = 11:12, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 11:12, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 11:12, reduction = "pca")

VizDimLoadings(IWD_1, dims = 13:14, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 13:14, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 13:14, reduction = "pca")

VizDimLoadings(IWD_1, dims = 15:16, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 15:16, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 15:16, reduction = "pca")

VizDimLoadings(IWD_1, dims = 17:18, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 2000 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 17:18, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 17:18, reduction = "pca")

VizDimLoadings(IWD_1, dims = 19:20, reduction = "pca")+
  ggtitle("IW0076-D: 500, 1000, 500 Variable Genes")+
  VizDimLoadings(IWD_2, dims = 19:20, reduction = "pca")+
  VizDimLoadings(IWD_3, dims = 19:20, reduction = "pca")


#heatmaps
h1 <- DimHeatmap(IWD_1, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(IWD_1, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(IWD_1, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(IWD_1, dims = 19:20, cells = 500, balanced = TRUE)

h1 <- DimHeatmap(IWD_2, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(IWD_2, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(IWD_2, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(IWD_2, dims = 19:20, cells = 500, balanced = TRUE)


h1 <- DimHeatmap(IWD_3, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(IWD_3, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(IWD_3, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(IWD_3, dims = 19:20, cells = 500, balanced = TRUE)


#clustering ----------------------------------------------------------------------

#pcs for each sample changing k.param
clustered <- FindNeighbors(IWD_1, dims = 1:18, k.param = 5)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_1, dims = 1:18, k.param = 10)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_1, dims = 1:18, k.param = 20)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_1, dims = 1:18, k.param = 50)
clustered <- FindClusters(clustered, resolution = 0.5)

#2

clustered <- FindNeighbors(IWD_2, dims = 1:18, k.param = 5)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_2, dims = 1:18, k.param = 10)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_2, dims = 1:18, k.param = 20)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_2, dims = 1:18, k.param = 50)
clustered <- FindClusters(clustered, resolution = 0.5)

#3

clustered <- FindNeighbors(IWD_3, dims = 1:18, k.param = 5)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_3, dims = 1:18, k.param = 10)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_3, dims = 1:18, k.param = 20)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_3, dims = 1:18, k.param = 50)
clustered <- FindClusters(clustered, resolution = 0.5)

#changing resolution

clustered <- FindNeighbors(IWD_1, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.3)

clustered <- FindNeighbors(IWD_1, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_1, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 1)

clustered <- FindNeighbors(IWD_1, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 1.5)

#2

clustered <- FindNeighbors(IWD_2, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.3)

clustered <- FindNeighbors(IWD_2, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_2, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 1.0)

clustered <- FindNeighbors(IWD_2, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 1.5)

#3

clustered <- FindNeighbors(IWD_3, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.3)

clustered <- FindNeighbors(IWD_3, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_3, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 1)

clustered <- FindNeighbors(IWD_3, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 1.5)

#changing n.trees

clustered <- FindNeighbors(IWD_1, dims = 1:18, n.trees = 15)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_1, dims = 1:18, n.trees = 30)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_1, dims = 1:18, n.trees = 50)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_1, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5, n.trees = 65)

#2

clustered <- FindNeighbors(IWD_2, dims = 1:18, n.trees = 15)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_2, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5, n.trees = 30)

clustered <- FindNeighbors(IWD_2, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5, n.trees = 50)

clustered <- FindNeighbors(IWD_2, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5, n.trees = 65)

#3

clustered <- FindNeighbors(IWD_3, dims = 1:18, n.trees = 15)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_3, dims = 1:18, n.trees = 30)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_3, dims = 1:18, n.trees = 50)
clustered <- FindClusters(clustered, resolution = 0.5)

clustered <- FindNeighbors(IWD_3, dims = 1:18, n.trees = 65)
clustered <- FindClusters(clustered, resolution = 0.5)


dim.list <- list(1:17, 1:18, 1:19, 1:20)
top10_list <- list(IWD_1_10, IWD_2_10, IWD_3_10)

umap.and.plots <- function(object, dims, top10) {
  object <- RunUMAP(object, dims = dims)
  dim_plot <- DimPlot(object, reduction = "umap", ncol = 2)
  feature_plot <- FeaturePlot(object, features = top10)
  print(dim_plot)
  print(feature_plot)
  return(list(dim_plot, feature_plot))
}

mapply(FUN = umap.and.plots, object = list(IWD_1), dims = dim.list, top10 = top10_list, SIMPLIFY = FALSE)

