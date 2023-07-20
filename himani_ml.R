#himani's ml 
library(Seurat)
library(dplyr)
library('umap')
umap.method="umap-learn"

#reading the obj (my sample: IWD)
IWD <- readRDS("~/Downloads/IWD1.rds")

#choosing diff number of variable genes
IWD_1 <-FindVariableFeatures(IWD, selection.method = "vst", nfeatures = 500)
IWD_2 <-FindVariableFeatures(IWD, selection.method = "vst", nfeatures = 1000)
IWD_3<-FindVariableFeatures(IWD, selection.method = "vst", nfeatures = 2000)

#creating the variable feature plot for each one
VariableFeaturePlot(IWD_1)
VariableFeaturePlot(IWD_2)
VariableFeaturePlot(IWD_3)

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

IWD_1 <- RunUMAP(IWD_1, dims = 1:5)
DimPlot(IWD_1, reduction = "umap") + labs(title = "Dims = 1:5") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_1 <- RunUMAP(IWD_1, dims = 1:10)
DimPlot(IWD_1, reduction = "umap")+ labs(title = "Dims = 1:10") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_1 <- RunUMAP(IWD_1, dims = 1:15)
DimPlot(IWD_1, reduction = "umap")+ labs(title = "Dims = 1:15") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
IWD_1 <- RunUMAP(IWD_1, dims = 1:20) 
DimPlot(IWD_1, reduction = "umap")+ labs(title = "Dims = 1:20") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
