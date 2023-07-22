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
top_genes <- list('KRT6A', "KRT14", "MUC5B", "PTPRC", "CD3D", "LYZ", "HLA-DRA",
                  "TPSAB1", "PECAM1", "PDGFRA", "PDGFRB", "RGS5", "TAGLN", "ACTA2")

umap.and.plots <- function(object, dims, nn, lr, top10) {
  object <- RunUMAP(object, dims = 1:18, learning.rate = lr)
  dim_plot <- DimPlot(object, reduction = "umap", ncol = 2)+
    ggtitle("500 Variable Plot dimensions")
  feature_plot <- FeaturePlot(object, features = top10)
  print(dim_plot + feature_plot)
}

#mapply(FUN = umap.and.plots, object = list(IWD_3), dims = dim.list, nn = nnlist, lr = lrlist, top10 = list(top_genes), SIMPLIFY = FALSE)

#top 3 umaps: 
umap1 <- RunUMAP(IWD_1, dims = 1:18, n.neighbors = 50)
dim_plot1 <- DimPlot(umap1, reduction = "umap", ncol = 2)+
  ggtitle("500 VF Nearest Neighbors = 20")
feature_plot <- FeaturePlot(umap1, features = top_genes)
#print(dim_plot + feature_plot)

umap2 <- RunUMAP(IWD_2, dims = 1:20)
dim_plot2 <- DimPlot(umap2, reduction = "umap", ncol = 2)+
  ggtitle("1000 VF Dims = 20")
feature_plot <- FeaturePlot(umap2, features = top_genes)
#print(dim_plot + feature_plot)

umap3 <- RunUMAP(IWD_3, dims = 1:18, n.neighbors = 50)
dim_plot3 <- DimPlot(umap3, reduction = "umap", ncol = 2)+
  ggtitle("2000 VF Nearest Neighbors = 20")
feature_plot <- FeaturePlot(umap3, features = top_genes)
#print(dim_plot + feature_plot)

#clustering ----------------------------------------------------------------------

#top 3 clusterings 

clustered <- FindNeighbors(umap1, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5)
umap <- RunUMAP(clustered, dims = 1:18)
plot2 <- DimPlot(umap, reduction = "umap", ncol = 2)+
  ggtitle("UMAP1 resolution = 0.5")

clustered <- FindNeighbors(umap2, dims = 1:18, n.trees = 30)
clustered <- FindClusters(clustered, resolution = 0.5)
umap <- RunUMAP(clustered, dims = 1:18)
plot3 <- DimPlot(umap, reduction = "umap", ncol = 2)+
  ggtitle("UMAP2 n.trees = 30")

clustered <- FindNeighbors(umap3, dims = 1:18, k.param = 5)
clustered <- FindClusters(clustered, resolution = 0.5)
umap <- RunUMAP(clustered, dims = 1:18)
plot1 <- DimPlot(umap, reduction = "umap", ncol = 2)+
  ggtitle("UMAP3 k.param = 5")

plot2+plot3+plot1


#top clustering

clustered <- FindNeighbors(umap1, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5)
umap <- RunUMAP(clustered, dims = 1:18)
plot2 <- DimPlot(umap, reduction = "umap", ncol = 2)+
  ggtitle("UMAP1 resolution = 0.5")
feature_plot <- FeaturePlot(umap, features = top_genes)
print(plot2 + feature_plot)

epi_pos_markers <- c('ATP1B3','KRT15','KRT5', 'KRT13','CD24','CNFN')
epi_neg_markers <- c('VIM','PTPRC','PDGFRA', 'VWF','CD3D','HLA-DPA1','CPA3','CLDN5','COL3A1')
VlnPlot(clustered, features = epi_pos_markers, pt.size = 0)
VlnPlot(clustered, features = epi_neg_markers, pt.size = 0)
clustered <- SetIdent(clustered, value = "seurat_clusters")
clustered <- RenameIdents(clustered, "13" = "Lymphocytes",
                          "11" = "Myeloid", 
                          "0" = "Epithelial",
                          "1" = "Epithelial",
                          "2" = "Epithelial",
                          "3" = "Epithelial",
                          "4" = "Epithelial",
                          "5" = "Epithelial",
                          "6" = "Epithelial",
                          "7" = "Epithelial",
                          "8" = "Epithelial",
                          "9" = "Epithelial",
                          "10" = "Epithelial",
                          "12" = "Fibroblasts")

umap <- RunUMAP(clustered, dims = 1:18)
plot2 <- DimPlot(umap, reduction = "umap", ncol = 2, label = TRUE)+
  ggtitle("UMAP1 resolution = 0.5")
feature_plot <- FeaturePlot(clustered, features = top_genes)
print(plot2 + feature_plot)
VlnPlot(clustered, features = top_genes, pt.size = 0, group.by = "seurat_clusters")

#finding key genes
tmo.markers <- FindAllMarkers(clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tmo.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
View(tmo.markers)

#top genes by cluster:
  #lympocytes: CCL5, IL7R
  #Fibroblasts: IGFBP7, AQP1
  #Myeloid: HLA-DRA, HLA-DPA1
  #Epithelial: PERP, S100A14
top_cluster <- list("CCL5", "IL7R", "IGFBP7", "AQP1", "HLA-DRA", "HLA-DPA1", "PERP", "S100A14")
VlnPlot(clustered, features = top_cluster, pt.size = 0)

DimHeatmap(clustered, cells = 500, balanced = TRUE)
DimHeatmap(IWD_1, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(IWD_1, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(IWD_1, dims = 19:20, cells = 500, balanced = TRUE)
