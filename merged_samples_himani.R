#himani's merged post QC clustering
library(Seurat)
library(dplyr)
library('umap')
umap.method="umap-learn"
library(ggplot2)
library(patchwork)

#reading the objs
HCD <- readRDS("~/Downloads/HCD_postQC.rds")
HCP <- readRDS("~/Downloads/HCP_postQC.rds")
IWD <- readRDS("~/Downloads/IWD_postQC.rds")
IWP <- readRDS("~/Downloads/IWP_postQC.rds")
P0D <- readRDS("~/Downloads/P0D_postQC.rds")
P0P <- readRDS("~/Downloads/P0P_postQC.rds")

#merging everything together
merged <- merge(HCD, y = c(HCP, IWD, IWP, P0D, P0P), add.cell.ids = c("HCD", "HCP", "IWD","IWP", "P0D", "P0P"))
merged <-FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <-NormalizeData(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

#umap before
umap <- RunUMAP(merged, dims = 1:18, n.neighbors = 50)
dim_plot <- DimPlot(umap, reduction = "umap", ncol = 2)+
  ggtitle("Merged Samples UMAP")
top_genes <- list('KRT6A', "KRT14", "MUC5B", "PTPRC", "CD3D", "LYZ", "HLA-DRA",
                  "TPSAB1", "PECAM1", "PDGFRA", "PDGFRB", "RGS5", "TAGLN", "ACTA2")
feature_plot <- FeaturePlot(umap, features = top_genes)
dim_plot+feature_plot

#clustering!!
clustered <- FindNeighbors(merged, dims = 1:18)
clustered <- FindClusters(clustered, resolution = 0.5)
umap <- RunUMAP(clustered, dims = 1:18)
plot2 <- DimPlot(umap, reduction = "umap", ncol = 2, label = TRUE)+
  ggtitle("Clustered UMAP")
feature_plot <- FeaturePlot(umap, features = top_genes)
print(plot2 + feature_plot)

#renaming idents
epi_pos_markers <- c('ATP1B3','KRT15','KRT5', 'KRT13','CD24','CNFN')
epi_neg_markers <- c('VIM','PTPRC','PDGFRA', 'VWF','CD3D','HLA-DPA1','CPA3','CLDN5','COL3A1')
VlnPlot(clustered, features = epi_pos_markers, pt.size = 0)
VlnPlot(clustered, features = epi_neg_markers, pt.size = 0)
clustered <- SetIdent(clustered, value = "seurat_clusters")
clustered <- RenameIdents(clustered, "11" = "Myeloid",
                          "15" = "Fibroblasts",
                          "14" = "Lymphocytes",
                          "13" = "Endothelial",
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
                          "12" = "Epithelial")

umap <- RunUMAP(clustered, dims = 1:18)
DimPlot(umap, reduction = "umap", ncol = 2, label = TRUE)+
  ggtitle("Merged UMAP (Annotated)")


#finding key genes
tmo.markers <- FindAllMarkers(clustered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tmo.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
View(tmo.markers)

#top genes by clusters (cluster markers)
tmo.markers%>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top10
VlnPlot(clustered, features = top10$gene, pt.size = 0)

#top two genes by cluster:
#Myeloid: HLA-DRA, HLA-DPA1
#Fibroblasts: DCN, MGP
#Lymphocytes: IL32, S100A4
#Endothelial: IGFBP7, SPARCL1
#Epithelial: KRT13, SPRR3
