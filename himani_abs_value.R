library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)

# Load data (filtered data) & create Seurat object
HCD.data <- Read10X_h5("~/Downloads/filtered_feature_bc_matrix.h5")
HCD <- CreateSeuratObject(counts = HCD.data, project = "HCE048-D")

HCP.data <- Read10X_h5("~/Downloads/filtered_feature_bc_matrix (1).h5")
HCP <- CreateSeuratObject(counts = HCP.data, project = "HCE048-P")

IWD.data <- Read10X_h5("~/Downloads/filtered_feature_bc_matrix (2).h5")
IWD <- CreateSeuratObject(counts = IWD.data, project = "IW0076-D")

IWP.data <- Read10X_h5("~/Downloads/filtered_feature_bc_matrix (3).h5")
IWP <- CreateSeuratObject(counts = IWP.data, project = "IW0076-P")

P0D.data <- Read10X_h5("~/Downloads/filtered_feature_bc_matrix (4).h5")
P0D <- CreateSeuratObject(counts = P0D.data, project = "P0630-D")

P0P.data <- Read10X_h5("~/Downloads/filtered_feature_bc_matrix (5).h5")
P0P <- CreateSeuratObject(counts = P0P.data, project = "P0630-P")

# visualizing it with a violin plot

#first calculating the percent mt values
HCD[["percent.mt"]] <- PercentageFeatureSet(HCD, pattern = "^MT-")
HCP[["percent.mt"]] <- PercentageFeatureSet(HCP, pattern = "^MT-")
IWD[["percent.mt"]] <- PercentageFeatureSet(IWD, pattern = "^MT-")
IWP[["percent.mt"]] <- PercentageFeatureSet(IWP, pattern = "^MT-")
P0D[["percent.mt"]] <- PercentageFeatureSet(P0D, pattern = "^MT-")
P0P[["percent.mt"]] <- PercentageFeatureSet(P0P, pattern = "^MT-")


#merging the data & creating the violin plot
merged <- merge(HCD, y = c(HCP, IWD, IWP, P0D, P0P), add.cell.ids = c("HCD", "HCP", "IWD","IWP", "P0D", "P0P"))
VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0)