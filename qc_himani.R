
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
cells<- VlnPlot(merged, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident")
genes <- VlnPlot(merged, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident")
pmt <- VlnPlot(merged, features = "percent.mt", pt.size = 0, group.by = "orig.ident")
cells + geom_hline(aes(yintercept = 3800), color = "black")
genes + geom_hline(aes(yintercept = 1100), color = "black")
pmt + geom_hline(aes(yintercept = 2), color = "black")+geom_hline(aes(yintercept = 20), color = "black")

#using absolute value QC strategy (identifying & removing low-quality cells)

# Define cutoff values

count_cutoff <- 3800
gene_cutoff <- 1100
percent_cutoff <- 2
percent_upper_cutoff <- 20

# Filter cells based on cutoffs
HCD <- subset(HCD, subset = nCount_RNA > count_cutoff)
HCP <- subset(HCP, subset = nCount_RNA > count_cutoff)
IWD <- subset(IWD, subset = nCount_RNA > count_cutoff)
IWP <- subset(IWP, subset = nCount_RNA > count_cutoff)
P0D <- subset(P0D, subset = nCount_RNA > count_cutoff)
P0P <- subset(P0P, subset = nCount_RNA > count_cutoff)

# Filter cells by gene
HCD <- subset(HCD, subset = nFeature_RNA > gene_cutoff)
HCP <- subset(HCP, subset = nFeature_RNA > gene_cutoff)
IWD <- subset(IWD, subset = nFeature_RNA > gene_cutoff)
IWP <- subset(IWP, subset = nFeature_RNA > gene_cutoff)
P0D <- subset(P0D, subset = nFeature_RNA > gene_cutoff)
P0P <- subset(P0P, subset = nFeature_RNA > gene_cutoff)

# Filter cells by percent
HCD <- subset(HCD, subset = percent.mt > percent_cutoff)
HCP <- subset(HCP, subset =percent.mt > percent_cutoff)
IWD <- subset(IWD, subset = percent.mt > percent_cutoff)
IWP <- subset(IWP, subset = percent.mt > percent_cutoff)
P0D <- subset(P0D, subset = percent.mt > percent_cutoff)
P0P <- subset(P0P, subset = percent.mt > percent_cutoff)
HCD <- subset(HCD, subset = percent.mt < percent_upper_cutoff)
HCP <- subset(HCP, subset = percent.mt < percent_upper_cutoff)
IWD <- subset(IWD, subset = percent.mt < percent_upper_cutoff)
IWP <- subset(IWP, subset = percent.mt < percent_upper_cutoff)
P0D <- subset(P0D, subset = percent.mt < percent_upper_cutoff)
P0P <- subset(P0P, subset = percent.mt < percent_upper_cutoff)

#updating the seurat object
HCD <- UpdateSeuratObject(HCD)
HCP <- UpdateSeuratObject(HCP)
IWD <- UpdateSeuratObject(IWD)
IWP <- UpdateSeuratObject(IWP)
P0D <- UpdateSeuratObject(P0D)
P0P <- UpdateSeuratObject(P0P)

#visualizing it after the changes using a violin plot

#first calculating the percent mt values
HCD[["percent.mt"]] <- PercentageFeatureSet(HCD, pattern = "^MT-")
HCP[["percent.mt"]] <- PercentageFeatureSet(HCP, pattern = "^MT-")
IWD[["percent.mt"]] <- PercentageFeatureSet(IWD, pattern = "^MT-")
IWP[["percent.mt"]] <- PercentageFeatureSet(IWP, pattern = "^MT-")
P0D[["percent.mt"]] <- PercentageFeatureSet(P0D, pattern = "^MT-")
P0P[["percent.mt"]] <- PercentageFeatureSet(P0P, pattern = "^MT-")

#merging the data & creating the violin plot
merged <- merge(HCD, y = c(HCP, IWD, IWP, P0D, P0P), add.cell.ids = c("HCD", "HCP", "IWD","IWP", "P0D", "P0P"))
cells<- VlnPlot(merged, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident")
genes <- VlnPlot(merged, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident")
pmt <- VlnPlot(merged, features = "percent.mt", pt.size = 0, group.by = "orig.ident")
cells + geom_hline(aes(yintercept = 3800), color = "black")
genes + geom_hline(aes(yintercept = 1100), color = "black")
pmt + geom_hline(aes(yintercept = 2), color = "black")+geom_hline(aes(yintercept = 20), color = "black")

