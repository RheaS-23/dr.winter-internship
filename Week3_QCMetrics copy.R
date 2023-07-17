library(dplyr)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(hdf5r)
library(SoupX)
library(ggplot2)
library(cowplot)

# Load data (filtered data) & create Seurat object
HCD.data <- Read10X("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/OneDrive_2_7-10-2023/HCE048-D_MP_cellranger6_01182023/outs/filtered_feature_bc_matrix/")
HCD <- CreateSeuratObject(counts = HCD.data, project = "HCE048-D")

HCP.data <- Read10X("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/OneDrive_2_7-10-2023/HCE048-P_MP_cellranger6_01182023/outs/filtered_feature_bc_matrix/")
HCP <- CreateSeuratObject(counts = HCP.data, project = "HCE048-P")

IWD.data <- Read10X("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/OneDrive_1_7-10-2023/IW0076-D_MP_cellranger6_01182023/outs/filtered_feature_bc_matrix/")
IWD <- CreateSeuratObject(counts = IWD.data, project = "IW0076-D")

IWP.data <- Read10X("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/OneDrive_1_7-10-2023/IW0076-P_MP_cellranger6_01182023/outs/filtered_feature_bc_matrix/")
IWP <- CreateSeuratObject(counts = IWP.data, project = "IW0076-P")

P0D.data <- Read10X("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/P0630-D filtered_feature_bc_matrix/")
P0D <- CreateSeuratObject(counts = P0D.data, project = "P0630-D")

P0P.data <- Read10X_h5("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/P0630-P filtered_feature_bc_matrix.h5")
P0P <- CreateSeuratObject(counts = P0P.data, project = "P0630-P")

# Load raw counts matrix
HCD_raw <- Read10X_h5("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/HCD_raw_feature_bc_matrix.h5")
HCP_raw <- Read10X_h5("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/HCP_raw_feature_bc_matrix.h5")
IWD_raw <- Read10X_h5("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/IWD_raw_feature_bc_matrix.h5")
IWP_raw <- Read10X_h5("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/IWP_raw_feature_bc_matrix.h5")
P0D_raw <- Read10X_h5("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/P0D_raw_feature_bc_matrix.h5")
P0P_raw <- Read10X_h5("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/P0P_raw_feature_bc_matrix.h5")


# Calculating percentage of mitochondrial genes per cell
HCD[["percent.mt"]] <- PercentageFeatureSet(HCD, pattern = "^MT-")
HCP[["percent.mt"]] <- PercentageFeatureSet(HCP, pattern = "^MT-")
IWD[["percent.mt"]] <- PercentageFeatureSet(IWD, pattern = "^MT-")
IWP[["percent.mt"]] <- PercentageFeatureSet(IWP, pattern = "^MT-")
P0D[["percent.mt"]] <- PercentageFeatureSet(P0D, pattern = "^MT-")
P0P[["percent.mt"]] <- PercentageFeatureSet(P0P, pattern = "^MT-")

# Calculating percentage of ribosomal genes per cell
HCD[["percent.rb"]] <- PercentageFeatureSet(HCD, pattern = "^RP[SL]")
HCP[["percent.rb"]] <- PercentageFeatureSet(HCP, pattern = "^RP[SL]")
IWD[["percent.rb"]] <- PercentageFeatureSet(IWD, pattern = "^RP[SL]")
IWP[["percent.rb"]] <- PercentageFeatureSet(IWP, pattern = "^RP[SL]")
P0D[["percent.rb"]] <- PercentageFeatureSet(P0D, pattern = "^RP[SL]")
P0P[["percent.rb"]] <- PercentageFeatureSet(P0P, pattern = "^RP[SL]")

med_calc <- function(obj) {
  # Calculating median value of each metric
  count_med <- median(obj@meta.data$nCount_RNA)
  feature_med <- median(obj@meta.data$nFeature_RNA)
  pct_mt_med <- median(obj@meta.data$percent.mt)
  pct_rb_med <- median(obj@meta.data$percent.rb)
  
  return(c(count_med, feature_med, pct_mt_med, pct_rb_med))
}

mad_calc <- function(obj) {
  count_med <- med_calc(obj)[1]
  feature_med <- med_calc(obj)[2]
  pct_mt_med <- med_calc(obj)[3]
  pct_rb_med <- med_calc(obj)[4]
  
  # Calculates absolute deviation (subtract median from each value)
  nCt_abs_dev <- abs(obj@meta.data$nCount_RNA - count_med)
  nFt_abs_dev <- abs(obj@meta.data$nFeature_RNA - feature_med)
  mt_abs_dev <- abs(obj@meta.data$percent.mt - pct_mt_med)
  rb_abs_dev <- abs(obj@meta.data$percent.rb - pct_rb_med)
  
  
  # Calculates the median absolute deviation
  nCt_mad <- mad(nCt_abs_dev, constant = 1)
  nFt_mad <- mad(nFt_abs_dev, feature_med)
  mt_mad <- mad(mt_abs_dev, pct_mt_med)
  rb_mad <- mad(rb_abs_dev, pct_rb_med)
  
  return(c(nCt_mad, nFt_mad, mt_mad, rb_mad))
}

mad_filter <- function(obj) {
  count_med <- med_calc(obj)[1]
  feature_med <- med_calc(obj)[2]
  pct_mt_med <- med_calc(obj)[3]
  pct_rb_med <- med_calc(obj)[4]
  
  nCt_mad <- mad_calc(obj)[1]
  nFt_mad <- mad_calc(obj)[2]
  mt_mad <- mad_calc(obj)[3]
  rb_mad <- mad_calc(obj)[4]
  
  # Filtering Cells if they fall within the range
  obj@meta.data[["nCount_Cutoff"]] <- ifelse(obj@meta.data$nCount_RNA >= (count_med - (3*nCt_mad)) & obj@meta.data$nCount_RNA <= (count_med + (3*nCt_mad)), TRUE, FALSE)
  obj@meta.data[["nFeature_Cutoff"]] <- ifelse(obj@meta.data$nFeature_RNA >= (feature_med - (3*nFt_mad)) & obj@meta.data$nFeature_RNA <= (feature_med + (3*nFt_mad)), TRUE, FALSE)
  obj@meta.data[["percent.mt_Cutoff"]] <- ifelse(obj@meta.data$percent.mt >= (pct_mt_med - (3*mt_mad)) & obj@meta.data$percent.mt <= (pct_mt_med + (3*mt_mad)), TRUE, FALSE)
  obj@meta.data[["percent.rb_Cutoff"]] <- ifelse(obj@meta.data$percent.rb >= (pct_rb_med - (3*rb_mad)) & obj@meta.data$percent.rb <= (pct_rb_med + (3*rb_mad)), TRUE, FALSE)
  
  # Subsetting the cells
  obj_mad <- subset(obj, nCount_Cutoff == TRUE & nFeature_Cutoff == TRUE & percent.mt_Cutoff == TRUE & percent.rb_Cutoff == TRUE)
  
  return(c(obj, obj_mad))
}

# amb_RNA <- function(obj, filt_mtx, raw_mtx) {
  # Create a SoupChannel object
  soup_channel <- SoupChannel(raw_mtx, filt_mtx)
  
  # Quick Cluster
  obj <- SCTransform(obj, verbose = F)
  obj <- RunPCA(obj, verbose = F)
  obj <- RunUMAP(obj, dims = 1:30, verbose = F)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = F)
  obj <- FindClusters(obj, verbose = T)
  
  obj_meta    <- obj@meta.data
  obj_umap    <- obj@reductions$umap@cell.embeddings
  soup_channel  <- setClusters(soup_channel, setNames(obj_meta$seurat_clusters, rownames(obj_meta)))
  soup_channel  <- setDR(soup_channel, obj_umap)
  
  # Calculate contamination fraction
  soup_channel_filt  <- autoEstCont(soup_channel)
  
  # Adjust the count matrix
  obj_adj  <- adjustCounts(soup_channel_filt, roundToInt = T)
  
  # Create Seurat object from adjusted matrix
  obj_amb <- CreateSeuratObject(counts = obj_adj)
  
  return(c(obj_amb))
}

# # Original Objects 
HCD_og <- mad_filter(HCD)[1]
HCP_og <- mad_filter(HCP)[1]
IWD_og <- mad_filter(IWD)[1]
IWP_og <- mad_filter(IWP)[1]
P0D_og <- mad_filter(P0D)[1]
P0P_og <- mad_filter(P0P)[1]
# # Modified MAD Objects
HCD_mad <- mad_filter(HCD)[2]
HCP_mad <- mad_filter(HCP)[2]
IWD_mad <- mad_filter(IWD)[2]
IWP_mad <- mad_filter(IWP)[2]
P0D_mad <- mad_filter(P0D)[2]
P0P_mad <- mad_filter(P0P)[2]



# # Ambient RNA Algorithm
# HCD_amb <- amb_RNA(HCD, HCD.data, HCD_raw)
# HCD_amb[[1]][["percent.mt"]] <- PercentageFeatureSet(HCD_amb[[1]], pattern = "^MT-")
# HCD_amb[[1]][["percent.rb"]] <- PercentageFeatureSet(HCD_amb[[1]], pattern = "^RP[SL]")
# saveRDS(HCD_amb[[1]], file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/HCD_SoupX.rds")
# 
# HCP_amb <- amb_RNA(HCP, HCP.data, HCP_raw)
# HCP_amb[[1]][["percent.mt"]] <- PercentageFeatureSet(HCP_amb[[1]], pattern = "^MT-")
# HCP_amb[[1]][["percent.rb"]] <- PercentageFeatureSet(HCP_amb[[1]], pattern = "^RP[SL]")
# saveRDS(HCP_amb[[1]], file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/HCP_SoupX.rds")
# 
# IWD_amb <- amb_RNA(IWD, IWD.data, IWD_raw)
# IWD_amb[[1]][["percent.mt"]] <- PercentageFeatureSet(IWD_amb[[1]], pattern = "^MT-")
# IWD_amb[[1]][["percent.rb"]] <- PercentageFeatureSet(IWD_amb[[1]], pattern = "^RP[SL]")
# saveRDS(IWD_amb[[1]], file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/IWD_SoupX.rds")
# 
# IWP_amb <- amb_RNA(IWP, IWP.data, IWP_raw)
# IWP_amb[[1]][["percent.mt"]] <- PercentageFeatureSet(IWP_amb[[1]], pattern = "^MT-")
# IWP_amb[[1]][["percent.rb"]] <- PercentageFeatureSet(IWP_amb[[1]], pattern = "^RP[SL]")
# saveRDS(IWP_amb[[1]], file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/IWP_SoupX.rds")
# 
# P0D_amb <- amb_RNA(P0D, P0D.data, P0D_raw)
# P0D_amb[[1]][["percent.mt"]] <- PercentageFeatureSet(P0D_amb[[1]], pattern = "^MT-")
# P0D_amb[[1]][["percent.rb"]] <- PercentageFeatureSet(P0D_amb[[1]], pattern = "^RP[SL]")
# saveRDS(P0D_amb[[1]], file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/P0D_SoupX.rds")
# 
# P0P_amb <- amb_RNA(P0P, P0P.data, P0P_raw)
# P0P_amb[[1]][["percent.mt"]] <- PercentageFeatureSet(P0P_amb[[1]], pattern = "^MT-")
# P0P_amb[[1]][["percent.rb"]] <- PercentageFeatureSet(P0P_amb[[1]], pattern = "^RP[SL]")
# saveRDS(P0P_amb[[1]], file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/P0P_SoupX.rds")


VlnPlot(HCD_og[[1]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), pt.size = 0.00001, ncol = 4)
VlnPlot(HCD_mad[[1]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), pt.size = 0.00001, ncol = 4)

# Plots
# HCD_amb[[1]] <- SetIdent(HCD_amb[[1]], value = "HCE048-D")
# VlnPlot(HCD_amb[[1]], features = c("nCount_RNA"), pt.size = 0.01)

# RidgePlot(HCD_og[[1]], features = c("nCount_RNA")) + geom_vline(xintercept = count_med) + 
#   geom_vline(xintercept = count_med - (3*nCt_mad)) + geom_vline(xintercept = count_med + (3*nCt_mad))
# RidgePlot(HCD_og[[1]], features = c("nFeature_RNA")) + geom_vline(xintercept = feature_med) + 
#   geom_vline(xintercept = feature_med - (2*nFt_mad)) + geom_vline(xintercept = feature_med + (2*nFt_mad))
 
# # Subsetting the cells
# HCD_soup <- subset(HCD, nCount_lower == TRUE)
# 
# VlnPlot(HCD_soup, features = c("nCount_RNA"), pt.size = 0.01)


cutoff <- function(obj) {
  feature_med4 <- med_calc(obj)[2]
  obj@meta.data[["nFeature_Cutoff"]] <- ifelse(obj@meta.data$nFeature_RNA >=(.45*feature_med4) &
                                                 obj@meta.data$nFeature_RNA <= (3*feature_med4), TRUE, FALSE)
  mt_med4 <- med_calc(obj)[4]
  obj@meta.data[["mt_Cutoff"]] <- ifelse(obj@meta.data$percent.mt >=(.2*mt_med4) &
                                           obj@meta.data$percent.mt <= (1.5*mt_med4), TRUE, FALSE)
  obj1 <- subset(obj, nFeature_Cutoff == TRUE | mt_Cutoff == TRUE)
  return(obj1)
}
# Subsetted
HCD1 <- cutoff(HCD)
HCP1 <- cutoff(HCP)
IWD1 <- cutoff(IWD)
IWP1 <- cutoff(IWP)
P0D1 <- cutoff(P0D)
P0P1 <- cutoff(P0P)

# Merging filtered samples
eso.big <- merge(HCD1, y = c(HCP1, IWD1, IWP1, P0D1, P0P1), add.cell.ids = c("HCD", "HCP", "IWD","IWP", "P0D", "P0P"), project = "Sample")

VlnPlot(eso.big, features = c("nFeature_RNA"), pt.size = 0, group.by = "orig.ident")
VlnPlot(eso.big, features = c("percent.mt"), pt.size = 0, group.by = "orig.ident")

# Merging unfiltered samples
og_com <- merge(HCD, y = c(HCP, IWD, IWP, P0D, P0P), add.cell.ids = c("HCD", "HCP", "IWD","IWP", "P0D", "P0P"), project = "Orig.Sample")
VlnPlot(og_com, features = c("nFeature_RNA", "percent.mt"), pt.size = 0, group.by = "orig.ident")

VlnPlot(eso.big, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0, ncol=3) + ggtitle("PostQC")

# Filtered Sample Plot
ggdraw() +
  draw_plot(VlnPlot(eso.big, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0, ncol = 3)) +
  draw_label("PostQC", x = 0.5, y = 0.01, fontface = "bold", hjust = 0.5)

# Unfiltered Sample Plot
ggdraw() +
  draw_plot(VlnPlot(og_com, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0, ncol = 3)) +
  draw_label("PreQC", x = 0.5, y = 0.01, fontface = "bold", hjust = 0.5)

