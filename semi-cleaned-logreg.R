# Test sets/Saving Files
# everything is basically the same in everything except for changing what you subset
library(Seurat)
library(dplyr)
library(caret)

# Load data
training <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/training_merged_obj.rds")

load("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/sc_eso_testSamples.RData")

# Normalize individual samples testing data
test_merged_obj <- NormalizeData(object = test_merged_obj, assay = "RNA")

# Find differentially expressed genes within the epithelial cell subset
markers <- FindMarkers(training, ident.1 = "Epithelial", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

markers %>% 
  top_n(n = 500, wt = avg_log2FC) -> top500

## Epithelial - 70% Training Data -------------------------------------------
# Subset Epithelial Cells
training_epi <- subset(training, annotation.cluster == "Epithelial")

# Getting sample and condition data as a data frame
condition_data <- data.frame("condition"= training_epi@meta.data$condition, "sample" = training_epi@meta.data$orig.ident)

# Remove repetitive rows to create keys matching samples to conditions
condition_data <- condition_data %>% distinct()

# Get pseudobulk values per sample
training_pseudobulk <- AverageExpression(training_epi, group.by = "orig.ident", slot = "data", assay = "RNA")

# Transpose -- genes are now columns
training_pseudobulk <- training_pseudobulk$RNA %>% t() 

# Convert into data frame 
training_pseudobulk <- training_pseudobulk %>% as.data.frame()

# Get aggregated expression of all genes across samples
gene_sums <- training_pseudobulk %>% colSums()

# Select genes with no expression across samples
drop <- which(gene_sums == 0)

# Remove genes with no expression in any sample from data set (lower memory requirement)
training_pseudobulk <- training_pseudobulk %>% dplyr::select(-all_of(drop))

# Add sample names as a column for binding condition data
training_pseudobulk["sample"] <- rownames(training_pseudobulk)

# Combine average expression data with condition data (into one data frame)
# full_join will splice together two data sets based on matching values in the "sample" column
training_pseudobulk <- full_join(training_pseudobulk, condition_data, by = "sample")

# Get rid of the "sample" column so it doesn't mess with your regression model
training_pseudobulk <- training_pseudobulk %>% dplyr::select(-sample)

# Get row indexes to assign to training data set/test data set
sample_indexes <- createDataPartition(y = training_pseudobulk$condition, p = 0.7)

# Select training data and test data set by row indexes
set.seed(47)
training_set <- training_pseudobulk[sample_indexes$Resample1,]
training_set$condition <- as.factor(training_set$condition)
test_set <- training_pseudobulk[-sample_indexes$Resample1,]
test_set$condition <- as.factor(test_set$condition)

saveRDS(training_set, file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_train_set.rds")
saveRDS(test_set, file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_set.rds")

## All Cells - 30% Testing Data -------------------------------------
# Getting sample and condition data as a data frame
condition_data <- data.frame("condition" = training@meta.data$condition, "sample" = training@meta.data$orig.ident)

# Remove repetitive rows to create keys matching samples to conditions
condition_data <- condition_data %>% distinct()

# Get pseudobulk values per sample
training_pseudobulk <- AverageExpression(training, group.by = "orig.ident", slot = "data", assay = "RNA")

# Transpose -- genes are now columns
training_pseudobulk <- training_pseudobulk$RNA %>% t() 

# Convert into data frame 
training_pseudobulk <- training_pseudobulk %>% as.data.frame()

# Select only top500 differential genes
training_pseudobulk <- training_pseudobulk %>% dplyr::select(all_of(rownames(top500)))

# Add sample names as a column for binding condition data
training_pseudobulk["sample"] <- rownames(training_pseudobulk)

# Combine average expression data with condition data (into one data frame)
# full_join will splice together two data sets based on matching values in the "sample" column
training_pseudobulk <- full_join(training_pseudobulk, condition_data, by = "sample")

# Get rid of the "sample" column so it doesn't mess with your regression model
training_pseudobulk <- training_pseudobulk %>% dplyr::select(-sample)

# Get row indexes to assign to training data set/test data set
sample_indexes <- createDataPartition(y = training_pseudobulk$condition, p = 0.7)

# Select training data and test data set by row indexes
set.seed(47)
test_set <- training_pseudobulk[-sample_indexes$Resample1,]
test_set$condition <- as.factor(test_set$condition)

saveRDS(test_set, file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_all1.rds")

## All Cells - Independent Samples -------------------------------------------

# Getting sample and condition data as a data frame
condition_data <- data.frame("condition" = test_merged_obj@meta.data$condition, "sample" = test_merged_obj@meta.data$orig.ident)

# Remove repetitive rows to create keys matching samples to conditions
condition_data <- condition_data %>% distinct()

# Get pseudobulk values per sample
training_pseudobulk <- AverageExpression(test_merged_obj, group.by = "orig.ident", slot = "data", assay = "RNA")

# Transpose -- genes are now columns
training_pseudobulk <- training_pseudobulk$RNA %>% t() 

# Convert into data frame 
training_pseudobulk <- training_pseudobulk %>% as.data.frame()

# Select only top500 differential genes
training_pseudobulk <- training_pseudobulk %>% dplyr::select(all_of(rownames(top500)))

# Add sample names as a column for binding condition data
training_pseudobulk["sample"] <- rownames(training_pseudobulk)

# Combine average expression data with condition data (into one data frame)
# full_join will splice together two data sets based on matching values in the "sample" column
training_pseudobulk <- full_join(training_pseudobulk, condition_data, by = "sample")

# Get rid of the "sample" column so it doesn't mess with your regression model
training_pseudobulk <- training_pseudobulk %>% dplyr::select(-sample)

training_pseudobulk$condition <- as.factor(training_pseudobulk$condition)

saveRDS(training_pseudobulk, file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_all.rds")

#----------------------------
# Select only top500 differential genes
training_set <- training_set %>% dplyr::select(all_of(rownames(top500)))

saveRDS(training_set, file = "/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_training.rds")






# All the training stuff
#--------------------------------------------------------------
library(Seurat)
library(dplyr)
library(caret)
library(glmnet)
library(nnet)
library(tidymodels)
library(SingleCellExperiment)
library(umap)
library(BiocSingular)
library(MLmetrics)
library(pROC)
library(data.table)
library(cvms)
library(datawizard)
## Set Up --------------------------------
# Load all data
training <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/training_merged_obj.rds")

# Subset Epithelial Cells
training_epi <- subset(training, annotation.cluster == "Epithelial")

# Load training and testing set - Building Model
training_set <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_train_set.rds")
test_set_epi <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_set.rds")

# Load Epithelial Cells
test_set_epi2 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_set1.rds")
test_epi <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_epi.rds")

# Load Myeloid Cells
test_my <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_my.rds")

# Load Lymphocytes 
test_lym <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_lym.rds")

# Load Pericytes 
test_per <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_per.rds")

# Load Endothelial 
test_end <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_end.rds")

# Load Mast Cells 
test_mc <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_mc.rds")

# Load Fibroblasts
test_fib <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_fib.rds")

# Load Distal Cells
test_dist <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_dist.rds")

# Load Proximal Cells
test_prox <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_prox.rds")

# Load All Cells
test_all <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_all.rds")

# Load HCE43-P
test_HCP <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_HCE43P.rds")

# Load HCE43-D
test_HCD <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_HCD.rds")

# Load IW0077-D
test_IWD <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_IWD.rds")

# Load IW0077-P
test_IWP <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_IWP.rds")

# Load P0541-D
test_P0D <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_P0D.rds")

# Load P0541-P
test_P0P <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week6_test_P0P.rds")



# Load Epithelial Cells
test_set_epi2 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_set1.rds")

# Load Immune Cells (try 2)
test_imm1 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_imm1.rds")


# Load Submucosal Glands
test_sub <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_sub.rds")

# Load Submucosal Glands (try2)
test_sub2 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_sub2.rds")

# Load Endothelial
test_end <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_end.rds")

# Load Endothelial (try2)
test_end1 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_end1.rds")

# Load Smooth Mucosal
test_smo <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_smo.rds")

# Load Smooth Mucosal (try2)
test_smo1 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_smo1.rds")

# Load Glandular Epithelial
test_gla <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_gla.rds")

# Load Glandular Epithelial (try2)
test_gla1 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_gla1.rds")

# Load Distal Cells
test_dist <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_dist.rds")

# Load Distal Cells (try2)
test_dist1 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_dist1.rds")

# Load Proximal Cells
test_prox <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_prox.rds")

# Load Proximal Cells (try2)
test_prox1 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_prox1.rds")

# # Load HC
# test_HC <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_HC.rds")

# Load All Cells
test_all <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_all.rds")

# Load All Cells (try1)
test_all1 <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/week5_test_all1.rds")


## DIFFERENTIAL GENES -----------------------------------------------

# Find differentially expressed genes within the epithelial cell subset
markers <- FindMarkers(training, ident.1 = "Epithelial", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

markers %>% 
  top_n(n = 500, wt = avg_log2FC) -> top500

training_diff <- training_set %>% dplyr::select(any_of(rownames(top500)), "condition")
test_diff <- test_set %>% dplyr::select(any_of(rownames(top500)), "condition")


set.seed(123)
log_model_diff <- train(condition~., data = training_diff, method = "multinom",
                        trControl = trainControl(method = "cv", number = 10), 
                        MaxNWts = 7000, tuneGrid = data.frame(decay = .002))

# Calculate Probabilities
pred_prob_diff <- predict(log_model_diff, newdata = test_diff, type = "prob")

# AUC
multiclass.roc(test_diff$condition, pred_prob_diff)

pred_prob_diff <- colnames(pred_prob_diff)[max.col(pred_prob_diff)]
pred_prob_diff <- factor(pred_prob_diff, levels = levels(test_diff$condition))

# Results/Stats
diff_matrix <- caret::confusionMatrix(pred_prob_diff, test_diff$condition)


# Validation w/ 30% epithelial cells
set.seed(123)
log_model_diff <- train(condition~., data = training_diff, method = "multinom",
                        trControl = trainControl(method = "cv", number = 10), 
                        MaxNWts = 7000, tuneGrid = data.frame(decay = .002))

# Calculate Probabilities
pred_prob_diff <- predict(log_model_diff, newdata = test_set_epi2, type = "prob")

# AUC
multiclass.roc(test_set_epi2$condition, pred_prob_diff)

pred_prob_diff <- colnames(pred_prob_diff)[max.col(pred_prob_diff)]
pred_prob_diff <- factor(pred_prob_diff, levels = levels(test_set_epi2$condition))

# Results/Stats
diff_matrix <- caret::confusionMatrix(pred_prob_diff, test_set_epi2$condition)


## Hyperparameter Tuning ---------------------------------------------------
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary,
                           search = "random")
hyp_fit_search <- train(condition~., data = training_diff, 
                        method = "multinom",
                        tuneLength = 30,
                        trControl = fitControl, MaxNWts = 7000)

## By condition
## TESTING by cell type -------------------------------------------------
# --- Epithelial Cells --------------

# Calculate Probabilities
pred_epi <- predict(log_model_diff, newdata = test_epi, type = "prob")

# AUC
multiclass.roc(test_epi$condition, pred_epi)

# Convert to match lengths
pred_epi <- colnames(pred_epi)[max.col(pred_epi)]
pred_epi <- factor(pred_epi, levels = levels(test_epi$condition))

# Results/Stats
epi_matrix <- caret::confusionMatrix(pred_epi, test_epi$condition)

# --- Myeloid Cells --------------

# Calculate Probabilities
pred_my <- predict(log_model_diff, newdata = test_my, type = "prob")

# AUC
multiclass.roc(test_my$condition, pred_my)

# Convert to match lengths
pred_my <- colnames(pred_my)[max.col(pred_my)]
pred_my <- factor(pred_my, levels = unique(test_my$condition))
test_my$condition <- factor(test_my$condition, levels = c("HC", "GERD", "SSc"))

# Results/Stats
my_matrix <- caret::confusionMatrix(pred_my, test_my$condition)

# --- Lymphocytes ---------------

# Calculate Probabilities
pred_lym <- predict(log_model_diff, newdata = test_lym, type = "prob")

# AUC
multiclass.roc(test_lym$condition, pred_lym)

# Convert to match lengths
pred_lym <- colnames(pred_lym)[max.col(pred_lym)]
pred_lym <- factor(pred_lym, levels = unique(test_lym$condition))

# Results/Stats
lym_matrix <- caret::confusionMatrix(pred_lym, test_lym$condition)

# --- Pericytes ---------------

# Make sure test_set has all three levels
test_per$condition <- factor(test_per$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_per <- predict(log_model_diff, newdata = test_per, type = "prob")

# AUC
multiclass.roc(test_per$condition, pred_per)

# Convert to match lengths
pred_per <- colnames(pred_per)[max.col(pred_per)]
pred_per <- factor(pred_per, levels = c("GERD", "HC", "SSc"))


# Results/Stats
per_matrix <- caret::confusionMatrix(pred_per, test_per$condition)


# --- Endothelial ---------------

# Make sure test_set has all three levels
test_end$condition <- factor(test_end$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_end <- predict(log_model_diff, newdata = test_end, type = "prob")

# AUC
multiclass.roc(test_end$condition, pred_end)

# Convert to match lengths
pred_end <- colnames(pred_end)[max.col(pred_end)]
pred_end <- factor(pred_end, levels = c("GERD", "HC", "SSc"))

# Results/Stats
end_matrix <- caret::confusionMatrix(pred_end, test_end$condition)



# --- Mast Cells ---------------

# Make sure test_set has all three levels
test_mc$condition <- factor(test_mc$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_mc <- predict(log_model_diff, newdata = test_mc, type = "prob")

# AUC
multiclass.roc(test_mc$condition, pred_mc)

# Convert to match lengths
pred_mc <- colnames(pred_mc)[max.col(pred_mc)]
pred_mc <- factor(pred_mc, levels = c("GERD", "HC", "SSc"))

# Results/Stats
mc_matrix <- caret::confusionMatrix(pred_mc, test_mc$condition)



# --- Fibroblasts ---------------

# Make sure test_set has all three levels
test_fib$condition <- factor(test_fib$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_fib <- predict(log_model_diff, newdata = test_fib, type = "prob")

# AUC
multiclass.roc(test_fib$condition, pred_fib)

# Convert to match lengths
pred_fib <- colnames(pred_fib)[max.col(pred_fib)]
pred_fib <- factor(pred_fib, levels = c("GERD", "HC", "SSc"))

# Results/Stats
fib_matrix <- caret::confusionMatrix(pred_fib, test_fib$condition)

# --- Distal Cells ---------------

# Make sure test_set has all three levels
test_dist$condition <- factor(test_dist$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_dist <- predict(log_model_diff, newdata = test_dist, type = "prob")

# AUC
multiclass.roc(test_dist$condition, pred_dist)

# Convert to match lengths
pred_dist <- colnames(pred_dist)[max.col(pred_dist)]
pred_dist <- factor(pred_dist, levels = c("GERD", "HC", "SSc"))

# Results/Stats
dist_matrix <- caret::confusionMatrix(pred_dist, test_dist$condition)

# --- Proximal Cells ---------------

# Make sure test_set has all three levels
test_prox$condition <- factor(test_prox$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_prox <- predict(log_model_diff, newdata = test_prox, type = "prob")

# AUC
multiclass.roc(test_prox$condition, pred_prox)

# Convert to match lengths
pred_prox <- colnames(pred_prox)[max.col(pred_prox)]
pred_prox <- factor(pred_prox, levels = c("GERD", "HC", "SSc"))

# Results/Stats
prox_matrix <- caret::confusionMatrix(pred_prox, test_prox$condition)


# --- HC Cells ---------------
test_HC <- copy(test_all)
levels(test_HC$condition)[levels(test_HC$condition)!='HC'] <- 'other'

# Calculate Probabilities
pred_HC <- predict(log_model_diff, newdata = test_HC, type = "prob")

# AUC
multiclass.roc(test_HC$condition, pred_HC)

# Convert to match lengths
pred_HC <- colnames(pred_HC)[max.col(pred_HC)]
pred_HC <- factor(pred_HC, levels = c("other", "HC"))
if (any(is.na(pred_HC))) {
  pred_HC[is.na(pred_HC)] <- "other"
}

# Results/Stats
HC_matrix <- caret::confusionMatrix(pred_HC, test_HC$condition)


# --- GERD Cells ---------------
test_GERD <- copy(test_all)
levels(test_GERD$condition)[levels(test_GERD$condition)!='GERD'] <- 'other'

# Calculate Probabilities
pred_GERD <- predict(log_model_diff, newdata = test_GERD, type = "prob")

# Convert to match lengths
pred_GERD <- colnames(pred_GERD)[max.col(pred_GERD)]
pred_GERD <- factor(pred_GERD, levels = c("other", "GERD"))
if (any(is.na(pred_GERD))) {
  pred_GERD[is.na(pred_GERD)] <- "other"
}

# Results/Stats
GERD_matrix <- caret::confusionMatrix(pred_GERD, test_GERD$condition)

# --- SSc Cells ---------------
test_SSc <- copy(test_all)
levels(test_SSc$condition)[levels(test_SSc$condition)!='SSc'] <- 'other'

# Calculate Probabilities
pred_SSc <- predict(log_model_diff, newdata = test_SSc, type = "prob")

# Convert to match lengths
pred_SSc <- colnames(pred_SSc)[max.col(pred_SSc)]
pred_SSc <- factor(pred_SSc, levels = c("other", "SSc"))
if (any(is.na(pred_SSc))) {
  pred_SSc[is.na(pred_SSc)] <- "other"
}

# Results/Stats
SSc_matrix <- caret::confusionMatrix(pred_SSc, test_SSc$condition)


# --- HCD  ---------------

# Make sure test_set has all three levels
test_HCD$condition <- factor(test_HCD$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_HCD <- predict(log_model_diff, newdata = test_HCD, type = "prob")

# AUC
all_auc <- multiclass.roc(test_HCD$condition, pred_HCD)

# Convert to match lengths
pred_HCD <- colnames(pred_HCD)[max.col(pred_HCD)]
pred_HCD <- factor(pred_HCD, levels = c("GERD", "HC", "SSc"))

# Results/Stats
HCD_matrix <- caret::confusionMatrix(pred_HCD, test_HCD$condition)


# --- All Cells ---------------

# Make sure test_set has all three levels
test_all$condition <- factor(test_all$condition, levels = c("GERD", "HC", "SSc"))

# Calculate Probabilities
pred_all <- predict(log_model_diff, newdata = test_all, type = "prob")

# AUC
all_auc <- multiclass.roc(test_all$condition, pred_all)

# Convert to match lengths
pred_all <- colnames(pred_all)[max.col(pred_all)]
pred_all <- factor(pred_all, levels = c("GERD", "HC", "SSc"))

# Results/Stats
all_matrix <- caret::confusionMatrix(pred_all, test_all$condition)

# Plotting confusion matrix
all_cm <- as.data.frame(all_matrix$table)
all_cm$Prediction <- factor(all_cm$Prediction, levels=levels(all_cm$Prediction))

ggplot(all_cm, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c("GERD","HC","SSc")) +
  scale_y_discrete(labels=c("GERD","HC","SSc")) +
  ggtitle("Individual Samples: All Cells")






# --- Endothelial --------

# Calculate Probabilities
pred_end1 <- predict(log_model_diff, newdata = test_end1, type = "prob")

# AUC
multiclass.roc(test_end1$condition, pred_end1)

# Convert to match lengths
pred_end1 <- colnames(pred_end1)[max.col(pred_end1)]
pred_end1 <- factor(pred_end1, levels = levels(test_end1$condition))

# Results/Stats
end_matrix <- caret::confusionMatrix(pred_end1, test_end1$condition)

test_imm <- test_imm %>% dplyr::select(any_of(rownames(top500)), "condition")
# Calculate Probabilities
pred_imm <- predict(log_model_diff, newdata = test_imm, type = "prob")

# AUC
multiclass.roc(test_diff$condition, pred_imm)

# Convert to match lengths
pred_imm <- colnames(pred_imm)[max.col(pred_imm)]
pred_imm <- factor(pred_imm, levels = levels(test_imm$condition))

# Results/Stats
imm_matrix <- caret::confusionMatrix(pred_imm, test_imm$condition)

# --- Smooth Mucosal --------------

# Calculate Probabilities
pred_smo1 <- predict(log_model_diff, newdata = test_smo1, type = "prob")

# AUC
multiclass.roc(test_smo1$condition, pred_smo1)

# Convert to match lengths
pred_smo1 <- colnames(pred_smo1)[max.col(pred_smo1)]
pred_smo1 <- factor(pred_smo1, levels = levels(test_smo1$condition))

# Results/Stats
smo_matrix <- caret::confusionMatrix(pred_smo1, test_smo1$condition)

test_imm <- test_imm %>% dplyr::select(any_of(rownames(top500)), "condition")
# Calculate Probabilities
pred_imm <- predict(log_model_diff, newdata = test_imm, type = "prob")

# AUC
multiclass.roc(test_diff$condition, pred_imm)

# Convert to match lengths
pred_imm <- colnames(pred_imm)[max.col(pred_imm)]
pred_imm <- factor(pred_imm, levels = levels(test_imm$condition))

# Results/Stats
imm_matrix <- caret::confusionMatrix(pred_imm, test_imm$condition)

# --- Glandular Epithelial ----

# Calculate Probabilities
pred_gla1 <- predict(log_model_diff, newdata = test_gla1, type = "prob")

# AUC
multiclass.roc(test_gla1$condition, pred_gla1)

# Convert to match lengths
pred_gla1 <- colnames(pred_gla1)[max.col(pred_gla1)]
pred_gla1 <- factor(pred_gla1, levels = levels(test_gla1$condition))

# Results/Stats
gla_matrix <- caret::confusionMatrix(pred_gla1, test_gla1$condition)

test_imm <- test_imm %>% dplyr::select(any_of(rownames(top500)), "condition")
# Calculate Probabilities
pred_imm <- predict(log_model_diff, newdata = test_imm, type = "prob")

# AUC
multiclass.roc(test_diff$condition, pred_imm)

# Convert to match lengths
pred_imm <- colnames(pred_imm)[max.col(pred_imm)]
pred_imm <- factor(pred_imm, levels = levels(test_imm$condition))

# Results/Stats
imm_matrix <- caret::confusionMatrix(pred_imm, test_imm$condition)

## TESTING by location ------------------------------------------
# --- Distal ---------
test_runumap_dist <- test_dist[,-26881]
test_umap_dist <- umap(test_runumap_dist, n_neighbors = 3)
test_umap_dist <- as.data.frame(test_umap_dist$layout)

# Rejoining condition column
test_umap_dist <- test_umap_dist %>% cbind(test_dist$condition)
colnames(test_umap_dist)[colnames(test_umap_dist) == "test_dist$condition"] <- "condition"


# Calculate Probabilities
pred_dist_umap <- predict(log_model_umap, newdata = test_umap_dist, type = "prob")

# AUC
multiclass.roc(test_umap_dist$condition, pred_dist_umap)

# Convert to match lengths
pred_dist_umap <- colnames(pred_dist_umap)[max.col(pred_dist_umap)]
pred_dist_umap <- factor(pred_dist_umap, levels = levels(test_umap_dist$condition))

# Results/Stats
dist_matrix <- caret::confusionMatrix(pred_dist_umap, test_umap_dist$condition)

# Distal (try2) -------
# Calculate Probabilities
pred_dist1 <- predict(log_model_diff, newdata = test_dist1, type = "prob")

# AUC
multiclass.roc(test_dist1$condition, pred_dist1)

# Convert to match lengths
pred_dist1 <- colnames(pred_dist1)[max.col(pred_dist1)]
pred_dist1 <- factor(pred_dist1, levels = levels(test_dist1$condition))

# Results/Stats
dist1_matrix <- caret::confusionMatrix(pred_dist1, test_dist1$condition)

# --- Proximal ---------
test_runumap_prox <- test_prox[,-26806]
test_umap_prox <- umap(test_runumap_prox, n_neighbors = 3)
test_umap_prox <- as.data.frame(test_umap_prox$layout)

# Rejoining condition column
test_umap_prox <- test_umap_prox %>% cbind(test_prox$condition)
colnames(test_umap_prox)[colnames(test_umap_prox) == "test_prox$condition"] <- "condition"

# Calculate Probabilities
pred_prox_umap <- predict(log_model_umap, newdata = test_umap_prox, type = "prob")

# AUC
multiclass.roc(test_umap_prox$condition, pred_prox_umap)

# Convert to match lengths
pred_prox_umap <- colnames(pred_prox_umap)[max.col(pred_prox_umap)]
pred_prox_umap <- factor(pred_prox_umap, levels = levels(test_umap_prox$condition))

# Results/Stats
prox_matrix <- caret::confusionMatrix(pred_prox_umap, test_umap_prox$condition)

# Proximal (try2) -------
# Calculate Probabilities
pred_prox1 <- predict(log_model_diff, newdata = test_prox1, type = "prob")

# AUC
multiclass.roc(test_prox1$condition, pred_prox1)

# Convert to match lengths
pred_prox1 <- colnames(pred_prox1)[max.col(pred_prox1)]
pred_prox1 <- factor(pred_prox1, levels = levels(test_prox1$condition))

# Results/Stats
prox1_matrix <- caret::confusionMatrix(pred_prox1, test_prox1$condition)

## TESTING by condition ------------------------------------------

# --- HC  ---------
test_HC <- copy(test_all1)

levels(test_HC$condition)[levels(test_HC$condition)!='HC'] <- 'other'

# Calculate Probabilities
pred_HC <- predict(log_model_diff, newdata = test_HC, type = "prob")

# AUC
multiclass.roc(test_HC$condition, pred_HC)

# Convert to match lengths
pred_HC <- colnames(pred_HC)[max.col(pred_HC)]
pred_HC <- factor(pred_HC, levels = levels(test_HC$condition))
if (any(is.na(pred_HC))) {
  pred_HC[is.na(pred_HC)] <- "other"
}

# Results/Stats
HC_matrix <- caret::confusionMatrix(pred_HC, test_HC$condition)

# --- GERD  ------------ 
test_GERD <- copy(test_all1)

levels(test_GERD$condition)[levels(test_GERD$condition)!='GERD'] <- 'other'

# Calculate Probabilities
pred_GERD <- predict(log_model_diff, newdata = test_GERD, type = "prob")

# AUC
multiclass.roc(test_GERD$condition, pred_GERD)

# Convert to match lengths
pred_GERD <- colnames(pred_GERD)[max.col(pred_GERD)]
pred_GERD <- factor(pred_GERD, levels = levels(test_GERD$condition))
if (any(is.na(pred_GERD))) {
  pred_GERD[is.na(pred_GERD)] <- "other"
}

# Results/Stats
GERD_matrix <- caret::confusionMatrix(pred_GERD, test_GERD$condition)

# --- SSc ------------ 
test_SSc <- copy(test_all1)

levels(test_SSc$condition)[levels(test_SSc$condition)!='SSc'] <- 'other'

# Calculate Probabilities
pred_SSc <- predict(log_model_diff, newdata = test_SSc, type = "prob")

# AUC
multiclass.roc(test_SSc$condition, pred_SSc)

# Convert to match lengths
pred_SSc <- colnames(pred_SSc)[max.col(pred_SSc)]
pred_SSc <- factor(pred_SSc, levels = levels(test_SSc$condition))
if (any(is.na(pred_SSc))) {
  pred_SSc[is.na(pred_SSc)] <- "other"
}

# Results/Stats
SSc_matrix <- caret::confusionMatrix(pred_SSc, test_SSc$condition)

## ALL CELLS --------------- 

# Calculate Probabilities
pred_all <- predict(log_model_diff, newdata = test_all1, type = "prob")

# AUC
multiclass.roc(test_all1$condition, pred_all)

# Convert to match lengths
pred_all <- colnames(pred_all)[max.col(pred_all)]
pred_all <- factor(pred_all, levels = levels(test_all1$condition))


# Results/Stats
all_matrix <- caret::confusionMatrix(pred_all, test_all1$condition)

## ALL CELLS - 30% (try2) --------------- 

# Calculate Probabilities
pred_all1 <- predict(log_model_diff, newdata = test_all1, type = "prob")

# AUC
multiclass.roc(test_all1$condition, pred_all)

# Convert to match lengths
pred_all1 <- colnames(pred_all1)[max.col(pred_all1)]
pred_all1 <- factor(pred_all1, levels = levels(test_all1$condition))

# Results/Stats
all1_matrix <- caret::confusionMatrix(pred_all1, test_all1$condition)


# Plotting confusion matrix
all1_cm <- as.data.frame(all1_matrix$table)
all1_cm$Prediction <- factor(all1_cm$Prediction, levels=levels(all1_cm$Prediction))

ggplot(all1_cm, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c("GERD","HC","SSc")) +
  scale_y_discrete(labels=c("GERD","HC","SSc")) +
  ggtitle("30% Held-out Data: All Cells")


