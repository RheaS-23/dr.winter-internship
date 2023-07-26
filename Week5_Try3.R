library(Seurat)
library(dplyr)
library(caret)
library(glmnet)
library(nnet)
library(tidymodels)

# Load data
training <- readRDS("/Users/ankithachalla/Desktop/RStudio/SingleCellProject/SingleCells/training_merged_obj.rds")

# Subset Epithelial Cells
training_epi <- subset(training, annotation.cluster == "Epithelial")

# Getting sample and condition data as a data frame
condition_data <- data.frame("condition" = training_epi@meta.data$condition, "sample" = training_epi@meta.data$orig.ident)

# Remove repetitive rows to create keys matching samples to conditions
condition_data <- condition_data %>% distinct()

# Get pseudobulk values per sample
training_pseudobulk <- AverageExpression(training_epi, group.by = "orig.ident", slot = "data", assay = "RNA")

# Transpose -- genes are now columns
training_pseudobulk <- training_pseudobulk$RNA %>% t() 

# Convert into data frame 
training_pseudobulk <- training_pseudobulk %>% as.data.frame()

# # Transpose data frame so genes are columns
# training_pseudobulk <- t(training_pseudobulk)

# Add sample names as a column of data for binding condition data
training_pseudobulk["sample"] <- rownames(training_pseudobulk)

# Combine average expression data with condition data
# full_join will splice together two data sets based on matching values in the "sample" column
training_pseudobulk <- full_join(training_pseudobulk, condition_data, by = "sample")

# Get rid of the "sample" column so it doesn't mess with your regression model
training_pseudobulk <- training_pseudobulk %>% select(-sample)

# Get row indexes to assign to training data set
sample_indexes <- createDataPartition(y = training_pseudobulk$condition, p = 0.7)

# Select training data set by rows 
training_set <- training_pseudobulk[sample_indexes$Resample1,]
test_set <- training_pseudobulk[-sample_indexes$Resample1,]

# Setting cross validation
control <- trainControl(method="cv", number=10)

# MAYBE THIS POTENTIALLY WORKS
model <- train(condition ~ ., data = training_set, method = "multinom", trControl = control)


