library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(scPred)
library(SingleCellExperiment)
library(caret)

# load data and create Seurat object
data_dir <- "../filtered_matrices_mex/hg19"
pbmc <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = pbmc, min.cells = 3, min.features = 200)

# add cell labels
labels <- read_tsv(file.path(data_dir, "68k_pbmc_barcodes_annotation.tsv"))
labels <- labels %>%
  select(barcodes, celltype) %>%
  filter(!is.na(celltype)) %>%
  column_to_rownames("barcodes")

seurat_obj <- subset(seurat_obj, cells = intersect(Cells(seurat_obj), rownames(labels)))
seurat_obj$celltype <- labels[Cells(seurat_obj), "celltype"]

# data preprocessing
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

labels <- seurat_obj$celltype

# create train and test sets
set.seed(42)
n <- ncol(seurat_obj)
fold_ids <- sample(rep(1:5, length.out = n))
names(fold_ids) <- colnames(seurat_obj)

true_all <- c()
pred_all <- c()
train_times <- c()
test_times <- c()
fold_ids <- sample(rep(1:5, length.out = n))
names(fold_ids) <- colnames(seurat_obj)

train_cells <- names(fold_ids)[fold_ids == 1]
test_cells  <- names(fold_ids)[fold_ids != 1]

seurat_train <- subset(seurat_obj, cells = train_cells)
seurat_test  <- subset(seurat_obj, cells = test_cells)

# train scPred model
seurat_train <- getFeatureSpace(seurat_train, "celltype")

start_train <- Sys.time()
seurat_train <- trainModel(seurat_train)
end_train <- Sys.time()
train_time <- as.numeric(difftime(end_train, start_train, units = "secs"))

# predict with scPred
seurat_test <- NormalizeData(seurat_test, normalization.method = "LogNormalize", scale.factor = 10000)

start_test <- Sys.time()
seurat_test <- scPredict(seurat_test, seurat_train)
end_test <- Sys.time()
test_time <- as.numeric(difftime(end_test, start_test, units = "secs"))

# save results
true_all <- seurat_test$celltype
pred_all <- seurat_test$scpred_prediction

write.csv(true_all, "scPred_True_Labels.csv", row.names = FALSE)
write.csv(pred_all, "scPred_Pred_Labels.csv", row.names = FALSE)
write.csv(train_time, "scPred_Training_Time.csv", row.names = FALSE)
write.csv(test_time, "scPred_Testing_Time.csv", row.names = FALSE)



# SingleR

library(SingleR)
library(scater)

set.seed(42)
n <- ncol(seurat_obj)
fold_ids <- sample(rep(1:5, length.out = n))
names(fold_ids) <- colnames(seurat_obj)

# create train and test sets
train_cells <- names(fold_ids)[fold_ids == 1]
test_cells  <- names(fold_ids)[fold_ids != 1]
seurat_train <- subset(seurat_obj, cells = train_cells)
seurat_test  <- subset(seurat_obj, cells = test_cells)


# change to SCE object
sce_train <- as.SingleCellExperiment(seurat_train)
sce_test  <- as.SingleCellExperiment(seurat_test)

# Log-normalization
sce_train <- logNormCounts(sce_train)
sce_test  <- logNormCounts(sce_test)

# get reference labels
ref_labels <- sce_train$celltype 

# run SingleR
start_test_time <- Sys.time()
pred <- SingleR(test = sce_test, ref = sce_train, labels = sce_train$celltype)
end_test_time <- Sys.time()
testing_time <- as.numeric(difftime(end_test_time, start_test_time, units = "secs"))

# save results
true_labels <- seurat_test$celltype
pred_labels <- pred$labels

write.csv(true_labels, "SingleR_True_Labels.csv", row.names = FALSE)
write.csv(pred_labels, "SingleR_Pred_Labels.csv", row.names = FALSE)
write.csv(data.frame(Testing_Time = testing_time), "SingleR_Testing_Times.csv", row.names = FALSE)

