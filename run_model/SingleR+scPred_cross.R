library(scPred)
library(Seurat)
library(SingleR)
library(scater)
library(SingleCellExperiment)
library(dplyr)
library(caret)
library(Matrix)

reference <- scPred::pbmc_1
query <- scPred::pbmc_2

# scPred

# Step 1: reference preprocessing
reference <- reference %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# Step 2: create scPred space + model training
reference <- getFeatureSpace(reference, "cell_type")

start_train <- Sys.time()
reference <- trainModel(reference)
end_train <- Sys.time()
training_time_scpred <- as.numeric(difftime(end_train, start_train, units = "secs"))

# Step 3: query normalization and prediction
query <- NormalizeData(query)

start_test <- Sys.time()
query <- scPredict(query, reference)
end_test <- Sys.time()
testing_time_scpred <- as.numeric(difftime(end_test, start_test, units = "secs"))

# getting the predicted labels
pred_labels_scpred <- query$scpred_prediction
true_labels <- query$cell_type

# save the results
write.csv(true_labels, "scPred_True_Labels_cross.csv", row.names = FALSE)
write.csv(pred_labels_scpred, "scPred_Pred_Labels_cross.csv", row.names = FALSE)
write.csv(data.frame(Training_Time = training_time_scpred), "scPred_Training_Times_cross.csv", row.names = FALSE)
write.csv(data.frame(Testing_Time = testing_time_scpred), "scPred_Testing_Times_cross.csv", row.names = FALSE)



# SingleR

# change to SingleCellExperiment
sce_ref <- as.SingleCellExperiment(reference)
sce_query <- as.SingleCellExperiment(query)

# logNorm
sce_ref <- logNormCounts(sce_ref)
sce_query <- logNormCounts(sce_query)
ref_labels <- sce_train$celltype

# SingleR prediction
start_test <- Sys.time()
pred <- SingleR(test = sce_query, ref = sce_ref, labels = sce_ref$cell_type)
end_test <- Sys.time()
testing_time_singleR <- as.numeric(difftime(end_test, start_test, units = "secs"))

# getting the predicted labels
pred_labels_singleR <- pred$labels
true_labels <- query$cell_type

# save the results
write.csv(true_labels, "SingleR_True_Labels_cross.csv", row.names = FALSE)
write.csv(pred_labels_singleR, "SingleR_Pred_Labels_cross.csv", row.names = FALSE)
write.csv(data.frame(Testing_Time = testing_time_singleR), "SingleR_Testing_Times_cross.csv", row.names = FALSE)
