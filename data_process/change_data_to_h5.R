library(Seurat)
library(SeuratDisk)
library(scPred)

SaveH5Seurat(scPred::pbmc_1, filename = "pbmc_1.h5Seurat")
SaveH5Seurat(scPred::pbmc_2, filename = "pbmc_2.h5Seurat")

Convert("pbmc_1.h5Seurat", dest = "h5ad", overwrite = TRUE)
Convert("pbmc_2.h5Seurat", dest = "h5ad", overwrite = TRUE)
