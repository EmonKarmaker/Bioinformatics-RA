library(Seurat)
library(SeuratDisk)
rds_obj<-readRDS('C:\Users\Emon Karmoker\Documents\ependymal_cells.rds')



# Using forward slashes
rds_obj <- readRDS('C:/Users/Emon Karmoker/Documents/ependymal_cells.rds')

df5_obj <- Read10X_h5(filename = '20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
                   use.names = TRUE,
                   unique.features=TRUE)
df5_obj[1:10,1:10]
seurat_hdf5 <- CreateSeuratObject(counts =df5_obj )



mtx_obj <- ReadMtx(mtx = 'C:/Users/Emon Karmoker/Documents/raw_feature_bc_matrix/matrix.mtx',
        features = 'C:/Users/Emon Karmoker/Documents/raw_feature_bc_matrix/features.tsv',
        cells = 'C:/Users/Emon Karmoker/Documents/raw_feature_bc_matrix/barcode.tsv')

# Specify the directory containing the data files
data_dir <- "C:/Users/Emon Karmoker/Documents/raw_feature_bc_matrix"

# Construct file paths for the barcode, matrix, and features files
barcode_file <- file.path(data_dir, "barcodes.tsv.gz")
matrix_file <- file.path(data_dir, "matrix.mtx.gz")
features_file <- file.path(data_dir, "features.tsv.gz")

# Use the ReadMtx function with specified file paths
mtx_obj <- ReadMtx(mtx = matrix_file, features = features_file, cells = barcode_file)

mtx_obj[1:10,1:10]

seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

Convert("adata_SS2_for_download.h5ad" ,dest="h5seurat", overwrite = TRUE )


seurat_anndata <- LoadH5Seurat("adata_SS2_for_download.")
