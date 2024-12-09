library(Seurat)
library(tidyverse)


nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m)

cts <-  nsclc.sparse.m$`Gene Expression`
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj

View(nsclc.seurat.obj@meta.data)
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)


VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Load the ggplot2 package
library(ggplot2)

# Now try to create the plot again
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)


nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)

nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)


plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)




nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(nsclc.seurat.obj)
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")
