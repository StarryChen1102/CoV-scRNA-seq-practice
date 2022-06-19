library(Seurat)
library(dplyr)
library(sctransform)
library(patchwork)

# read data

ileum.data <- read.table("ileum/ileum_raw_UMIcounts.txt", header = TRUE, sep = "\t")
row.names(ileum.data) <- ileum.data$GENE
ileum.data <- subset(ileum.data, select = -c(GENE))
ileum.data <- t(ileum.data)

cell.info <- read.csv("ileum/ileum_cell_info.txt", sep = "\t")
rownames(cell.info) <- cell.info$UniqueCell_ID
cell.info <- subset(cell.info, select = -c(UniqueCell_ID))
rownames(cell.info) <- gsub("-", ".", rownames(cell.info))

ileum.selected <- ileum.data[cell.info$Sample_ID == "Ileum-1" | cell.info$Sample_ID == "Ileum-2",]
#ileum.selected <- ileum.selected[seq(1, nrow(ileum.selected.1), 3),]

ileum.selected <- t(ileum.selected)

# create Seurat object

ileum <- CreateSeuratObject(counts = ileum.selected, project = "ileum cells", min.cells = 3, min.features = 200)

# quality control

VlnPlot(ileum, features = c("nFeature_RNA", "nCount_RNA"))
plot1 <- FeatureScatter(ileum, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

ileum <- subset(ileum, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)

plot1 <- FeatureScatter(ileum, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

# normalization

ileum <- SCTransform(ileum)

# PCA

ileum <- RunPCA(ileum, features = rownames(ileum.selected))
ElbowPlot(ileum)
VizDimLoadings(ileum, dims = 1:5, reduction = "pca")
DimHeatmap(ileum, dims = 1:5, cells = 100, balanced = TRUE)

# UMAP and t-SNE

ileum <- RunUMAP(ileum, dims = 1:8)
ileum <- RunTSNE(ileum)

plot1 <- DimPlot(ileum, reduction = "pca")
plot2 <- DimPlot(ileum, reduction = "umap")
plot3 <- DimPlot(ileum, reduction = "tsne")
plot1 + plot2 + plot3

# clustering

ileum <- FindNeighbors(ileum, dims = 1:8)
ileum <- FindClusters(ileum)

plot1 <- DimPlot(ileum, reduction = "pca")
plot2 <- DimPlot(ileum, reduction = "umap")
plot3 <- DimPlot(ileum, reduction = "tsne")
plot1 + plot2 + plot3

# find markers

ileum.markers <- FindAllMarkers(ileum, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "scale.data")
ileum.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_diff)

ileum.features <- c("STMN1","HMGB2","DMBT1","ASCL2","RGMB","ALDOB","CA7","SPIB","SPINK4","TFF3","CHGA","NEUROD1","SEPP1","APOC3")
ileum.markers[ileum.features,]

DotPlot(ileum, features = ileum.features, group.by = "seurat_clusters") + RotatedAxis()

VlnPlot(ileum, features = c("ACE2","ANPEP","ENPEP","DPP4"), slot = "scale.data")

# cell type annotation and visualization

new.cluster.ids <- c("PROs","ECs","ECs","PROs","ECs","SCs","ECs","PROs","ECs","TAs","SCs","PROs4","Gs","Gs","SCs","PCs","EECs")
names(new.cluster.ids) <- levels(ileum)
ileum <- RenameIdents(ileum, new.cluster.ids)

plot2 <- DimPlot(ileum, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot2

VlnPlot(ileum, features = "ACE2", slot = "scale.data")
VlnPlot(ileum, features = "ANPEP", slot = "scale.data")
VlnPlot(ileum, features = "ENPEP", slot = "scale.data")
VlnPlot(ileum, features = "DPP4", slot = "scale.data")

DotPlot(ileum, features = ileum.features) + RotatedAxis()

# count the number of cells

table(Idents(ileum))
