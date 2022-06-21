library(Seurat)
library(dplyr)
library(sctransform)
library(patchwork)

kn1.data <- Read10X(data.dir = "kidney1")
kn2.data <- Read10X(data.dir = "kidney2")
kn3.data <- Read10X(data.dir = "kidney3")

kn1 <- CreateSeuratObject(counts = kn1.data, project = "kidney1", min.cells = 3, min.features = 200)
kn2 <- CreateSeuratObject(counts = kn2.data, project = "kidney2", min.cells = 3, min.features = 200)
kn3 <- CreateSeuratObject(counts = kn3.data, project = "kidney3", min.cells = 3, min.features = 200)

kn1[["percent.mt"]] <- PercentageFeatureSet(kn1, pattern = "^MT-")
kn2[["percent.mt"]] <- PercentageFeatureSet(kn2, pattern = "^MT-")
kn3[["percent.mt"]] <- PercentageFeatureSet(kn3, pattern = "^MT-")

VlnPlot(kn1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
kn1 <- subset(kn1, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & percent.mt < 25)

VlnPlot(kn2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
kn2 <- subset(kn2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & percent.mt < 25)

VlnPlot(kn3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
kn3 <- subset(kn3, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & percent.mt < 25)

plot1 <- FeatureScatter(kn1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(kn2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(kn3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 + plot3

kn.list <- list(kn1, kn2, kn3)
kn.list <- lapply(X = kn.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = kn.list, nfeatures = 3000)
kn.list <- PrepSCTIntegration(object.list = kn.list, anchor.features = features)

kn.anchors <- FindIntegrationAnchors(object.list = kn.list, normalization.method = "SCT", anchor.features = features)
kn.combined.sct <- IntegrateData(anchorset = kn.anchors, normalization.method = "SCT")

kn.combined.sct <- RunPCA(kn.combined.sct, verbose = FALSE)
ElbowPlot(kn.combined.sct)
VizDimLoadings(kn.combined.sct, dims = 1:2, reduction = "pca")
DimHeatmap(kn.combined.sct, dims = 1:15, cells = 100, balanced = TRUE)

kn.combined.sct <- RunUMAP(kn.combined.sct, reduction = "pca", dims = 1:10)
kn.combined.sct <- RunTSNE(kn.combined.sct)

kn.combined.sct <- FindNeighbors(kn.combined.sct, dims = 1:10)
kn.combined.sct <- FindClusters(kn.combined.sct)
kn.markers <- FindAllMarkers(kn.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "scale.data")
kn.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_diff)

kn.features = c("SLC13A3","GPX3","SLC22A7","DCXR","LYZ","CD14","CD3D","NKG7","DEFB1","CD79A","CD79B","ATP6V1G3","KRT8","KRT18")
kn.markers[kn.features,]
kn.markers[kn.markers$cluster == 13,]

plot1 <- DimPlot(kn.combined.sct, reduction = "pca")
plot2 <- DimPlot(kn.combined.sct, reduction = "umap")
plot3 <- DimPlot(kn.combined.sct, reduction = "tsne")
plot1 + plot2 + plot3

p1 <- DimPlot(kn.combined.sct, reduction = "pca", group.by = "orig.ident")
p2 <- DimPlot(kn.combined.sct, reduction = "umap", group.by = "orig.ident")
p3 <- DimPlot(kn.combined.sct, reduction = "tsne", group.by = "orig.ident")
p1 + p2 + p3

DotPlot(kn.combined.sct, features = kn.features, group.by = "seurat_clusters") + RotatedAxis()
FeaturePlot(kn.combined.sct,features = kn.features, reduction = "umap", slot = "scale.data")

new.cluster.ids <- c("0","1","Distal Ts","3","4","5","6","7","Proximal STs", "Proximal CTs","10","Proximal Ts","Monocytes","NK&T cells","Collecting DIs","15","B cells","Glomerular PEs")
names(new.cluster.ids) <- levels(kn.combined.sct)
kn.combined.sct <- RenameIdents(kn.combined.sct, new.cluster.ids)

plot2 <- DimPlot(kn.combined.sct, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot2

VlnPlot(kn.combined.sct, features = "ANPEP", slot = "scale.data")
VlnPlot(kn.combined.sct, features = "ENPEP", slot = "scale.data")
VlnPlot(kn.combined.sct, features = "DPP4", slot = "scale.data")

kn.combined.sct1 <- subset(kn.combined.sct, idents = c("Distal Ts","Proximal STs","Proximal CTs","Proximal Ts","Monocytes","NK&T cells","Collecting DIs","B cells","Glomerular PEs"))
VlnPlot(kn.combined.sct1, features = "ANPEP", slot = "scale.data")
VlnPlot(kn.combined.sct1, features = "ENPEP", slot = "scale.data")
VlnPlot(kn.combined.sct1, features = "DPP4", slot = "scale.data")
plot2 <- DimPlot(kn.combined.sct1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot2
DotPlot(kn.combined.sct1, features = kn.features) + RotatedAxis()

table(Idents(kn.combined.sct))
