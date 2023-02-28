library(Seurat)
P01_Extinct_Clusters <- readRDS("P01_Extinct_Clusters.rds")
meta<-read.delim("P1_Meta.txt", header = T, row.names = 1)
P01_Extinct_Clusters <- AddMetaData(
  object = P01_Extinct_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)

P02_Extinct_Clusters <- readRDS("P02_Extinct_Clusters.rds")
meta<-read.delim("P2_Meta.txt", header = T, row.names = 1)
P02_Extinct_Clusters <- AddMetaData(
  object = P02_Extinct_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)

P06_Extinct_Clusters <- readRDS("P06_Extinct_Clusters.rds")
meta<-read.delim("P6_Meta.txt", header = T, row.names = 1)
P06_Extinct_Clusters <- AddMetaData(
  object = P06_Extinct_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)

P09_Extinct_Clusters <- readRDS("P09_Extinct_Clusters.rds")
meta<-read.delim("P9_Meta.txt", header = T, row.names = 1)
P09_Extinct_Clusters <- AddMetaData(
  object = P09_Extinct_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)

P10_Persist_Clusters <- readRDS("P10_Extinct_Clusters.rds")
meta<-read.delim("P10_Meta.txt", header = T, row.names = 1)
P10_Persist_Clusters <- AddMetaData(
  object = P10_Persist_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)

P11_Persist_Clusters <- readRDS("P11_Extinct_Clusters.rds")
meta<-read.delim("P11_Meta.txt", header = T, row.names = 1)
P11_Persist_Clusters <- AddMetaData(
  object = P11_Persist_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)

P14_Persist_Clusters <- readRDS("P14_Extinct_Clusters.rds")
meta<-read.delim("P14_Meta.txt", header = T, row.names = 1)
P14_Persist_Clusters <- AddMetaData(
  object = P14_Persist_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)

P15_Persist_Clusters <- readRDS("P15_Extinct_Clusters.rds")
meta<-read.delim("P15_Meta.txt", header = T, row.names = 1)
P15_Persist_Clusters <- AddMetaData(
  object = P15_Persist_Clusters,
  metadata = meta,
  col.name = 'Treatment'
)


P01_Extinct_Clusters$Patients <- 'P01.Chemosensitive'
P02_Extinct_Clusters$Patients <- 'P02.Chemosensitive'
P06_Extinct_Clusters$Patients <- 'P06.Chemosensitive'
P09_Extinct_Clusters$Patients <- 'P09.Chemosensitive'
P10_Persist_Clusters$Patients <- 'P10.Chemoresistant'
P11_Persist_Clusters$Patients <- 'P11.Chemoresistant'
P14_Persist_Clusters$Patients <- 'P14.Chemoresistant'
P15_Persist_Clusters$Patients <- 'P15.Chemoresistant'



P01_Extinct_Clusters$Response <- 'Chemosensitive'
P02_Extinct_Clusters$Response <- 'Chemosensitive'
P06_Extinct_Clusters$Response <- 'Chemosensitive'
P09_Extinct_Clusters$Response <- 'Chemosensitive'
P10_Persist_Clusters$Response <- 'Chemoresistant'
P11_Persist_Clusters$Response <- 'Chemoresistant'
P14_Persist_Clusters$Response <- 'Chemoresistant'
P15_Persist_Clusters$Response <- 'Chemoresistant'





combined <- merge(
  x = P01_Extinct_Clusters,
  y = list(P02_Extinct_Clusters,P06_Extinct_Clusters,P09_Extinct_Clusters, 
           P10_Persist_Clusters, P11_Persist_Clusters, P14_Persist_Clusters, 
           P15_Persist_Clusters),
  add.cell.ids = c("B1", "B2" ,"B3", "B4"),
  merge.data = TRUE,
  project = "BC_MET"
)


combined <- combined <- subset(combined, subset = Treatment == Pre)


combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
combined <- subset(combined, subset = nFeature_RNA > 200 & percent.mt < 20)
combined

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(combined, split.by = "Response")



# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})



# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
anchors.combined <- IntegrateData(anchorset = anchors)



# specify that we will perform downstream analysis on the corrected data note that the
DefaultAssay(anchors.combined) <- "integrated"




# Run the standard workflow for visualization and clustering
anchors.combined <- ScaleData(anchors.combined, verbose = FALSE)
anchors.combined <- RunPCA(anchors.combined, npcs = 30, verbose = FALSE)
anchors.combined <- RunUMAP(anchors.combined, reduction = "pca", dims = 1:7)
anchors.combined <- FindNeighbors(anchors.combined, reduction = "pca", dims = 1:7)
anchors.combined <- FindClusters(anchors.combined, resolution = 0.4)
# Visualization
anchors.combined@meta.data

DimPlot(anchors.combined, reduction = "umap")
DimPlot(anchors.combined, reduction = "umap", group.by = "Response")








