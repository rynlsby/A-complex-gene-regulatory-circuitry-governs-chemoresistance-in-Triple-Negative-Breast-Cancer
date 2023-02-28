library(SingleCellExperiment)
library(SCopeLoomR)

obj.se <- as.SingleCellExperiment(obj)

# DGEM (Digital gene expression matrix)
dgem <- counts(obj.se)
dim(dgem)
class(dgem)
head(colnames(dgem))  #should contain the Cell ID/name
# Known cell information/annotation  
cell.info <- colData(obj.se)
cell.info$nGene <- colSums(dgem>0)
head(cell.info)
# Default embedding (e.g. umap or PCA coordinates)



default.umap <- Embeddings(obj, reduction = "umap")
default.umap.name <- "umap on full expression matrix"
head(default.umap)


file.name <- "PRE_ALL.loom"
loom<-build_loom(
  file.name=file.name,
  dgem=dgem,
  title="Pre ALL",
  default.embedding=default.umap,
  default.embedding.name=default.umap.name
)
close_loom(loom)

loom <- open_loom(file.path = file.name, mode = "r+")
loom<-add_col_attr(
  loom=loom,
  key = "Response",
  value=cell.info$Response,
  as.annotation=T
)

close_loom(loom)

