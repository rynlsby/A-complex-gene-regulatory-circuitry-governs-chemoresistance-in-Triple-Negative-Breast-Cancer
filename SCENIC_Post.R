# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
#library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

library(SCopeLoomR)
loom <- open_loom("PRE_ALL.loom")
# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
cellClusters <- get_cell_annotation(loom)
close_loom(loom)



### To check whether it was read properly:
length(regulons);  head(names(regulons))
regulonAUC
length(regulonAucThresholds)
plot(embeddings$`umap on full expression matrix`)


cellClusters <- cellClusters[,-(2),drop=FALSE]
head(cellClusters)


selectedResolution <- "Chemoresitant" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,"Response"]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))


# plot:
options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

viewTable(topRegulators, options = list(pageLength = 10))
write.csv(topRegulators, "topRegulators.csv")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters)
rss<- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellClusters$Cell.type)
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

options(repr.plot.width=15, repr.plot.height=15) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = c("Chemoresistant", "Chemosensitive"),n = 50) # cluster ID
plotRSS(rss, )
# List of embeddings available:
cat(names(embeddings), sep="\n")

# Overview of these embeddings (see below for details)
regulonsToPlot <- "SP1(+)"
options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2, ceiling(length(names(embeddings))/2)))
for (selectedEmbedding in names(embeddings))
  AUCell::AUCell_plotTSNE(embeddings[[selectedEmbedding]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, 
                          sub=selectedEmbedding)


selectedEmbedding <- embeddings[["umap on full expression matrix"]]

tfsToPlot <- c("TFAP2C", "TFAP2A", "SP1") 
regulonsToPlot <-c("SP1(+)")  
  
options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(3,3))
# Plot expression:
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log[tfsToPlot,], plots=c("Expression"), cex = .5)
# Plot regulon activity:
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, thresholds = 0.2)
dev.off()
options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(selectedEmbedding, .5)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=6, drawlabels=FALSE)  
  
regulonsToPlot <- "SP1(+)"
options(repr.plot.width=10, repr.plot.height=4) # To set the figure size in Jupyter
par(mfrow=c(1,3))
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, 
                        regulonAUC[regulonsToPlot,], thresholds = regulonAucThresholds[regulonsToPlot],
                        plots=c("AUC", "histogram", "binary"), cex = .5, thresholds = 0.15)  
aucellApp <- AUCell_createViewerApp(auc=regulonAUC,
                                    thresholds=regulonAucThresholds,
                                    tSNE=selectedEmbedding, 
                                    exprMat=exprMat_log)
savedSelections <- shiny::runApp(aucellApp)  
