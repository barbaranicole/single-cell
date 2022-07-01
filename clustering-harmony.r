library(Seurat)
#packageVersion("Seurat")
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(matrixStats)
library(patchwork)
library(clustree)
library(pheatmap)
library(scran)
library(harmony)
library(cluster)
library("RColorBrewer")

cs <- readRDS("merge_norm2.rds")
print(cs)

cs <- RunPCA(object = cs,
             assay = "RNA",
             npcs = 50, #no mas de 50 
             rev.pca = FALSE,
             weight.by.var = TRUE,
             verbose = TRUE,
             ndims.print = 1:5,
             nfeatures.print = 50,
             reduction.key = "PC_",
             seed.use = 42,
             approx = TRUE)

DimPlot(cs, reduction = "pca", group.by = "sample")

options(repr.plot.height = 8, repr.plot.width = 10)

VizDimLoadings(cs, dims = 1:2, reduction = "pca")  #genes que explican cada unno de los componentes
# son genes que contribuyen mas a la varianza en ese componente 
#tomar estos genes y compararlos con un database de vias o referencia de funciones - deberian estar relacionados 
DimHeatmap(object = cs,
           dims = 1:2,
           nfeatures = 30,
           cells = NULL,
           reduction = "pca",
           slot = "scale.data",
           assays = "RNA",
           combine = TRUE)

options(repr.plot.height = 10, repr.plot.width = 8)
VizDimLoadings(cs, dims = 15:20, reduction = "pca")

DimHeatmap(object = cs,
           dims = 10:15,
           nfeatures = 30,
           cells = NULL,
           reduction = "pca",
           slot = "scale.data",
           assays = "RNA",
           combine = TRUE)

ElbowPlot(cs, ndims = 50)

cs <- RunHarmony(cs, "sample")

cs <- FindNeighbors(cs, reduction = "harmony", dims = 1:20, k.param = 20)

cs <- FindClusters(cs, graph.name = "RNA_snn", 
                         resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), algorithm = 1)

options(repr.plot.height = 12, repr.plot.width = 10)
clustree(cs@meta.data, prefix = "RNA_snn_res.")

ResolutionList <- grep("_snn_res", colnames(cs@meta.data), value = TRUE)
ResolutionList <- ResolutionList[!ResolutionList %in% 'RNA_snn_res.0']

#can not be use with resolution 0
options(repr.plot.height = 12, repr.plot.width = 9)
dist.pca <- dist(cs@reductions$harmony@cell.embeddings[,c(1:20)])

for (r in ResolutionList){
    
    pca_clus <- cs@meta.data[r][, 1]    
    sil <- silhouette(as.numeric(pca_clus), dist.pca)
    
    #set colors
    n <- length(unique(cs@meta.data[r][, 1]))
    
    #pdf(paste0(r, "_sil.pdf"), width=7, height=7)
    plot(sil, border = NA, main=(paste('Silhoutte of ', r)), col=1:n)

    #dev.off()
}

cs <- RunUMAP(cs, reduction = "harmony", dims = 1:20, n.neighbors = 15, min.dist = 0.1, seed.use = 42)

options(repr.plot.height = 6, repr.plot.width = 8)
pca_umap <- DimPlot(cs, reduction = "umap", group.by = "sample")
pca_umap

options(repr.plot.height = 7, repr.plot.width = 20)
pca_umap <- DimPlot(cs, reduction = "umap", group.by = c("RNA_snn_res.0.3", "RNA_snn_res.0.4", "RNA_snn_res.0.5"), 
                    label = TRUE, label.size = 5)
pca_umap

library(scales)

hue_pal()(11)

cs <- FindClusters(cs, resolution = 0.3, algorithm = 1)
cells_to_highlight <- CellsByIdentities(object = cs, idents = c(0))
pca_umap0 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#F8766D")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(1))
pca_umap1 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#DB8E00")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(2))
pca_umap2 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                    cells.highlight = cells_to_highlight, cols.highlight = "#AEA200")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(3))
pca_umap3 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#64B200")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(4))
pca_umap4 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00BD5C")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(5))
pca_umap5 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00C1A7")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(6))
pca_umap6 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00BADE")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(7))
pca_umap7 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00A6FF")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(8))
pca_umap8 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#B385FF")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(9))
pca_umap9 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#EF67EB")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(10))
pca_umap10 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#FF63B6")

options(repr.plot.height = 12, repr.plot.width = 20)

pca_umap0 + pca_umap1 + pca_umap2 + pca_umap3 +pca_umap4 + pca_umap5 + 
pca_umap6 + pca_umap7 + pca_umap8 + pca_umap9 + pca_umap10

hue_pal()(13)

cs <- FindClusters(cs, resolution = 0.4, algorithm = 1)
cells_to_highlight <- CellsByIdentities(object = cs, idents = c(0))
pca_umap0 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#F8766D")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(1))
pca_umap1 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#E18A00")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(2))
pca_umap2 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                    cells.highlight = cells_to_highlight, cols.highlight = "#BE9C00")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(3))
pca_umap3 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#8CAB00")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(4))
pca_umap4 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#24B700")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(5))
pca_umap5 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00BE70")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(6))
pca_umap6 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00C1AB")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(7))
pca_umap7 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00BBDA")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(8))
pca_umap8 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00ACFC")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(9))
pca_umap9 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#8B93FF")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(10))
pca_umap10 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#D575FE")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(11))
pca_umap11 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#F962DD")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(12))
pca_umap12 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#FF65AC")

options(repr.plot.height = 15, repr.plot.width = 20)

pca_umap0 + pca_umap1 + pca_umap2 + pca_umap3 +pca_umap4 + pca_umap5 + 
pca_umap6 + pca_umap7 + pca_umap8 + pca_umap9 + pca_umap10 + pca_umap11 + pca_umap12

hue_pal()(14)

cs <- FindClusters(cs, resolution = 0.5, algorithm = 1)
cells_to_highlight <- CellsByIdentities(object = cs, idents = c(0))
pca_umap0 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#F8766D")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(1))
pca_umap1 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#E38900")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(2))
pca_umap2 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                    cells.highlight = cells_to_highlight, cols.highlight = "#C49A00")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(3))
pca_umap3 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#99A800")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(4))
pca_umap4 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#53B400")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(5))
pca_umap5 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00BC56")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(6))
pca_umap6 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00C094")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(7))
pca_umap7 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00BFC4")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(8))
pca_umap8 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#00B6EB")


cells_to_highlight <- CellsByIdentities(object = cs, idents = c(9))
pca_umap9 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#06A4FF")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(10))
pca_umap10 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#A58AFF")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(11))
pca_umap11 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#DF70F8")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(12))
pca_umap12 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#FB61D7")

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(13))
pca_umap13 <- DimPlot(cs, reduction = "umap", label = TRUE, label.size = 5, 
                     cells.highlight = cells_to_highlight, cols.highlight = "#FF66A8")

options(repr.plot.height = 15, repr.plot.width = 20)

pca_umap0 + pca_umap1 + pca_umap2 + pca_umap3 +pca_umap4 + pca_umap5 + pca_umap6 + 
pca_umap7 + pca_umap8 + pca_umap9 + pca_umap10 + pca_umap11 + pca_umap12 + pca_umap13

cells_to_highlight <- CellsByIdentities(object = cs, idents = c(7))
options(repr.plot.height = 7, repr.plot.width = 20)

pca_umap <- DimPlot(cs, reduction = "umap",  group.by = c("RNA_snn_res.0.3", "RNA_snn_res.0.4", "RNA_snn_res.0.5"),
                    label = TRUE, label.size = 5, cells.highlight = cells_to_highlight)
pca_umap

cs <- FindClusters(cs, resolution = 0.3, algorithm = 1)
allcluster.markers03 <- FindAllMarkers(cs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cs <- FindClusters(cs, resolution = 0.4, algorithm = 1)
allcluster.markers04 <- FindAllMarkers(cs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cs <- FindClusters(cs, resolution = 0.5, algorithm = 1)
allcluster.markers05 <- FindAllMarkers(cs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

allcluster.markers05 <- FindAllMarkers(cs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_number <- 5
top2_markers <- allcluster.markers03 %>% group_by(cluster) %>% top_n(n = top_number, wt = avg_log2FC)

markers.to.plot <- unique(top2_markers$gene)

options(repr.plot.height = 10, repr.plot.width = 20)
DotPlot(cs, features = markers.to.plot, dot.scale = 8) + RotatedAxis()

top_number <- 5
top2_markers <- allcluster.markers04 %>% group_by(cluster) %>% top_n(n = top_number, wt = avg_log2FC)

markers.to.plot <- unique(top2_markers$gene)

options(repr.plot.height = 10, repr.plot.width = 20)
DotPlot(cs, features = markers.to.plot, dot.scale = 8) + RotatedAxis()

top_number <- 5
top2_markers <- allcluster.markers05 %>% group_by(cluster) %>% top_n(n = top_number, wt = avg_log2FC)

markers.to.plot <- unique(top2_markers$gene)

options(repr.plot.height = 10, repr.plot.width = 20)
DotPlot(cs, features = markers.to.plot, dot.scale = 8) + RotatedAxis()

options(repr.plot.height = 10, repr.plot.width = 18)
top10 <- allcluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(cs, features = top10$gene) + NoLegend()

jobname <- "merge2_clust_harmony"
saveRDS(cs, paste0(jobname, ".rds"))

#subset(data1.markers, subset = cluster == 0 & p_val < 0.05 & avg_log2FC > 0.5)

write.table(allcluster.markers, "C:/Users/bgbar/Documents/Work/lladser/merge2_clusters_03.txt", sep="\t")

cs@meta.data

cs <- readRDS("merge_clust_harmony.rds")
print(cs)

c_n <- unique(cs@meta.data$seurat_clusters)

for (n in c_n){
    print(n)
    print(length(WhichCells(cs, idents = c(n))))
    }

test <- GetAssayData(cs, slot = "counts")[, WhichCells(cs, ident = c(12))]

test

library(data.table)
data_to_write_out <- as.data.frame(as.matrix(test))
fwrite(x = data_to_write_out, file = "test_cluster12.csv")

merge_count <- data.frame(ncount = colSums(test))
merge_count %>% arrange(-ncount) %>% mutate(rank = row_number()) -> merge_count

gg_depth <- ggplot(merge_count, aes(ncount)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 1000) +
    geom_vline(xintercept = ceiling(quantile(merge_count$ncount, 0.5)),
               linetype = "dashed", color = "black") +
    ggplot2::annotate("text",
                      label = paste("mean ~", 
                                    ceiling(mean(merge_count$ncount))), 
                      x = mean(merge_count$ncount),
                      y = max(hist(merge_count$ncount, 
                                   breaks = seq(from = 0, 
                                                to = max(merge_count$ncount),
                                                by = max(merge_count$ncount)/1000),
                                   plot = FALSE)$count),
                      hjust = 0) +
    theme_classic() +
    
    ggtitle("Barcode count depth")


gg50<- ggplot(merge_count, aes(ncount)) +

    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 500) +
    scale_x_continuous(limits = c(-.5, ceiling(quantile(merge_count$ncount, 0.5))), n.breaks = 10) +
    #xlim(c(-.5, ceiling(quantile(pac1$ncount, 0.5)))) +

    theme_classic() +
    ggtitle("Count depth - Quantile 50")


gg_log <- ggplot(merge_count, aes(ncount)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
    geom_vline(xintercept = ceiling(quantile(merge_count$ncount, 0.5)),
               linetype = "dashed", color = "black") +
    ggplot2::annotate("text",
                      label = paste("mean ~", 
                                    ceiling(mean(merge_count$ncount))), 
                      x = mean(merge_count$ncount),
                      y = max(hist(log10(merge_count$ncount), 
                                   breaks = seq(from = 0, 
                                                to = max(log10(merge_count$ncount)),
                                                by = max(log10(merge_count$ncount))/100),
                                   plot = FALSE)$count),
                      hjust = 0) +
    theme_classic() +
    scale_x_log10() +
    ggtitle("Log10 Barcode count depth")

options(repr.plot.height = 12, repr.plot.width = 16)
plot_grid(plotlist = list(gg_depth, gg50, gg_log), nrow = 3)
