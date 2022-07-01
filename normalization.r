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

cs <- readRDS("merge_filt2.rds")
print(cs)

raw_mat <- cs@assays$RNA@counts 
print(dim(raw_mat))
sample_barcodes <- sample(colnames(raw_mat), 100)
print(head(sample_barcodes))
sub_raw_mat <- data.matrix(raw_mat[, sample_barcodes])
srm <- data.frame(melt(sub_raw_mat))
print(head(srm))

#observar el numero de transcritos detectadaos para una celula para ver las diferencias de como fueron secuenciadas 100 celulas aleatorias 
gg_raw <- ggplot(srm, aes(x = Var2, y = value)) +
    geom_boxplot(aes(color = Var2)) +
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()


options(repr.plot.height = 12, repr.plot.width = 10)
plot(gg_raw)

#normalizamos tomando cada umi de cada gen y dividido por el total de transcritos por cada celula y 
#se multiplica por un escalar 
print(dim(cs@assays$RNA@data))
cs <- NormalizeData(object = cs,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000,
                    verbose = TRUE)
print(dim(cs@assays$RNA@data))

#se ve el efecto de la normalizacion dentro de las 100 celulas 
norm_mat <- cs@assays$RNA@data
print(dim(norm_mat))
sub_norm_mat <- data.matrix(norm_mat[, sample_barcodes])
snm <- data.frame(melt(sub_norm_mat))
print(head(snm))

gg_norm <- ggplot(snm, aes(x = Var2, y = value)) +
    geom_boxplot(aes(color = Var2)) +
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#las celulas se encuentran en los mismo intervalos 
plot(gg_norm)

#dependencia de la varianza con la media - se calcula ambos de los valores ya normailzados y se plotea la dependencia
qc <- data.frame('var' = rowVars(data.matrix(norm_mat)),
                 'mean' = rowMeans(norm_mat),
                 'id' = rownames(norm_mat))

print(head(qc))

gg_mv <- ggplot(qc, aes(x = mean, y = var)) +
        geom_point(size = 0.1) +
        geom_smooth() +
        theme_classic() 

#genes poco secuenciados con varianza pequena 
plot(gg_mv)
#cada punto es un gen y la var y media es de los transcritos en todas las celulas del dataset 

#se encuentra la funcion de la curva y se ve la distancia que se encuentran desde la funcion de la dependencia de la varianza con la media 
#es una varianza informativa, cuales son los genes que con una detecccion similar se alejan mas de la varianza esperada
cs <- FindVariableFeatures(cs, 
                           selection.method = "vst", 
                           nfeatures = 2000)
#la varianza del gen solo depende de la profundidad de la secuenciacion

ntop_genes <- 30

top <- head(VariableFeatures(cs), ntop_genes)

gg_var <- VariableFeaturePlot(cs)

gg_var_lab <- LabelPoints(plot = gg_var, 
                          points = top, 
                          repel = TRUE)

plot(gg_var_lab)

gene_var <- modelGeneVar(norm_mat)
fdr_thr <- 0.1

gene_var <- gene_var %>%
    data.frame %>%
    mutate(hvg = ifelse(FDR < fdr_thr, TRUE, FALSE)) %>%
    arrange(FDR)

print(head(gene_var))
print(table(gene_var$hvg))
ntop <- table(gene_var$hvg)[2]
print(ntop)

#se asigna un p value, para cada nivel de deteccion dependiedo de cuanto alejado este de la funcion se puede asignar un p value. Se asume que la distribucion es normal 
top_hvg <- rownames(head(gene_var, ntop))

scs <- cs
scs@assays$RNA@meta.features$vst.variable <- rownames(scs@assays$RNA@data) %in% top_hvg

gg_var <- VariableFeaturePlot(scs)

gg_var_lab_scran <- LabelPoints(plot = gg_var,
                          points = top_hvg, 
                          repel = TRUE)

plot(gg_var_lab_scran)

gg_tvar_scran <- ggplot(gene_var, aes(x = mean, y = total, color = hvg)) +
    geom_point() +
    scale_color_manual(values=c("black", "red")) +
    theme_classic()

gg_bvar_scran <- ggplot(gene_var, aes(x = mean, y = bio, color = hvg)) +
    geom_point() +
    scale_color_manual(values=c("black", "red")) +
    theme_classic()

plot_grid(plotlist = list(gg_tvar_scran, gg_bvar_scran), nrow = 1)

#centrar la distribucion de los genes en 0 
cs <- ScaleData(object = cs,
                features = rownames(cs),
                vars.to.regress = NULL,
                model.use = "linear",
                do.scale = TRUE,
                do.center = TRUE,
                scale.max = 10,
                verbose = TRUE)

sample_genes <- sample(rownames(norm_mat), 50)
print(head(sample_genes))
sub_norm_gmat <- data.matrix(norm_mat[sample_genes, ])
sngm <- data.frame(melt(sub_norm_gmat))
print(head(sngm))

gg_gnorm <- ggplot(sngm, aes(x = Var1, y = value)) +
    geom_boxplot(aes(color = Var1)) +
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

#cada punto es la expresion de el gen en cada barcode 
plot(gg_gnorm)

scale_mat <- cs@assays$RNA@scale.data
sub_scale_gmat <- data.matrix(scale_mat[sample_genes, ])
ssgm <- data.frame(melt(sub_scale_gmat))
print(head(ssgm))

gg_gscale <- ggplot(ssgm, aes(x = Var1, y = value)) +
    geom_boxplot(aes(color = Var1)) +
    theme_classic() +
    theme(legend.position = "none") +
    coord_flip()

plot(gg_gscale)

print(dim(scale_mat))

jobname <- "merge_norm2"
saveRDS(cs, paste0(jobname, ".rds"))

