library(Seurat)
#packageVersion("Seurat")
library(ggplot2)
library(dplyr)
library(cowplot)

library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)

pac1 <- Read10X(data.dir = "C:/Users/bgbar/Documents/Work/lladser/outs_pac2/filtered_feature_bc_matrix/")
pac1 <- pac1$`Gene Expression`

dim(pac1)

pac1 <- data.frame(ncount = colSums(pac1))
pac1 %>% arrange(-ncount) %>% mutate(rank = row_number()) -> pac1

gg_depth <- ggplot(pac1, aes(ncount)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 1000) +
    geom_vline(xintercept = ceiling(quantile(pac1$ncount, 0.5)),
               linetype = "dashed", color = "black") +
    ggplot2::annotate("text",
                      label = paste("mean ~", 
                                    ceiling(mean(pac1$ncount))), 
                      x = mean(pac1$ncount),
                      y = max(hist(pac1$ncount, 
                                   breaks = seq(from = 0, 
                                                to = max(pac1$ncount),
                                                by = max(pac1$ncount)/1000),
                                   plot = FALSE)$count),
                      hjust = 0) +
    theme_classic() +
    
    ggtitle("Barcode count depth")


gg50<- ggplot(pac1, aes(ncount)) +

    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 500) +
    scale_x_continuous(limits = c(-.5, ceiling(quantile(pac1$ncount, 0.5))), n.breaks = 10) +
    #xlim(c(-.5, ceiling(quantile(pac1$ncount, 0.5)))) +

    theme_classic() +
    ggtitle("Count depth - Quantile 50")


gg_log <- ggplot(pac1, aes(ncount)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
    geom_vline(xintercept = ceiling(quantile(pac1$ncount, 0.5)),
               linetype = "dashed", color = "black") +
    ggplot2::annotate("text",
                      label = paste("mean ~", 
                                    ceiling(mean(pac1$ncount))), 
                      x = mean(pac1$ncount),
                      y = max(hist(log10(pac1$ncount), 
                                   breaks = seq(from = 0, 
                                                to = max(log10(pac1$ncount)),
                                                by = max(log10(pac1$ncount))/100),
                                   plot = FALSE)$count),
                      hjust = 0) +
    theme_classic() +
    scale_x_log10() +
    ggtitle("Log10 Barcode count depth")

options(repr.plot.height = 12, repr.plot.width = 16)
plot_grid(plotlist = list(gg_depth, gg50, gg_log), nrow = 3)

sc <- Read10X(data.dir = "C:/Users/bgbar/Documents/Work/lladser/outs_pac2/filtered_feature_bc_matrix/")
sc <- sc$`Gene Expression`
qc <- data.frame(ncount = colSums(sc))

qc[["nfeat"]] <- colSums(sc != 0)

gg_024 <- ggplot(qc, aes(nfeat)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
    geom_vline(xintercept = ceiling(quantile(qc$nfeat, 0.5)),
               linetype = "dashed", color = "black") +
    ggplot2::annotate("text",
                      label = paste("mean ~", 
                                    ceiling(mean(qc$nfeat))), 
                      x = mean(qc$nfeat),
                      y = max(hist(qc$nfeat, 
                                   breaks = seq(from = 0, 
                                                to = max(qc$nfeat),
                                                by = max(qc$nfeat)/100),
                                   plot = FALSE)$count,
                             na.rm = TRUE),
                      hjust = 0) +
    theme_classic() +
    ggtitle("Gene count depth")

gg50_024 <- ggplot(qc, aes(nfeat)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 500) +
    scale_x_continuous(limits = c(-.5, ceiling(quantile(qc$nfeat, 0.5))), n.breaks = 20) +
    #xlim(c(-.5, ceiling(quantile(qc$nfeat, 0.5)))) +
    theme_classic() +
    ggtitle("Gene depth - Quantile 50")


options(repr.plot.height = 9, repr.plot.width = 16)
plot_grid(plotlist = list(gg_024, gg50_024), nrow = 2)

sc <- Read10X(data.dir = "C:/Users/bgbar/Documents/Work/lladser/outs_pac2/filtered_feature_bc_matrix/")
sc <- sc$`Gene Expression`
qc <- data.frame(ncount = colSums(sc))

qc[["nfeat"]] <- colSums(sc != 0)

gg_e <- ggplot(qc, aes(x = ncount, y = nfeat)) +
    geom_point(stat = "identity", alpha = 0.2, size = 0.5) +
    theme_classic() +
    ggtitle("Genes-detected vs. transcript depth")

options(repr.plot.height = 10, repr.plot.width = 10)
plot(gg_e)

sc <- Read10X(data.dir = "C:/Users/bgbar/Documents/Work/lladser/outs_pac2/filtered_feature_bc_matrix/")
sc <- sc$`Gene Expression`
qc <- data.frame(ncount = colSums(sc))


mtgenes <- rownames(sc)[grep("^MT-", rownames(sc))]
qc$mtcount <- colSums(sc[mtgenes, ])
qc$pct_mtcount <- (qc$mtcount/qc$ncount)*100


gg_mt <- ggplot(qc, aes(x = ncount, y = pct_mtcount)) +
        geom_point(color = "red", alpha = .2, size = 0.5) +
        geom_density_2d(color = "black", alpha = 0.5) +
        theme_classic() +
        scale_y_continuous(n.breaks = 10) +
        ggtitle("% Mitochondrial-counts vs. depth")


den_mt <- ggplot(qc, aes(pct_mtcount)) +
        geom_density(fill = "red", color = "white", alpha = 0.8) +
        geom_vline(xintercept = median(qc$pct_mtcount),
        linetype = "dashed", color = "black") +
        ggplot2::annotate("text",
        label = paste("mean ~", 
        round(mean(qc$pct_mtcount), digits=2), "%"), 
            x = mean(qc$pct_mtcount), 
            y = max(hist(qc$pct_mtcount,
                         breaks = seq(from = 0, 
                                      to = max(qc$pct_mtcount),
                                      by = max(qc$pct_mtcount)/1000),
                         plot = FALSE)$density),
                          hjust = 1,
                          vjust = 0) +
        theme_classic() +
        scale_x_continuous(n.breaks = 10) +
        coord_flip() +
        ggtitle("Mitochondrial-counts density")


rbgenes <- rownames(sc)[grep("^RPS|^RPL", rownames(sc))]
qc$rbcount <- colSums(sc[rbgenes, ])
qc$pct_rbcount <- (qc$rbcount/qc$ncount)*100


gg_rb <- ggplot(qc, aes(x = ncount, y = pct_rbcount)) +
        geom_point(color = "red", alpha = .2, size = 0.5) +
        geom_density_2d(color = "black", alpha = 0.5) +
        theme_classic() +
        ggtitle("% Ribosomal-counts vs. depth")


den_rb <- ggplot(qc, aes(pct_rbcount)) +
        geom_density(fill = "red", color = "white", alpha = 0.8) +
        geom_vline(xintercept = median(qc$pct_rbcount),
        linetype = "dashed", color = "black") +
        ggplot2::annotate("text",
        label = paste("mean ~", 
        round(mean(qc$pct_rbcount), digits=2), "%"), 
            x = mean(qc$pct_rbcount), 
            y = max(hist(qc$pct_rbcount,
                         breaks = seq(from = 0, 
                                      to = max(qc$pct_rbcount),
                                      by = max(qc$pct_rbcount)/1000),
                         plot = FALSE)$density),
                          hjust = 1,
                          vjust = 0) +
        theme_classic() +
        coord_flip() +
        ggtitle("Ribosomal-counts density")

qc$ncount

qc$pct_mtcount

options(repr.plot.height = 8, repr.plot.width = 10)
plot_grid(plotlist = list(gg_mt, den_mt), ncol = 2)

options(repr.plot.height = 12, repr.plot.width = 14)
plot_grid(plotlist = list(gg_mt, den_mt, gg_rb, den_rb), ncol = 2, nrow = 2)

sc <- Read10X(data.dir = "C:/Users/bgbar/Documents/Work/lladser/outs_pac2/filtered_feature_bc_matrix/")
sc <- sc$`Gene Expression`
qc <- data.frame(ncount = colSums(sc))

features <- data.frame('feat_count_depth' = rowSums(sc),
                       'id' = rownames(sc))

librarysize <- sum(features$feat_count_depth)
features$pctlib <- (features$feat_count_depth/librarysize)
libsizecell <- colSums(sc)
ntop_genes <- 50

topfeat <- head(arrange(features, -feat_count_depth), ntop_genes)$id
topexp <- sc[topfeat, ]
pctlibsizecell <- lapply(seq(ncol(topexp)), function(i) (topexp[, i]/libsizecell[i])*100)
names(pctlibsizecell) <- colnames(topexp)
           
top <- do.call(rbind,
                 lapply(names(pctlibsizecell), function(cell) {
                     data.frame('pctlibcellsize'=pctlibsizecell[[cell]],
                                'id'=names(pctlibsizecell[[cell]]),
                                'cell'=cell)
                     })
                 )

top <- merge(top, features, by='id')
    
gtop <- ggplot(top, aes(x = id, y = pctlibcellsize)) +
                geom_jitter(aes(group = id), alpha=.3, color = "grey", shape=".") +
                geom_violin(alpha=.3, fill = "grey") +
                scale_x_discrete(limits=rev(head(arrange(features, -feat_count_depth), ntop_genes)$id)) +
                coord_flip() +
                theme_classic() +
                theme(axis.text=element_text(size=5)) +
                ggtitle(paste0("Top ", ntop_genes, " deepest genes"),
                        subtitle = paste0("They account for ~",
                                          ceiling(sum(features[topfeat,]$pctlib)*100),
                                          "% of total-counts"))

plot(gtop)

#devtools::install_github("constantAmateur/SoupX",ref='devel')

library(SoupX)

sample_dir <- "C:/Users/bgbar/Documents/Work/lladser/outs_pac2/"
soup_sc = load10X(sample_dir)
soup_sc = autoEstCont(soup_sc)
sc = adjustCounts(soup_sc)

mtgenes <- rownames(sc)[grep("^MT-", rownames(sc))]
cqc <- data.frame("ncount" = colSums(sc), "mtcount" = colSums(sc[mtgenes, ]))
cqc$pct_mtcount <- (cqc$mtcount/cqc$ncount)*100

good_barcodes <- rownames(cqc[cqc$pct_mtcount < 20, ])
print(length(good_barcodes))
print(head(good_barcodes))

print(dim(sc))
sc <- sc[, good_barcodes]
print(dim(sc))

project <- "pac_2"
min_genes <- 300
cs <- CreateSeuratObject(counts = sc, project = project, min.cells = 3, min.features = min_genes)
print(head(cs@meta.data))

cs

jobname <- "pac2_filt"
saveRDS(cs, paste0(jobname, ".rds"))

## MERGE ##
samples <- c("outs_pac1", "outs_pac2")

merge_pac <- do.call(cbind, 
             lapply(samples, function(s) {
                 sample_name <-  gsub(pattern = "outs_",replacement = "", s)
                 sample_dir <- paste0("C:/Users/bgbar/Documents/Work/lladser/", s, "/filtered_feature_bc_matrix")
                 counts <- Read10X(data.dir = sample_dir)
                 sc <- counts$`Gene Expression`
                 colnames(sc) <- paste(colnames(sc), sample_name, sep = "-")
                 sc
             }))

merge_pac

merge_count <- data.frame(ncount = colSums(merge_pac))
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


merge_feat <- data.frame(ncount = colSums(merge_pac))

merge_feat[["nfeat"]] <- colSums(merge_pac != 0)

gg_024 <- ggplot(merge_feat, aes(nfeat)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
    geom_vline(xintercept = ceiling(quantile(merge_feat$nfeat, 0.5)),
               linetype = "dashed", color = "black") +
    ggplot2::annotate("text",
                      label = paste("mean ~", 
                                    ceiling(mean(merge_feat$nfeat))), 
                      x = mean(merge_feat$nfeat),
                      y = max(hist(merge_feat$nfeat, 
                                   breaks = seq(from = 0, 
                                                to = max(merge_feat$nfeat),
                                                by = max(merge_feat$nfeat)/100),
                                   plot = FALSE)$count,
                             na.rm = TRUE),
                      hjust = 0) +
    theme_classic() +
    ggtitle("Gene count depth")

gg50_024 <- ggplot(merge_feat, aes(nfeat)) +
    geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 500) +
    scale_x_continuous(limits = c(-.5, ceiling(quantile(merge_feat$nfeat, 0.5))), n.breaks = 20) +
    #xlim(c(-.5, ceiling(quantile(qc$nfeat, 0.5)))) +
    theme_classic() +
    ggtitle("Gene depth - Quantile 50")


options(repr.plot.height = 9, repr.plot.width = 16)
plot_grid(plotlist = list(gg_024, gg50_024), nrow = 2)

qc <- data.frame(ncount = colSums(merge_pac))

mtgenes <- rownames(merge_pac)[grep("^MT-", rownames(merge_pac))]
qc$mtcount <- colSums(merge_pac[mtgenes, ])
qc$pct_mtcount <- (qc$mtcount/qc$ncount)*100

gg_mt <- ggplot(qc, aes(x = ncount, y = pct_mtcount)) +
        geom_point(color = "red", alpha = .2, size = 0.5) +
        geom_density_2d(color = "black", alpha = 0.5) +
        theme_classic() +
        scale_y_continuous(n.breaks = 10) +
        ggtitle("% Mitochondrial-counts vs. depth")

den_mt <- ggplot(qc, aes(pct_mtcount)) +
        geom_density(fill = "red", color = "white", alpha = 0.8) +
        geom_vline(xintercept = median(qc$pct_mtcount),
        linetype = "dashed", color = "black") +
        ggplot2::annotate("text",
        label = paste("mean ~", 
        round(mean(qc$pct_mtcount), digits=2), "%"), 
            x = mean(qc$pct_mtcount), 
            y = max(hist(qc$pct_mtcount,
                         breaks = seq(from = 0, 
                                      to = max(qc$pct_mtcount),
                                      by = max(qc$pct_mtcount)/1000),
                         plot = FALSE)$density),
                          hjust = 1,
                          vjust = 0) +
        theme_classic() +
        scale_x_continuous(n.breaks = 10) +
        coord_flip() +
        ggtitle("Mitochondrial-counts density")


rbgenes <- rownames(merge_pac)[grep("^RPS|^RPL", rownames(merge_pac))]
qc$rbcount <- colSums(merge_pac[rbgenes, ])
qc$pct_rbcount <- (qc$rbcount/qc$ncount)*100

gg_rb <- ggplot(qc, aes(x = ncount, y = pct_rbcount)) +
        geom_point(color = "red", alpha = .2, size = 0.5) +
        geom_density_2d(color = "black", alpha = 0.5) +
        theme_classic() +
        ggtitle("% Ribosomal-counts vs. depth")

den_rb <- ggplot(qc, aes(pct_rbcount)) +
        geom_density(fill = "red", color = "white", alpha = 0.8) +
        geom_vline(xintercept = median(qc$pct_rbcount),
        linetype = "dashed", color = "black") +
        ggplot2::annotate("text",
        label = paste("mean ~", 
        round(mean(qc$pct_rbcount), digits=2), "%"), 
            x = mean(qc$pct_rbcount), 
            y = max(hist(qc$pct_rbcount,
                         breaks = seq(from = 0, 
                                      to = max(qc$pct_rbcount),
                                      by = max(qc$pct_rbcount)/1000),
                         plot = FALSE)$density),
                          hjust = 1,
                          vjust = 0) +
        theme_classic() +
        coord_flip() +
        ggtitle("Ribosomal-counts density")

options(repr.plot.height = 12, repr.plot.width = 14)
plot_grid(plotlist = list(gg_mt, den_mt, gg_rb, den_rb), ncol = 2, nrow = 2)

options(repr.plot.height = 8, repr.plot.width = 10)
plot_grid(plotlist = list(gg_mt, den_mt), ncol = 2)

mtgenes <- rownames(merge_pac)[grep("^MT-", rownames(merge_pac))]
cqc <- data.frame("ncount" = colSums(merge_pac), "mtcount" = colSums(merge_pac[mtgenes, ]))
cqc$pct_mtcount <- (cqc$mtcount/cqc$ncount)*100

good_barcodes <- rownames(cqc[cqc$pct_mtcount < 10, ])
print(length(good_barcodes))
print(head(good_barcodes))

print(dim(merge_pac))
merge_pac <- merge_pac[, good_barcodes]
print(dim(merge_pac))

#remove TRC genes
TRAV <- rownames(merge_pac)[grep("^TRAV", rownames(merge_pac))]
TRAJ <- rownames(merge_pac)[grep("^TRAJ", rownames(merge_pac))]
TRBV <- rownames(merge_pac)[grep("^TRBV", rownames(merge_pac))]
TRBJ <- rownames(merge_pac)[grep("^TRBJ", rownames(merge_pac))]
TRDV <- rownames(merge_pac)[grep("^TRDV", rownames(merge_pac))]
TRDJ <- rownames(merge_pac)[grep("^TRDJ", rownames(merge_pac))]
TRGV <- rownames(merge_pac)[grep("^TRGV", rownames(merge_pac))]
TRGJ <- rownames(merge_pac)[grep("^TRGJ", rownames(merge_pac))]

merge_pac <- merge_pac[-(which(rownames(merge_pac) %in% c(TRAV, TRAJ, TRBV, TRBJ, TRDV, TRDJ, TRGV,TRGJ))),]

merge_pac

min_genes <- 300
cs <- CreateSeuratObject(counts = merge_pac, min.cells = 3, min.features = min_genes)
print(head(cs@meta.data))

#add samples to metadata
meta <- cs@meta.data
meta$barcode <- rownames(meta)
b <- meta$barcode
meta$sample <- unlist(regmatches(b, regexec("pac\\d+", b)))
cs@meta.data <- meta[colnames(cs), ]

print(head(cs@meta.data))

print(cs)

jobname <- "merge_filt2"
saveRDS(cs, paste0(jobname, ".rds"))
