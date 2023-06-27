rm(list = ls())
options(stringsAsFactors = FALSE)

### loading in packages needed ----------------------------------------------------------------------------
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(psych)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)
library(reshape2)
library(ggsignif)
library(cellcall)

### functions ---------------------------------------------------------------------------------------------
DoubletEst <- function (object, pN = 0.25, ndims = 30) {
  for (pkg in c("Seurat", "DoubletFinder")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste(pkg, " package needed for this function to work. Please install it.", 
                 sep = ""), call. = FALSE)
    }
  }
  library(Seurat)
  library(DoubletFinder)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object, features = rownames(object))
  object <- RunPCA(object)
  object <- FindNeighbors(object, dims = 1:ndims)
  object <- FindClusters(object)
  bcmvn <- find.pK(summarizeSweep(paramSweep_v3(object, PCs = 1:30), 
                                  GT = FALSE))
  homotypic.prop <- modelHomotypic(object@meta.data$seurat_clusters)
  nExp_poi <- round((0.076 * ncol(object)/10000) * length(colnames(object)))
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  object <- doubletFinder_v3(object, PCs = 1:10, pN = pN, 
                             pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
  doublet <- colnames(object)[object[[paste0("DF.classifications_", 
                                             as.character(pN), "_", as.character(mpK), "_", as.character(nExp_poi))]] == "Doublet"]
  return(doublet)
}

detect_noise_genes <- function(obs){
  cycle_genes <- intersect(c("UBE2C","HMGB2", "HMGN2", "TUBA1B", "MKI67",
                             "CCNB1", "TUBB", "TOP2A", "TUBB4B"), rownames(obs))
  cycle_score <-  as.numeric(colMeans(as.matrix(obs@assays$RNA@data)[cycle_genes,]))
  
  library(WGCNA)
  gene_data <- data.frame(cycling = cycle_score, 
                          t(as.matrix(obs@assays$RNA@data)),check.names = F)
  cor_data <- cor(gene_data, method = 'pearson')[ ,c('cycling'), drop = F]
  genes <- rownames(cor_data)[cor_data[,'cycling'] > 0.3]
  
  return(genes)
} 

trans2ROE <- function(pro.data) {
  for (i in 1:5) {
    pro.data.2 <- pro.data
    total <- sum(pro.data[(2*i-1):(2*i), 'total_num'])
    for (j in 2:(ncol(pro.data)-1)) {
      sub_num <- sum(pro.data[(2*i-1):(2*i), j])
      pro.data.2[2*i-1, j] <- pro.data.2[2*i-1, 'total_num']*sub_num/total
      pro.data.2[2*i, j] <- pro.data.2[2*i, 'total_num']*sub_num/total
      
      pro.data[2*i-1, j] <- pro.data[2*i-1, j]/pro.data.2[2*i-1, j]
      pro.data[2*i, j] <- pro.data[2*i, j]/pro.data.2[2*i, j]
    }
    rm(total)
  }
  
  pro.data[c('cc4_0', 'cc5_0', 'cc6_0', 'cc7_0', 'cc8_0'), 'Pathology'] <- 'pre-treatment'
  pro.data[c('cc4_1', 'cc5_1', 'cc6_1', 'cc7_1', 'cc8_1'), 'Pathology'] <- 'post-treatment'
  pro.data$Pathology <- factor(pro.data$Pathology, levels = c('pre-treatment', 'post-treatment'))
  pro.data$Patient <- as.character(pro.data$Patient)
  pro.data$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
  
  return(pro.data)
}

### reading post-treatment samples and discarding doublets -----------------------------------------------
post_mat_1 <- Read10X(data.dir = '/data3/analysis/cc/cc4_1/data/filtered_feature_bc_matrix/')
post_mat_2 <- Read10X(data.dir = '/data3/analysis/cc/cc5_1/data/filtered_feature_bc_matrix/')
post_mat_3 <- Read10X(data.dir = '/data3/analysis/cc/cc6_1/data/filtered_feature_bc_matrix/')
post_mat_4 <- Read10X(data.dir = '/data3/analysis/cc/cc7_1/data/filtered_feature_bc_matrix/')
post_mat_5 <- Read10X(data.dir = '/data3/analysis/cc/cc8_1/data/filtered_feature_bc_matrix/')

cc4_1 <- CreateSeuratObject(post_mat_1, min.features = 200)
cc5_1 <- CreateSeuratObject(post_mat_2, min.features = 200)
cc6_1 <- CreateSeuratObject(post_mat_3, min.features = 200)
cc7_1 <- CreateSeuratObject(post_mat_4, min.features = 200)
cc8_1 <- CreateSeuratObject(post_mat_5, min.features = 200)

cc4_1@meta.data$Patient <- 'cc4_1'
cc5_1@meta.data$Patient <- 'cc5_1'
cc6_1@meta.data$Patient <- 'cc6_1'
cc7_1@meta.data$Patient <- 'cc7_1'
cc8_1@meta.data$Patient <- 'cc8_1'

cc4_1_doublet <- DoubletEst(cc4_1)
cc5_1_doublet <- DoubletEst(cc5_1)
cc6_1_doublet <- DoubletEst(cc6_1)
cc7_1_doublet <- DoubletEst(cc7_1)
cc8_1_doublet <- DoubletEst(cc8_1)

cc4_1 <- cc4_1[, !colnames(cc4_1) %in% cc4_1_doublet]
cc5_1 <- cc5_1[, !colnames(cc5_1) %in% cc5_1_doublet]
cc6_1 <- cc6_1[, !colnames(cc6_1) %in% cc6_1_doublet]
cc7_1 <- cc7_1[, !colnames(cc7_1) %in% cc7_1_doublet]
cc8_1 <- cc8_1[, !colnames(cc8_1) %in% cc8_1_doublet]

post.list <- c(cc4_1, cc5_1, cc6_1, cc7_1, cc8_1)
for (i in 1:length(post.list)) {
  post.list[[i]]@meta.data$mt_ratio <- PercentageFeatureSet(post.list[[i]], pattern = '^MT-')
  post.list[[i]] <- post.list[[i]][, post.list[[i]]@meta.data$mt_ratio < 20]
}

### colors setting -------------------------------------------------------------------------------------
color.main <- c("#1B9E77", "#8B2323", "#8B008B", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#D95F02", "#fb832d")
names(color.main) <- levels(Idents(cc))
color.sample <- c(brewer.pal(8, 'Set2'), '#c72e29', '#016392', '#be9c2e', 'purple', 'navy')
names(color.sample) <- levels(cc@meta.data$Patient)
color.pathol <- c('#377EB8', 'firebrick3', 'forestgreen')
names(color.pathol) <- levels(as.factor(cc@meta.data$Pathology))
color.t <- brewer.pal(12, 'Paired')
color.t[11] <- '#666666' 
names(color.t) <- levels(Idents(t.cell))
color.epi <- brewer.pal(5, 'Set1')
names(color.epi) <- levels(Idents(epi))
color.mye <- brewer.pal(7, 'Set1')
names(color.mye) <- levels(Idents(mye.cell))
color.fib <- brewer.pal(6, 'Paired')
names(color.fib) <- levels(Idents(fib))
color.ec <- brewer.pal(5, 'Set1')
names(color.ec) <- levels(Idents(endo))
color.patient.2 <- c( '#c72e29', '#016392', '#be9c2e', 'purple', 'navy')
names(color.patient.2) <- c('cc4', 'cc5', 'cc6', 'cc7', 'cc8')
color.cycle <- c('#458B00', '#EE3B3B', '#8B008B')
names(color.cycle) <- c('G1', 'S', 'G2M')

### loading data of normal and pre-treatment patients -----------------------------------------------------------------
load('/data3/analysis/cc/pre.RData')
colnames(obs.integrated@meta.data)
DimPlot(obs.integrated, group.by = 'patient')
pre_tr <- obs.integrated[, obs.integrated@meta.data$patient %in% c('N1', 'N2', 'N3', 'cc4_0', 'cc5_0', 'cc6_0', 'cc7_0', 'cc8_0')]
pre_tr@meta.data$Patient <- pre_tr@meta.data$patient
pre_tr@meta.data <- pre_tr@meta.data[, c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'mt_ratio', 'Patient')]
pre.list <- SplitObject(pre_tr, split.by = 'Patient')

### integrating and processing data -----------------------------------------------------------------------------------
object.list <- c(post.list, pre.list)
genes.inter <- rownames(cc4_1)
for (i in 2:length(object.list)) {
  genes.inter <- intersect(genes.inter, rownames(object.list[[i]]))
}
for (i in 1:length(object.list)) {
  object.list[[i]] <- object.list[[i]][genes.inter, ]
  object.list[[i]] <- NormalizeData(object.list[[i]])
  object.list[[i]] <- FindVariableFeatures(object.list[[i]]
}

anchors <- FindIntegrationAnchors(object.list)
cc <- IntegrateData(anchors, dims = 1:30)

DefaultAssay(cc) <- 'RNA'
colnames(cc@meta.data)
cc@meta.data <- cc@meta.data[, c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Patient', 'mt_ratio')]
cc <- subset(cc, nCount_RNA < 100000 & nFeature_RNA < 7500)

DefaultAssay(cc) <- 'integrated'
noise.genes <- detect_noise_genes(cc)
cc <- ScaleData(cc, features = setdiff(rownames(cc), noise.genes))
cc <- RunPCA(cc, features = setdiff(VariableFeatures(cc), noise.genes))
ElbowPlot(cc, ndims = 50)
cc <- FindNeighbors(cc, dims = 1:20)
cc <- FindClusters(cc, resolution = 0.09)
cc <- RunUMAP(cc, dims = 1:20, seed.use = 1, n.neighbors = 25, min.dist = 0.6)

# adding annotations for Pathology
cc@meta.data$Pathology[cc@meta.data$Patient %in% c('N1', 'N2', 'N3')] <- 'normal'
cc@meta.data$Pathology[cc@meta.data$Patient %in% c('cc4_0', 'cc5_0', 'cc6_0', 'cc7_0', 'cc8_0')] <- 'pre-treatment'
cc@meta.data$Pathology[cc@meta.data$Patient %in% c('cc4_1', 'cc5_1', 'cc6_1', 'cc7_1', 'cc8_1')] <- 'post-treatment'
cc@meta.data$Patient <- factor(cc@meta.data$Patient, levels = c('NC', 'YXJN', 'YXJT', 'cc4_0', 'cc4_1', 'cc5_0', 'cc5_1', 'cc6_0', 'cc6_1', 'cc7_0', 'cc7_1', 'cc8_0', 'cc8_1'))

DimPlot(cc, group.by = 'Pathology', cells.highlight = colnames(cc)[cc@meta.data$Pathology == 'normal'], sizes.highlight = 0.2) +
  DimPlot(cc, group.by = 'Pathology', cells.highlight = colnames(cc)[cc@meta.data$Pathology == 'pre-treatment'], cols.highlight = 'darkblue', sizes.highlight = 0.2) +
  DimPlot(cc, group.by = 'Pathology', cells.highlight = colnames(cc)[cc@meta.data$Pathology == 'post-treatment'], cols.highlight = 'purple', sizes.highlight = 0.2)

# naming main clusters
new.ident_main <- c('Mye', 'T', 'Epi', 'vCAF', 'Endo', 'Plasma', 'mCAF', 'DC', 'B')
names(new.ident_main) <- levels(Idents(cc))
cc <- RenameIdents(cc, new.ident_main)
Idents(cc) <- factor(Idents(cc), levels = c('Mye', 'T', 'Plasma', 'B', 'DC', 'Endo', 'Epi', 'vCAF', 'mCAF'))
cc@meta.data$Cluster <- Idents(cc)
cc@meta.data$Pathology <- factor(cc@meta.data$Pathology, levels = c('normal', 'pre-treatment', 'post-treatment'))
DimPlot(cc, label = TRUE, label.size = 4, cols = color.main) + DimPlot(cc, group.by = 'Patient', cols = color.sample) + 
  DimPlot(cc, group.by = 'Pathology', cols = color.pathol)

# finding markers for main clusters
DefaultAssay(cc) <- 'RNA'
cc <- ScaleData(cc, features = rownames(cc))
markers.all <- FindAllMarkers(cc)
top10.all <- markers.all %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)

avgexp.all <- AverageExpression(cc, features = rev(unique(rev(top10.all$gene))), assays = 'RNA')
avgexp.all <- t(scale(t(avgexp.all$RNA)))
annot.all <- data.frame(row.names = colnames(avgexp.all), Cluster = colnames(avgexp.all))
pheatmap(avgexp.all, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annot.all, annotation_colors = list(Cluster = color.main), show_colnames = FALSE)

# quality presentation 
cc <- PercentageFeatureSet(cc, pattern = '^MT-', col.name = 'mt_ratio')
qc.data <- as.data.frame(cc@meta.data[, c('nCount_RNA', 'nFeature_RNA', 'mt_ratio', 'Cluster', 'Patient')])
head(qc.data)
qc.data$Sample <- qc.data$Patient
ggplot(data = qc.data, aes(x = Sample, y = nCount_RNA)) + geom_boxplot(aes(fill = Sample), outlier.size = 1) +
  scale_fill_manual(values = color.sample) +
  theme(panel.background = element_blank(), axis.line = element_line(size = 0.5))
ggplot(data = qc.data, aes(x = Sample, y = nFeature_RNA)) + geom_boxplot(aes(fill = Sample), outlier.size = 1) +
  scale_fill_manual(values = color.patient) +
  theme(panel.background = element_blank(), axis.line = element_line(size = 0.5))
ggplot(data = qc.data, aes(x = Sample, y = mt_ratio)) + geom_boxplot(aes(fill = Sample), outlier.size = 1) +
  scale_fill_manual(values = color.patient) +
  theme(panel.background = element_blank(), axis.line = element_line(size = 0.5))

# R(o/e) for main clusters
levels(cc@meta.data$Cluster)
pro.cc.data <- dcast(as.data.frame(table(cc@meta.data[, c('Patient', 'Cluster')])), Patient~Cluster, value.var = 'Freq')
pro.cc.data <- pro.cc.data[4:13, ]
rownames(pro.cc.data) <- levels(pro.cc.data$Patient)[4:13]
for (i in 1:nrow(pro.cc.data)) {
  pro.cc.data[i, 'total'] <- sum(pro.cc.data[i, 2:10])
  for (j in 2:ncol(pro.cc.data)) {
    pro.cc.data[i, j] <- pro.cc.data[i, j]/pro.cc.data[i, 'total']
  }
}
pro.cc.data[c('cc4_0', 'cc5_0', 'cc6_0', 'cc7_0', 'cc8_0'), 'Pathology'] <- 'pre-treatment'
pro.cc.data[c('cc4_1', 'cc5_1', 'cc6_1', 'cc7_1', 'cc8_1'), 'Pathology'] <- 'post-treatment'
pro.cc.data$Patient <- as.character(pro.cc.data$Patient)
pro.cc.data$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
pro.cc.data$Pathology <- factor(pro.cc.data$Pathology, levels = c('pre-treatment', 'post-treatment'))

pro.cc.data.with.facs <- pro.cc.data
for (i in 1:nrow(pro.cc.data)) {
  pro.cc.data.with.facs[i, 'mCAF'] <- pro.cc.data[i, 'mCAF']/sum(pro.cc.data[i, c('mCAF', 'vCAF')]) * facs_ratio[i, 'Fib cells'] * 0.01
  pro.cc.data.with.facs[i, 'vCAF'] <- pro.cc.data[i, 'vCAF']/sum(pro.cc.data[i, c('mCAF', 'vCAF')]) * facs_ratio[i, 'Fib cells'] * 0.01
  pro.cc.data.with.facs[i, 'DC'] <- pro.cc.data[i, 'DC']/sum(pro.cc.data[i, c('DC', 'B', 'Plasma', 'T', 'Mye')]) * facs_ratio[i, 'Immune cells'] * 0.01 
  pro.cc.data.with.facs[i, 'B'] <- pro.cc.data[i, 'B']/sum(pro.cc.data[i, c('DC', 'B', 'Plasma', 'T', 'Mye')]) * facs_ratio[i, 'Immune cells'] * 0.01
  pro.cc.data.with.facs[i, 'Plasma'] <- pro.cc.data[i, 'Plasma']/sum(pro.cc.data[i, c('DC', 'B', 'Plasma', 'T', 'Mye')]) * facs_ratio[i, 'Immune cells'] * 0.01
  pro.cc.data.with.facs[i, 'T'] <- pro.cc.data[i, 'T']/sum(pro.cc.data[i, c('DC', 'B', 'Plasma', 'T', 'Mye')]) * facs_ratio[i, 'Immune cells'] * 0.01
  pro.cc.data.with.facs[i, 'Mye'] <- pro.cc.data[i, 'Mye']/sum(pro.cc.data[i, c('DC', 'B', 'Plasma', 'T', 'Mye')]) * facs_ratio[i, 'Immune cells'] * 0.01
  pro.cc.data.with.facs[i, 'Epi'] <- facs_ratio[i, 'Epi cells'] * 0.01
  pro.cc.data.with.facs[i, 'Endo'] <- facs_ratio[i, 'Endo cells'] *0.01
} 

pro.cc.data.with.facs$total_num <- apply(pro.cc.data.with.facs[, 2:10], 1, sum)
for (i in 1:5) {
  pro.data.2 <- pro.cc.data.with.facs
  total <- sum(pro.cc.data.with.facs[(2*i-1):(2*i), 'total_num'])
  for (j in 2:(ncol(pro.cc.data.with.facs)-3)) {
    sub_num <- sum(pro.cc.data.with.facs[(2*i-1):(2*i), j])
    pro.data.2[2*i-1, j] <- pro.data.2[2*i-1, 'total_num']*sub_num/total
    pro.data.2[2*i, j] <- pro.data.2[2*i, 'total_num']*sub_num/total
    
    pro.cc.data.with.facs[2*i-1, j] <- pro.cc.data.with.facs[2*i-1, j]/pro.data.2[2*i-1, j]
    pro.cc.data.with.facs[2*i, j] <- pro.cc.data.with.facs[2*i, j]/pro.data.2[2*i, j]
  }
}
pro.cc.data.with.facs[c('cc4_0', 'cc5_0', 'cc6_0', 'cc7_0', 'cc8_0'), 'Pathology'] <- 'pre-treatment'
pro.cc.data.with.facs[c('cc4_1', 'cc5_1', 'cc6_1', 'cc7_1', 'cc8_1'), 'Pathology'] <- 'post-treatment'
pro.cc.data.with.facs$Pathology <- factor(pro.cc.data.with.facs$Pathology, levels = c('pre-treatment', 'post-treatment'))
pro.cc.data.with.facs$Patient <- as.character(pro.cc.data.with.facs$Patient)
pro.cc.data.with.facs$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
pro.cc.data.with.facs.melt <- melt(pro.cc.data.with.facs, id.vars = c('Patient', 'Pathology'), measure.vars = levels(cc@meta.data$Cluster))

ggplot(data = pro.cc.data.with.facs.melt, aes(x = Pathology, y = value)) + 
  geom_point(aes(color = Patient), size = 2) + 
  geom_line(aes(group = Patient), color = 'grey') +
  facet_wrap(~variable, ncol = 5, scales = 'free') + 
  theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.line = element_line(size = 0.5), axis.ticks.x = element_blank()) + 
  geom_boxplot(aes(y = value, colour = Pathology), fill = NA, size = 0.5, outlier.colour = NA) + 
  scale_color_manual(values = c(color.patient.2, 'pre-treatment' = 'darkred', 'post-treatment' = 'darkgreen')) + 
  ylab(label = 'R(O/E)') +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), size = 0.5, test = 'wilcox.test', margin_top = 0.1, vjust = 1.6)


### for Epi -------------------------------------------------------------------------------------------------------------
epi <- cc[, Idents(cc) %in% 'Epi']

DefaultAssay(epi) <- 'integrated'
epi <- ScaleData(epi, features = rownames(epi))
epi <- RunPCA(epi)
ElbowPlot(epi, ndims = 50)
epi <- RunUMAP(epi, dims = 1:10, n.epochs = 800, n.neighbors = 50)
epi <- FindNeighbors(epi, dims = 1:10)
epi <- FindClusters(epi, resolution = 0.3)
DimPlot(epi, label = TRUE, pt.size = 1, split.by = 'Pathology', cols = color.epi)

new.ident.epi <- c('Epi1', 'Epi2', 'Epi3', 'Epi4', 'Epi5')
names(new.ident.epi) <- levels(new.ident.epi)
epi <- RenameIdents(epi, new.ident.epi)

# R(o/e) for epi
pro.epi.data <- dcast(as.data.frame(table(epi@meta.data[, c('Patient', 'Cluster')])), Patient~Cluster, value.var = 'Freq')
pro.epi.data <- pro.epi.data[4:13, ]
rownames(pro.epi.data) <- levels(epi@meta.data$Patient)[4:13]
pro.epi.data$Patient <- as.character(pro.epi.data$Patient)
pro.epi.data$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
pro.epi.data$total_num <- apply(pro.epi.data[2:6], 1, sum)
pro.epi.data <- trans2ROE(pro.epi.data)
pro.epi.data.melt <- melt(pro.epi.data, id.vars = c('Patient', 'Pathology'), measure.vars = levels(epi@meta.data$Cluster))
ggplot(data = pro.epi.data.melt, aes(x = Pathology, y = value)) + 
  geom_point(aes(color = Patient), size = 3) + 
  geom_line(aes(group = Patient), color = 'grey') +
  facet_wrap(~variable, ncol = 7, scales = 'free') + 
  theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.line = element_line(size = 0.5), axis.ticks.x = element_blank()) + 
  geom_boxplot(aes(y = value, colour = Pathology), fill = NA, size = 0.5, outlier.colour = NA) + 
  scale_color_manual(values = c(color.patient.2, 'pre-treatment' = 'darkred', 'post-treatment' = 'darkgreen')) + 
  ylab(label = 'R(O/E)') +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), size = 0.5, test = 'wilcox.test', margin_top = 0.1, vjust = 1.6)

# GO enrichment for pre and post Epi
epi.small <- epi[, epi@meta.data$Pathology %in% c('pre-treatment', 'post-treatment')]
Idents(epi.small) <- 'Pathology'
DefaultAssay(epi.small) <- 'RNA'
markers.epi.small <- FindAllMarkers(epi.small, min.pct = 0.25)

markers.epi.pre <- markers.epi.small[markers.epi.small$cluster %in% 'pre-treatment', ] %>% top_n(100, wt = avg_logFC) %>% pull(gene)
markers.epi.pre.entrez <- bitr(markers.epi.pre, fromType = 'SYMBOL', toType = 'ENTREZID', orgDb = org.Hs.eg.db) %>% pull(ENTREZID)
go.epi.pre <- enrichGO(markers.epi.pre.entrez, keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)

markers.epi.post <- markers.epi.small[markers.epi.small$cluster %in% 'post-treatment', ] %>% top_n(100, wt = avg_logFC) %>% pull(gene)
markers.epi.post.entrez <- bitr(markers.epi.post, fromType = 'SYMBOL', toType = 'ENTREZID', orgDb = org.Hs.eg.db) %>% pull(ENTREZID)
go.epi.post <- enrichGO(markers.epi.post.entrez, keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)

# enrichment score for antigen processing and presentation
genes.antpropresentation <- strsplit(go.epi.post@result$geneID[go.epi.post@result$Description %in% 'antigen processing and presentation'], split = '/')
contri.epi.antpp <- data.frame(row.names = colnames(epi.small),
                              Mean = apply(as.data.frame(t(epi.small@assays$RNA@data[genes.antpropresentation[[1]], ])), 1, mean),
                              Pathology = epi.small@meta.data$Pathology, 
                              Cluster = epi.small@meta.data$Cluster)
ggplot(contri.epi.antpp, aes(x = Pathology, y = Mean)) + geom_violin(aes(fill = Pathology)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) + facet_wrap(.~Cluster, scales = 'free') +
  scale_fill_manual(values = color.pathol[2:3]) + labs(y = 'Expression', title = 'Antigen Processiong and Presentation') +
  theme(panel.background = element_blank(), axis.line = element_line(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), vjust = 2)

# MHC signature genes for pre and post
mhc.1.genes <- c('HLA-A', 'HLA-B', 'HLA-C')
mhc.2.genes <- setdiff(rownames(epi)[substr(rownames(epi), 0, 5) %in% 'HLA-D'], 'HLA-DQB2')
epi.mhc  <- epi[c(mhc.1.genes, mhc.2.genes), epi@meta.data$Pathology %in% c('pre-treatment', 'post-treatment')] ## there's no expression of BTLA in NK
markers.epi.subcluster <- data.frame(row.names = rownames(epi.small))
for (i in 1:length(levels(epi.small@meta.data$Cluster))) {
  cluster <- levels(epi@meta.data$Cluster)[i]
  epi.new <- epi.small[, epi.small@meta.data$Cluster %in% cluster]
  Idents(epi.new) <- 'Pathology'
  markers.epi.new <- FindAllMarkers(epi.new, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 1)
  markers.epi.new <- markers.epi.new[1:(nrow(markers.epi.new)/2), ]
  markers.epi.new$bonferroni <- p.adjust(markers.epi.new$p_val, method = 'bonferroni')
  markers.epi.new$bonferroni <- -log10(markers.epi.new$bonferroni)
  for (i in 1:nrow(markers.epi.new)) {
    if (markers.epi.new[i, 'avg_log2FC'] > 0) {
      markers.epi.new[i, 'bonferroni'] <- -1*markers.epi.new[i, 'bonferroni']
    }
  }
  markers.epi.new <- markers.epi.new[c(mhc.1.genes, mhc.2.genes), ]
  markers.epi.subcluster <- cbind(markers.epi.subcluster, markers.epi.new[, 'bonferroni', drop = FALSE])
  rm(cluster, epi.new, markers.epi.new)
}
colnames(markers.epi.subcluster) <- levels(epi@meta.data$Cluster)

Idents(epi.small) <- 'Pathology'
markers.epi.new <- FindAllMarkers(epi.small, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 1)
markers.epi.new <- markers.epi.new[1:(nrow(markers.epi.new)/2), ]
markers.epi.new$bonferroni <- p.adjust(markers.epi.new$p_val, method = 'bonferroni')
markers.epi.new$bonferroni <- -log10(markers.epi.new$bonferroni)
for (i in 1:nrow(markers.epi.new)) {
  if (markers.epi.new[i, 'avg_log2FC'] > 0) {
    markers.epi.new[i, 'bonferroni'] <- -1*markers.epi.new[i, 'bonferroni']
  }
}
markers.epi.new <- markers.epi.new[c(mhc.1.genes, mhc.2.genes), ]
markers.epi <- markers.epi.new[, 'bonferroni', drop = FALSE]
colnames(markers.epi) <- 'Epi'
markers.epi <- cbind(markers.epi, markers.epi.subcluster)
markers.epi <- markers.epi[c('HLA-B', 'HLA-C', 'HLA-A', 
                             'HLA-DRB5', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA2', 'HLA-DQA1', 'HLA-DOA', 'HLA-DQB1', 'HLA-DRB1', 'HLA-DMB', 'HLA-DRA', 'HLA-DMA', 'HLA-DOB'), ]
range(markers.epi[, 2:6])
markers.epi[markers.epi > 8] <- 8
markers.epi[markers.epi < -8] <- -8
annot.epi <- data.frame(row.names = colnames(markers.epi),
                        Cluster = colnames(markers.epi))
annot.epi.row <- data.frame(row.names = rownames(markers.epi),
                            Ident = c(rep('MHC I', 3), rep('MHC II', 12)))
pheatmap(markers.epi, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annot.epi, show_colnames = FALSE, 
         color = colorRampPalette(c('#104E8B', 'white', 'firebrick3'))(200), gaps_row = 3, border_color = 'darkgrey',
         annotation_colors = list(Cluster = c('Epi' = 'navy', color.epi[1:5])), annotation_row = annot.epi.row)

# cell cycle status of Epis 	 
load('/data3/R/modules/cellcycle_genes_human.RData')
epi.new <- epi
epi.new <- CellCycleScoring(epi.new, s.features = s.genes_human, g2m.features = g2m.genes_human)
epi.new@meta.data$Phase <- factor(epi.new@meta.data$Phase, levels = c('G1', 'S', 'G2M'))
levels(epi.new@meta.data$Cluster)
epi.phase <- epi.new@meta.data[, c('Path_Cluster', 'Phase')]
epi.new@meta.data$Path_Cluster <- factor(epi.new@meta.data$Path_Cluster, levels = c(paste0(rep('normal_', 5), levels(epi@meta.data$Cluster)),
                                                                                    paste0(rep('pre-treatment_', 5), levels(epi@meta.data$Cluster)),
                                                                                    paste0(rep('post-treatment_', 5), levels(epi@meta.data$Cluster))))
ggplot(data = epi.phase, aes(fill = Phase, x = Path_Cluster)) + geom_bar(stat = 'count', position = 'fill') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), panel.background = element_blank(), panel.border = element_rect(size = 1, fill = NA)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = 'Percentage', x = '') + scale_fill_manual(values = color.cycle)

# CNV analysis for Epi
matrix.epi.cnv <- epi@assays$RNA@counts
copykat.epi <- copykat(rawmat = matrix.epi.cnv, id.type = 'S', cell.line = 'no', ngene.chr = 5, win.size = 25, sam.name = 'cc.epi',
                       KS.cut = 0.15, distance = 'euclidean', n.cores = 2)
save(copykat.epi, file = '/data3/analysis/cc/copykatepi.RData') 

predic.epi <- data.frame(copykat.epi$prediction)
epi.cnv <- epi[, rownames(predic.epi)]
epi.cnv@meta.data[, 'CNV_predict'] <- predic.epi[, 'copykat.pred']
epi.cnv@meta.data$CNV_predict <- factor(epi.cnv@meta.data$CNV_predict, levels = c('diploid', 'aneuploid'))
DimPlot(epi.cnv, group.by = 'CNV_predict', cols = c('grey', 'red'), order = TRUE) + DimPlot(epi.cnv, cols = color.epi, group.by = 'Cluster')
predic.epi$Cluster <- epi.cnv@meta.data$Cluster
predic.epi$Pathology <- epi.cnv@meta.data$Pathology
epi.cnv@meta.data$Path_cluster <- paste0(epi.cnv@meta.data$Pathology, '_', epi.cnv@meta.data$Cluster)
predic.epi$Path_cluster <- epi.cnv@meta.data$Path_cluster
Idents(epi.cnv) <- 'Path_cluster'
Idents(epi.cnv) <- factor(Idents(epi.cnv), levels = paste0(c(rep('normal_', 5), rep('pre-treatment_', 5), rep('post-treatment_', 5)),
                                                           rep(levels(epi.cnv@meta.data$Cluster), 3)))
predic.epi.cluster <- data.frame(row.names = levels(Idents(epi.cnv)), Cluster = rep(levels(epi.cnv@meta.data$Cluster), 3), 
                                 Pathology = c(rep('normal', 5), rep('pre-treatment', 5), rep('post-treatment', 5)))
predic.epi.cluster$Pathology <- factor(predic.epi.cluster$Pathology, levels = levels(epi.cnv@meta.data$Pathology))
for (i in 1:length(rownames(predic.epi.cluster))) {
  predic.epi.cluster[i, 'num_aneuploid'] <- nrow(predic.epi[predic.epi$copykat.pred %in% c('aneuploid') & predic.epi$Path_cluster %in% rownames(predic.epi.cluster)[i],])
  predic.epi.cluster[i, 'percent_aneuploid'] <- 100*predic.epi.cluster$num_aneuploid[i]/ncol(epi.cnv[, epi.cnv@meta.data$Path_cluster %in% rownames(predic.epi.cluster)[i]])
  predic.epi.cluster[i, 'percent_aneuploid'] <- round(predic.epi.cluster[i, 'percent_aneuploid'], 2)
}
ggplot(predic.epi.cluster, aes(x = Cluster, y = percent_aneuploid, fill = Cluster)) + geom_bar(stat = 'identity', width = 0.75) + 
  scale_fill_manual(values = color.epi[1:5]) + facet_grid(.~Pathology) +
  theme(panel.background = element_blank(), axis.line = element_line(color = 'black'), panel.border = element_rect(size = 1, fill = NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  geom_text(mapping = aes(label = percent_aneuploid, y = percent_aneuploid + 1)) + labs(x = element_blank())
  
### for Mye --------------------------------------------------------------------------------------------------------------------
mye.cell <- cc[, Idents(cc) %in% c('Mye', 'DC')]

DefaultAssay(mye.cell) <- 'integrated'
mye.cell <- ScaleData(mye.cell)
mye.cell <- FindVariableFeatures(mye.cell)
mye.cell <- RunPCA(mye.cell, features = VariableFeatures(mye.cell))
ElbowPlot(mye.cell, ndims = 50)
mye.cell <- RunUMAP(mye.cell, dims = 1:10, min.dist = 0.2, local.connectivity = 2)
mye.cell <- FindNeighbors(mye.cell, dims = 1:10)
mye.cell <- FindClusters(mye.cell, resolution = 0.25)
new.ident.mye <- c('IFIT2_M-MDSC', 'CCL20_Mac', 'APOE_Mac', 'FCGR3A_Mac', 'FCN1_M-MDSC', 'CD1C_DC', 'pDC')
names(new.ident.mye) <- levels(Idents(mye.cell))
mye.cell <- RenameIdents(mye.cell, new.ident.mye)

DimPlot(mye.cell, cols = color.mye, pt.size = 1) + 
  DimPlot(mye.cell, cols = color.mye, cells = colnames(mye.cell)[mye.cell@meta.data$Pathology %in% 'pre-treatment'], pt.size = 1) +
  DimPlot(mye.cell, cols = color.mye, cells = colnames(mye.cell)[mye.cell@meta.data$Pathology %in% 'post-treatment'], pt.size = 1)

# proportion change in mye.cell
mat.bar.ratio.mye <- data.frame(row.names = levels(mye.cell), Cluster = levels(mye.cell))
for (i in 1:length(levels(mye.cell))) {
  mat.bar.ratio.mye[i, 'pre'] <- nrow(mye.cell@meta.data[mye.cell@meta.data$Cluster %in% levels(mye.cell)[i] & 
                                                       mye.cell@meta.data$Pathology %in% 'pre-treatment', ])/nrow(mye.cell@meta.data[mye.cell@meta.data$Pathology %in% 'pre-treatment', ])
  mat.bar.ratio.mye[i, 'post'] <- nrow(mye.cell@meta.data[mye.cell@meta.data$Cluster %in% levels(mye.cell)[i] & 
                                                        mye.cell@meta.data$Pathology %in% 'post-treatment', ])/nrow(mye.cell@meta.data[mye.cell@meta.data$Pathology %in% 'post-treatment', ])
}
mat.bar.ratio.mye.melt <- melt(mat.bar.ratio.mye, id.vars = 'Cluster')
mat.bar.ratio.mye.melt$Cluster <- factor(mat.bar.ratio.mye.melt$Cluster, levels = levels(mye.cell))
x.end.mye <- c()
y.end.mye <- c()
for (i in 1:(length(levels(mye.cell))-1)) {
  x.end.mye[i] <- sum(mat.bar.ratio.mye$pre[(8-i):7])
  y.end.mye[i] <- sum(mat.bar.ratio.mye$post[(8-i):7])
}
ggplot(mat.bar.ratio.mye.melt, aes(x = variable, y = value)) + 
  geom_bar(aes(fill = Cluster), position = 'fill', stat = 'identity', width = 0.6) +
  scale_fill_manual(values = color.mye) +
  geom_segment(data = data.frame(), aes(x = rep(1.3, 6), y = x.end.mye, xend = rep(1.7, 6), yend = y.end.mye)) +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(x = '', y = 'Proportion', title = 'Proportion change in Mye')
  
# heatmap of MHC molecules in Mac and Mac subclusters: post vs. pre with bonferroni correction
mye.cell.mac  <- mye.cell[c(m1.signature, m2.signature, mhc.1.genes, mhc.2.genes), 
                          mye.cell@meta.data$Pathology %in% c('pre-treatment', 'post-treatment') & mye.cell@meta.data$Cluster %in% c('FCGR3A_Mac', 'APOE_Mac', 'CCL20_Mac')] 
markers.mac.subcluster <- data.frame(row.names = rownames(mye.cell.mac))
for (i in 1:3) {
  cluster <- levels(mye.cell.small@meta.data$Cluster)[i]
  mac.new <- mye.cell.mac[, mye.cell.mac@meta.data$Cluster %in% cluster]
  Idents(mac.new) <- 'Pathology'
  markers.mac.new <- FindAllMarkers(mac.new, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 1)
  markers.mac.new <- markers.mac.new[1:(nrow(markers.mac.new)/2), ]
  markers.mac.new$bonferroni <- p.adjust(markers.mac.new$p_val, method = 'bonferroni')
  markers.mac.new$bonferroni <- -log10(markers.mac.new$bonferroni)
  for (i in 1:nrow(markers.mac.new)) {
    if (markers.mac.new[i, 'avg_logFC'] > 0) {
      markers.mac.new[i, 'bonferroni'] <- -1*markers.mac.new[i, 'bonferroni']
    }
  }
  markers.mac.new <- markers.mac.new[c(m1.signature, m2.signature, mhc.1.genes, mhc.2.genes), ]
  markers.mac.subcluster <- cbind(markers.mac.subcluster, markers.mac.new[, 'bonferroni', drop = FALSE])
  rm(cluster, mac.new, markers.mac.new)
}
colnames(markers.mac.subcluster) <- levels(mye.cell.mac@meta.data$Cluster)[1:3]

Idents(mye.cell.mac) <- 'Pathology'
markers.mac.new <- FindAllMarkers(mye.cell.mac, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 1)
markers.mac.new <- markers.mac.new[1:(nrow(markers.mac.new)/2), ]
markers.mac.new$bonferroni <- p.adjust(markers.mac.new$p_val, method = 'bonferroni')
markers.mac.new$bonferroni <- -log10(markers.mac.new$bonferroni)
for (i in 1:nrow(markers.mac.new)) {
  if (markers.mac.new[i, 'avg_logFC'] > 0) {
    markers.mac.new[i, 'bonferroni'] <- -1*markers.mac.new[i, 'bonferroni']
  }
}
markers.mac.new <- markers.mac.new[rownames(markers.mac.subcluster), ]
markers.mac <- markers.mac.new[, 'bonferroni', drop = FALSE]
colnames(markers.mac) <- 'Mac'
rm(markers.mac.new)
markers.mac <- cbind(markers.mac, markers.mac.subcluster)

quantile(markers.mac[, 1])
range(markers.mac[, 2:4])
markers.mac[markers.mac > 10] <- 10
markers.mac[markers.mac < -10] <- -10
annot.mac <- data.frame(row.names = colnames(markers.mac),
                        Cluster = colnames(markers.mac))
markers.mac <- markers.mac[c(m1.signature[c(3, 5, 9, 1, 2, 8, 6, 4, 7, 10)], 
                             m2.signature[c(2, 10, 11, 3, 12, 9, 1, 14, 8, 4:5, 7, 13, 15:18, 19, 6)], 
                             mhc.1.genes[c(2, 1, 3)], 
                             mhc.2.genes[c(2, 6, 1, 3:4, 11, 12, 9, 7, 8, 10)]), ]
pheatmap(markers.mac, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annot.mac, show_colnames = FALSE, 
         color = colorRampPalette(c('#104E8B', 'white', 'firebrick3'))(100), gaps_row = c(10, 29, 32), border_color = 'darkgrey',
         annotation_colors = list(Cluster = c('Mac' = 'navy', color.mye[1:3])), annotation_row = annot.mye.row)

# R(o/e) for Mye subclusters post vs pre
pro.mye.data <- dcast(as.data.frame(table(mye.cell@meta.data[, c('Patient', 'Cluster')])), Patient~Cluster, value.var = 'Freq')
pro.mye.data <- pro.mye.data[4:13, ]
rownames(pro.mye.data) <- levels(mye.cell@meta.data$Patient)[4:13]
pro.mye.data$Patient <- as.chara`cter(pro.mye.data$Patient)
pro.mye.data$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
pro.mye.data$total_num <- apply(pro.mye.data[2:8], 1, sum)
pro.mye.data <- trans2ROE(pro.mye.data)
pro.mye.data.melt <- melt(pro.mye.data, id.vars = c('Patient', 'Pathology'), measure.vars = levels(mye.cell@meta.data$Cluster))
ggplot(data = pro.mye.data.melt, aes(x = Pathology, y = value)) + 
  geom_point(aes(color = Patient), size = 3) + 
  geom_line(aes(group = Patient), color = 'grey') +
  facet_wrap(~variable, ncol = 4, scales = 'free') + 
  theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.line = element_line(), axis.ticks.x = element_blank()) + 
  geom_boxplot(aes(y = value, colour = Pathology), fill = NA, size = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = c(color.patient.2, 'pre-treatment' = 'darkred', 'post-treatment' = 'darkgreen')) + 
  ylab(label = 'R(O/E)') +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), size = 0.5, vjust = 1.5)

# score for genesets m1/m2/mhc2
mye.cell <- AddModuleScore(mye.cell, features = list(m1.signature), name = 'M1.signature')
mye.cell <- AddModuleScore(mye.cell, features = list(m2.signature), name = 'M2.signature')
mye.cell <- AddModuleScore(mye.cell, features = list(mhc.2.genes), name = 'MHC2.signature')
mat.mye.genesets <- mye.cell@meta.data[mye.cell@meta.data$Pathology %in% c('pre-treatment', 'post-treatment'), c('M1.signature1', 'M2.signature1', 'MHC2.signature1', 'Cluster', 'Pathology')]
mat.mye.genesets.melt <- melt(mat.mye.genesets, id.vars = c('Cluster', 'Pathology'), measure.vars = c('M1.signature1', 'M2.signature1', 'MHC2.signature1'))
ggplot(mat.mye.genesets.melt, aes(x = Pathology, y = value, fill = Pathology)) + geom_violin(scale = 'width') + scale_fill_manual(values = color.pathol[2:3]) +
  geom_boxplot(width = 0.2, fill = 'white') + geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), vjust = 2, color = 'red') +
  facet_wrap(Cluster~variable, scales = 'free') + 
  theme(panel.background = element_blank(), axis.line = element_line(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = '', y = 'Score', title = 'Genesets Scores')
  
### for T&NK ----------------------------------------------------------------------------------------------------------------
t.cell <- cc[, Idents(cc) %in% 'T']
dim(t.cell)

DefaultAssay(t.cell) <- 'integrated'
t.cell <- ScaleData(t.cell)
t.cell <- FindVariableFeatures(t.cell)
t.cell <- RunPCA(t.cell, features = VariableFeatures(t.cell))
ElbowPlot(t.cell, ndims = 50)
t.cell <- FindNeighbors(t.cell, dims = 1:10)
t.cell <- FindClusters(t.cell, resolution = 0.6)
new.ident.t <- c('CD4_naive', 'CCL5_CD8_T', 'Th17', 'TNFRSF9low_Treg', 'exhausted_CD8_T', 'CD56_NK', 'IFIT_T', 'TNFRSF9high_Treg', 'GZMK_CD8_T',
                 'proliferating_T', 'CD16_NK', 'JUN_T') 
names(new.ident.t) <- levels(Idents(t.cell))
t.cell <- RenameIdents(t.cell, new.ident.t)

DimPlot(t.cell, cols = color.t, pt.size = 1) + 
  DimPlot(t.cell, cells = colnames(t.cell)[t.cell@meta.data$Pathology %in% 'pre-treatment'], cols = color.t, pt.size = 1) +
  DimPlot(t.cell, cells = colnames(t.cell)[t.cell@meta.data$Pathology %in% 'post-treatment'], cols = color.t, pt.size = 1)

# proportion change in T
mat.bar.ratio.t <- data.frame(row.names = levels(t.cell),
                                Cluster = levels(t.cell))
for (i in 1:length(levels(t.cell))) {
  mat.bar.ratio.t[i, 'pre'] <- nrow(t.cell@meta.data[t.cell@meta.data$Cluster %in% levels(t.cell)[i] & 
                                                      t.cell@meta.data$Pathology %in% 'pre-treatment', ])/nrow(t.cell@meta.data[t.cell@meta.data$Pathology %in% 'pre-treatment', ])
  mat.bar.ratio.t[i, 'post'] <- nrow(t.cell@meta.data[t.cell@meta.data$Cluster %in% levels(t.cell)[i] & 
                                                       t.cell@meta.data$Pathology %in% 'post-treatment', ])/nrow(t.cell@meta.data[t.cell@meta.data$Pathology %in% 'post-treatment', ])
}
mat.bar.ratio.t.melt <- melt(mat.bar.ratio.t, id.vars = 'Cluster')
mat.bar.ratio.t.melt$Cluster <- factor(mat.bar.ratio.t.melt$Cluster, levels = levels(t.cell))
x.end.t <- c()
y.end.t <- c()
for (i in 1:11) {
  x.end.t[i] <- sum(mat.bar.ratio.t$pre[(13-i):12])
  y.end.t[i] <- sum(mat.bar.ratio.t$post[(13-i):12])
}
ggplot(mat.bar.ratio.t.melt, aes(x = variable, y = value)) + 
  geom_bar(aes(fill = Cluster), position = 'fill', stat = 'identity', width = 0.6) +
  scale_fill_manual(values = color.t) +
  geom_segment(data = data.frame(), aes(x = rep(1.3, 11), y = x.end.t, xend = rep(1.7, 11), yend = y.end.t)) +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(x = '', y = 'Proportion', title = 'Proportion change in T')
  
# GO for NK
t.nk.small <- t.cell[, Idents(t.cell) %in% c('CD16_NK', 'CD56_NK') & t.nk.small@meta.data$Pathology %in%c('pre-treatment', 'post-treatment')]
Idents(t.nk.small) <- 'Pathology'
markers.nk <- FindAllMarkers(t.nk.small)
markers.nk.pathol <- markers.nk %>% group_by(cluster)
nrow(markers.nk.pathol)
markers.nk.pathol <- markers.nk.pathol[2456:4910, ]
markers.nk.pathol$change <- ifelse(markers.nk.pathol$p_val < 1e-35 & abs(markers.nk.pathol$avg_logFC) > 0.5, 
                                   ifelse(markers.nk.pathol$avg_logFC > 0.5, 'Up', 'Down'), 'Stable')
ggplot(markers.nk.pathol, aes(x = avg_logFC, y = -log10(p_val), color = change)) + geom_point(size = 1.5) + 
  geom_text_repel(data = subset(markers.nk.pathol, !markers.nk.pathol$change %in% 'Stable'), aes(label = gene), 
                  size = 3, max.overlaps = 40, show.legend = FALSE) +
  scale_color_manual(values = list('Up' = 'forestgreen', 'Down' = 'firebrick3', 'Stable' = 'grey')) +
  geom_vline(xintercept = c(-0.5, 0.5), lty = 2) + geom_hline(yintercept = 35, lty = 2) + 
  labs(title = 'Post vs. Pre for NK') + theme(plot.title = element_text(hjust = 0.5, face = 'bold')) + #ylim(0, 200) +
  theme(panel.background = element_blank(), axis.line = element_line()) + guides()

markers.nk.pre <- markers.nk.small[markers.nk.small$cluster %in% 'pre-treatment', ] %>% top_n(100, wt = avg_logFC) %>% pull(gene)
markers.nk.pre.entrez <- bitr(markers.nk.pre, fromType = 'SYMBOL', toType = 'ENTREZID', orgDb = org.Hs.eg.db) %>% pull(ENTREZID)
go.nk.pre <- enrichGO(markers.nk.pre.entrez, keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)

markers.nk.post <- markers.nk.small[markers.nk.small$cluster %in% 'post-treatment', ] %>% top_n(100, wt = avg_logFC) %>% pull(gene)
markers.nk.post.entrez <- bitr(markers.nk.post, fromType = 'SYMBOL', toType = 'ENTREZID', orgDb = org.Hs.eg.db) %>% pull(ENTREZID)
go.nk.post <- enrichGO(markers.nk.post.entrez, keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)

# cytotoxicity&inhibitory signature genes in NK
cytotoxic.markers <- c('NKG7', 'GZMH', 'GNLY', 'IFNG', 'GZMB', 'GZMK', 'PRF1', 'GZMA') ## reorder the genes
inhibitory.markers <- c('CTLA4', 'HAVCR2', 'TIGIT', 'LAG3', 'PDCD1', 'BTLA') ## reorder the genes
Idents(t.cell) <- 'Pathology'
t.small.temp  <- t.cell[c(cytotoxic.markers, inhibitory.markers), 
                              (Idents(t.cell) %in% c('pre-treatment', 'post-treatment')) & (!t.cell@meta.data$Cluster %in% c('CD16_NK', 'CD56_NK'))]
markers.t.subcluster <- data.frame(row.names = rownames(t.small.temp))
levels(t.small.temp@meta.data$Cluster)
for (i in 3:12) {
  cluster <- levels(t.small.temp@meta.data$Cluster)[i]
  nk.new <- t.small.temp[, t.small.temp@meta.data$Cluster %in% cluster]
  markers.t.new <- FindAllMarkers(nk.new, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 2)
  markers.t.new <- markers.t.new[1:(nrow(markers.t.new)/2), ]
  markers.t.new$bonferroni <- p.adjust(markers.t.new$p_val, method = 'bonferroni')
  markers.t.new$bonferroni <- -log10(markers.t.new$bonferroni)
  for (i in 1:nrow(markers.t.new)) {
    if (markers.t.new[i, 'avg_logFC'] > 0) {
      markers.t.new[i, 'bonferroni'] <- -1*markers.t.new[i, 'bonferroni']
    }
  }
  markers.t.new <- markers.t.new[c(cytotoxic.markers, inhibitory.markers), ]
  markers.t.subcluster <- cbind(markers.t.subcluster, markers.t.new[, 'bonferroni', drop = FALSE])
  rm(cluster, nk.new, markers.t.new)
}
colnames(markers.t.subcluster) <- levels(Idents(t.small.temp))[3:12]

DefaultAssay(t.small.temp) <- 'RNA'
markers.t.except.nk <- data.frame(row.names = rownames(t.small.temp))
markers.t.except.nk <- FindAllMarkers(t.small.temp, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 2)
markers.t.except.nk <- markers.t.except.nk[1:(nrow(markers.t.except.nk)/2), ]
markers.t.except.nk$bonferroni <- p.adjust(markers.t.except.nk$p_val, method = 'bonferroni')
markers.t.except.nk$bonferroni <- -log10(markers.t.except.nk$bonferroni)
for (i in 1:nrow(markers.t.except.nk)) {
  if (markers.t.except.nk[i, 'avg_logFC'] > 0) {
    markers.t.except.nk[i, 'bonferroni'] <- -1*markers.t.except.nk[i, 'bonferroni']
  }
}
markers.t.except.nk <- markers.t.except.nk[, 'bonferroni', drop = FALSE]
markers.t.except.nk <- markers.t.except.nk[c(cytotoxic.markers, inhibitory.markers), , drop = FALSE]
markers.t.except.nk <- cbind(markers.t.except.nk, markers.t.subcluster)

colnames(markers.t.except.nk) <- c('T', levels(t.small.temp@meta.data$Cluster)[3:12])
range(markers.t.except.nk[, 2:11])  

markers.t.except.nk[markers.t.except.nk > 5] <- 5
markers.t.except.nk[markers.t.except.nk < -5] <- -5
annot.t <- data.frame(row.names = colnames(markers.t.except.nk),
                       Cluster = colnames(markers.t.except.nk))
annot.t.row <- data.frame(row.names = rownames(markers.t.except.nk),
                          Ident = c(rep('Cytoxicity Signature', 8), rep('Inhibitory Signature', 6)))
pheatmap(markers.t.except.nk, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annot.t, show_colnames = FALSE, 
         color = colorRampPalette(c('#104E8B', 'white', 'firebrick3'))(200), gaps_row = 8, border_color = 'darkgrey',
         annotation_colors = list(Cluster = c('T' = 'navy', color.t[3:12])), annotation_row = annot.t.row)

# survival analysis for NK subpopulation markers
caselist <- getCaseLists(CGDS('http://www.cbioportal.org/'), cancerStudy = 'cesc_tcga')[1, 1]
datatable(getGeneticProfiles(CGDS('http://www.cbioportal.org/'), 'cesc_tcga'))
genetic.pro.name <- getGeneticProfiles(CGDS('http://www.cbioportal.org/'), 'cesc_tcga')[4, 1]
DT::datatable(getCancerStudies(CGDS('http://www.cbioportal.org/', 'cesc_tcga')))

genetic.pro.data.prf1 <- getProfileData(CGDS('http://www.cbioportal.org/'),
                                       genes = 'GZMB',
                                       geneticProfiles = genetic.pro.name, caseList = caselist)
genetic.pro.data.prf1 <- na.omit(genetic.pro.data.prf1)

clinical.data <- getClinicalData(CGDS('http://www.cbioportal.org/'), caseList = caselist)
clinical.data <- data.frame(row.names = rownames(clinical.data), OS_MONTHS = clinical.data$OS_MONTHS, OS_STATUS = clinical.data$OS_STATUS)
clinical.data <- clinical.data[rownames(genetic.pro.data.prf1), ]

clinical.data <- cbind(clinical.data, PRF1 = genetic.pro.data.prf1[, 'GZMB', drop = FALSE])
head(clinical.data)
clinical.data$OS_STATUS <- as.character(clinical.data$OS_STATUS)
clinical.data$OS_STATUS[clinical.data$OS_STATUS %in% '1:DECEASED'] <- '1'
clinical.data$OS_STATUS[clinical.data$OS_STATUS %in% '0:LIVING'] <- '0'
clinical.data$OS_STATUS <- as.factor(clinical.data$OS_STATUS)

sroc <- survivalROC::survivalROC(Stime = clinical.data$OS_MONTHS, status = clinical.data$OS_STATUS,
                                 marker = clinical.data$GZMB, predict.time = 1825, method = 'KM')
head(sroc)
sroc$cut.values[which.max(sroc$TP- sroc$FP)]

clinical.data$level[clinical.data$GZMB >= sroc$cut.values[which.max(sroc$TP- sroc$FP)]] <- 'High'
clinical.data$level[clinical.data$GZMB < sroc$cut.values[which.max(sroc$TP- sroc$FP)]] <- 'Low'
quantile(clinical.data$GZMB)

attach(clinical.data)
surv.prf1 <- Surv(OS_MONTHS, OS_STATUS == '1')
surv.prf1_level <- survfit(surv.prf1~level)
ggsurvplot(surv.prf1_level, conf.int = F, risk.table = FALSE, ncenso.plot = TRUE, data = clinical.data, pval = TRUE, palette = c('red', 'blue'),
           ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line())) + 
  labs(title = 'GZMB')
detach(clinical.data)

# R(o/e)
pro.t.data <- dcast(as.data.frame(table(t.cell@meta.data[, c('Patient', 'Cluster')])), Patient~Cluster, value.var = 'Freq')
pro.t.data <- pro.t.data[4:13, ]
rownames(pro.t.data) <- levels(t.cell@meta.data$Patient)[4:13]
pro.t.data$Patient <- as.character(pro.t.data$Patient)
pro.t.data$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
pro.t.data$total_num <- apply(pro.t.data[2:13], 1, sum)
pro.t.data <- trans2ROE(pro.t.data)
pro.t.data.melt <- melt(pro.t.data, id.vars = c('Patient', 'Pathology'), measure.vars = levels(t.cell@meta.data$Cluster))
ggplot(data = pro.t.data.melt, aes(x = Pathology, y = value)) + 
  geom_point(aes(color = Patient), size = 3) + 
  geom_line(aes(group = Patient), color = 'grey') +
  facet_wrap(~variable, ncol = 6, scales = 'free') + 
  theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.line = element_line(), axis.ticks.x = element_blank()) + 
  geom_boxplot(aes(y = value, colour = Pathology), fill = NA, outlier.shape = NA) + 
  scale_color_manual(values = c(color.patient.2, 'pre-treatment' = 'darkred', 'post-treatment' = 'darkgreen')) + 
  ylab(label = 'R(O/E)') +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), size = 0.5, vjust = 1.5)

# expression of inhibitory and cytotoxic genesets in T except NKs:post vs. pre
cytotoxic.markers <- c('NKG7', 'GZMH', 'GNLY', 'IFNG', 'GZMB', 'GZMK', 'PRF1', 'GZMA') ## reorder the genes
inhibitory.markers <- c('CTLA4', 'HAVCR2', 'TIGIT', 'LAG3', 'PDCD1', 'BTLA') ## reorder the genes
Idents(t.cell) <- 'Pathology'
t.small.temp  <- t.cell[c(cytotoxic.markers, inhibitory.markers), 
                              (Idents(t.cell) %in% c('pre-treatment', 'post-treatment')) & (!t.cell@meta.data$Cluster %in% c('CD16_NK', 'CD56_NK'))]
markers.t.subcluster <- data.frame(row.names = rownames(t.small.temp))
levels(t.small.temp@meta.data$Cluster)
for (i in 3:12) {
  cluster <- levels(t.small.temp@meta.data$Cluster)[i]
  nk.new <- t.small.temp[, t.small.temp@meta.data$Cluster %in% cluster]
  markers.t.new <- FindAllMarkers(nk.new, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 2)
  markers.t.new <- markers.t.new[1:(nrow(markers.t.new)/2), ]
  markers.t.new$bonferroni <- p.adjust(markers.t.new$p_val, method = 'bonferroni')
  markers.t.new$bonferroni <- -log10(markers.t.new$bonferroni)
  for (i in 1:nrow(markers.t.new)) {
    if (markers.t.new[i, 'avg_logFC'] > 0) {
      markers.t.new[i, 'bonferroni'] <- -1*markers.t.new[i, 'bonferroni']
    }
  }
  markers.t.new <- markers.t.new[c(cytotoxic.markers, inhibitory.markers), ]
  markers.t.subcluster <- cbind(markers.t.subcluster, markers.t.new[, 'bonferroni', drop = FALSE])
  rm(cluster, nk.new, markers.t.new)
}
colnames(markers.t.subcluster) <- levels(Idents(t.small.temp))[3:12]

markers.t.except.nk <- data.frame(row.names = rownames(t.small.temp))
markers.t.except.nk <- FindAllMarkers(t.small.temp, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 2)
markers.t.except.nk <- markers.t.except.nk[1:(nrow(markers.t.except.nk)/2), ]
markers.t.except.nk$bonferroni <- p.adjust(markers.t.except.nk$p_val, method = 'bonferroni')
markers.t.except.nk$bonferroni <- -log10(markers.t.except.nk$bonferroni)
for (i in 1:nrow(markers.t.except.nk)) {
  if (markers.t.except.nk[i, 'avg_logFC'] > 0) {
    markers.t.except.nk[i, 'bonferroni'] <- -1*markers.t.except.nk[i, 'bonferroni']
  }
}
markers.t.except.nk <- markers.t.except.nk[, 'bonferroni', drop = FALSE]
markers.t.except.nk <- markers.t.except.nk[c(cytotoxic.markers, inhibitory.markers), , drop = FALSE]
markers.t.except.nk <- cbind(markers.t.except.nk, markers.t.subcluster)

colnames(markers.t.except.nk) <- c('T', levels(t.small.temp@meta.data$Cluster)[3:12])
range(markers.t.except.nk[, 2:11])  

markers.t.except.nk[markers.t.except.nk > 5] <- 5
markers.t.except.nk[markers.t.except.nk < -5] <- -5
annot.t <- data.frame(row.names = colnames(markers.t.except.nk),
                       Cluster = colnames(markers.t.except.nk))
annot.t.row <- data.frame(row.names = rownames(markers.t.except.nk),
                          Ident = c(rep('Cytoxicity Signature', 8), rep('Inhibitory Signature', 6)))
pheatmap(markers.t.except.nk, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annot.t, show_colnames = FALSE, 
         color = colorRampPalette(c('#104E8B', 'white', 'firebrick3'))(200), gaps_row = 8, border_color = 'darkgrey',
         annotation_colors = list(Cluster = c('T' = 'navy', color.t[3:12])), annotation_row = annot.t.row)

### for CAF ---------------------------------------------------------------------------------
fib <- cc[, Idents(cc) %in% c('mCAF', 'vCAF')]
dim(fib) # 10278 cells

DefaultAssay(fib) <- 'integrated'
fib <- ScaleData(fib, features = rownames(fib))
fib <- RunPCA(fib)
ElbowPlot(fib, ndims = 50)
fib <- RunUMAP(fib, dims = 1:15, seed.use = 1, min.dist = 0.2)
fib <- FindNeighbors(fib, dims = 1:15)
fib <- FindClusters(fib, resolution = 0.6)
fib.ident <- c('mCAF', 'apCAF2', 'vCAF2', 'vCAF1', 'apCAF1', 'apCAF1', 'vCAF2', 'vCAF2', 'iCAF')
names(fib.ident) <- levels(Idents(fib))
fib <- RenameIdents(fib, fib.ident)
Idents(fib) <- factor(Idents(fib), levels = c('mCAF', 'iCAF', 'apCAF1', 'apCAF2', 'vCAF1', 'vCAF2'))
fib@meta.data$Cluster <- Idents(fib)

# proportion change in CAF
mat.bar.ratio.caf <- data.frame(row.names = levels(fib),
                                Cluster = levels(fib))
for (i in 1:length(levels(fib))) {
  mat.bar.ratio.caf[i, 'pre'] <- nrow(fib@meta.data[fib@meta.data$Cluster %in% levels(fib)[i] & 
                                                           fib@meta.data$Pathology %in% 'pre-treatment', ])/nrow(fib@meta.data[fib@meta.data$Pathology %in% 'pre-treatment', ])
  mat.bar.ratio.caf[i, 'post'] <- nrow(fib@meta.data[fib@meta.data$Cluster %in% levels(fib)[i] & 
                                                            fib@meta.data$Pathology %in% 'post-treatment', ])/nrow(fib@meta.data[fib@meta.data$Pathology %in% 'post-treatment', ])
}
mat.bar.ratio.caf.melt <- melt(mat.bar.ratio.caf, id.vars = 'Cluster')
mat.bar.ratio.caf.melt$Cluster <- factor(mat.bar.ratio.caf.melt$Cluster, levels = levels(fib))
x.end.caf <- c()
y.end.caf <- c()
for (i in 1:(length(levels(fib))-1)) {
  x.end.caf[i] <- sum(mat.bar.ratio.caf$pre[(7-i):6])
  y.end.caf[i] <- sum(mat.bar.ratio.caf$post[(7-i):6])
}
ggplot(mat.bar.ratio.caf.melt, aes(x = variable, y = value)) + 
  geom_bar(aes(fill = Cluster), position = 'fill', stat = 'identity', width = 0.6) +
  scale_fill_manual(values = color.fib) +
  geom_segment(data = data.frame(), aes(x = rep(1.3, 5), y = x.end.caf, xend = rep(1.7, 5), yend = y.end.caf)) +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(x = '', y = 'Proportion', title = 'Proportion change in CAF')
  
# R(o/e) for CAF
pro.caf.data <- dcast(as.data.frame(table(fib@meta.data[, c('Patient', 'Cluster')])), Patient~Cluster, value.var = 'Freq')
pro.caf.data <- pro.caf.data[4:13, ]
rownames(pro.caf.data) <- levels(fib@meta.data$Patient)[4:13]
pro.caf.data$Patient <- as.character(pro.caf.data$Patient)
pro.caf.data$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
pro.caf.data$total_num <- apply(pro.caf.data[2:7], 1, sum)
pro.caf.data <- trans2ROE(pro.caf.data)
pro.caf.data.melt <- melt(pro.caf.data, id.vars = c('Patient', 'Pathology'), measure.vars = levels(fib@meta.data$Cluster))
ggplot(data = pro.caf.data.melt, aes(x = Pathology, y = value)) + 
  geom_point(aes(color = Patient), size = 3) + 
  geom_line(aes(group = Patient), color = 'grey') +
  facet_wrap(~variable, ncol = 3, scales = 'free') + 
  theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.line = element_line(), axis.ticks.x = element_blank()) + 
  geom_boxplot(aes(y = value, colour = Pathology), fill = NA, size = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = c(color.patient.2, 'pre-treatment' = 'darkred', 'post-treatment' = 'darkgreen')) + 
  ylab(label = 'R(O/E)') +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), size = 0.5, vjust = 1.5)
  
# Heatmap of genesets in CAF
angio.genes <- c('VEGFA', 'PDGFA', 'PDGFB', 'ANGPT2', 'EGFL6')
inflam.genes <- c('CXCL1', 'CXCL14', 'IL6', 'CCL2', 'CXCL16', 'CXCL12', 'CXCL2', 'CXCL9', 'CXCL10', 'IFNG')
mhc.2.genes <- c('HLA-DPB1', 'HLA-DMA', 'HLA-DRB1', 'HLA-DRA', 'HLA-DQA1', 'HLA-DMB', 'HLA-DQB1', 'HLA-DOA', 'HLA-DQA2', 'HLA-DPA1', 'HLA-DRB5', 'HLA-DOB')
collagen.genes <- c('COL1A1', 'COL1A2', 'COL3A1', 'COL5A1', 'COL5A2', 'COL6A1', 'COL6A2', 'COL6A3', 'COL4A1')
avgexpre.caf.sets <- AverageExpression(fib, features = unique(intersect(c(collagen.genes, inflam.genes, mhc.2.genes,angio.genes), rownames(fib))),
                                       assays = 'RNA')
avgexpre.caf.sets$RNA <- t(scale(t(avgexpre.caf.sets$RNA)))
pheatmap(avgexpre.caf.sets$RNA, cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, annotation_col = annot.caf,
         annotation_colors = list(Cluster = color.fib), color = colorRampPalette(c('#104E8B','white', 'orange', 'firebrick3'))(100),
         gaps_row = c(9, 19, 31))
		 
# Heatmap for bonferroni for CAF post vs. pre
fib.small  <- fib[c(mhc.1.genes, mhc.2.genes), fib@meta.data$Pathology %in% c('pre-treatment', 'post-treatment')] ## there's no expression of BTLA in NK
markers.fib.subcluster <- data.frame(row.names = rownames(fib.small))
for (i in 1:length(levels(fib.small@meta.data$Cluster))) {
  cluster <- levels(fib@meta.data$Cluster)[i]
  fib.new <- fib.small[, fib.small@meta.data$Cluster %in% cluster]
  Idents(fib.new) <- 'Pathology'
  markers.fib.new <- FindAllMarkers(fib.new, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 1)
  markers.fib.new <- markers.fib.new[1:(nrow(markers.fib.new)/2), ]
  markers.fib.new$bonferroni <- p.adjust(markers.fib.new$p_val, method = 'bonferroni')
  markers.fib.new$bonferroni <- -log10(markers.fib.new$bonferroni)
  for (i in 1:nrow(markers.fib.new)) {
    if (markers.fib.new[i, 'avg_logFC'] > 0) {
      markers.fib.new[i, 'bonferroni'] <- -1*markers.fib.new[i, 'bonferroni']
    }
  }
  markers.fib.new <- markers.fib.new[c(mhc.1.genes, mhc.2.genes), ]
  markers.fib.subcluster <- cbind(markers.fib.subcluster, markers.fib.new[, 'bonferroni', drop = FALSE])
  rm(cluster, fib.new, markers.fib.new)
}
colnames(markers.fib.subcluster) <- levels(fib@meta.data$Cluster)

Idents(fib.small) <- 'Pathology'
markers.fib.new <- FindAllMarkers(fib.small, min.pct = 0, logfc.threshold = 0, assay = 'RNA', return.thresh = 1)
markers.fib.new <- markers.fib.new[1:(nrow(markers.fib.new)/2), ]
markers.fib.new$bonferroni <- p.adjust(markers.fib.new$p_val, method = 'bonferroni')
markers.fib.new$bonferroni <- -log10(markers.fib.new$bonferroni)
for (i in 1:nrow(markers.fib.new)) {
  if (markers.fib.new[i, 'avg_logFC'] > 0) {
    markers.fib.new[i, 'bonferroni'] <- -1*markers.fib.new[i, 'bonferroni']
  }
}
markers.fib.new <- markers.fib.new[c(mhc.1.genes, mhc.2.genes), ]
markers.fib <- markers.fib.new[, 'bonferroni', drop = FALSE]
colnames(markers.fib) <- 'fib'
rm(markers.fib.new)
markers.fib <- cbind(markers.fib, markers.fib.subcluster)
markers.fib <- markers.fib[c('HLA-B', 'HLA-C', 'HLA-A', 
                             'HLA-DRB5', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA2', 'HLA-DQA1', 'HLA-DOA', 'HLA-DQB1', 'HLA-DRB1', 'HLA-DMB', 'HLA-DRA', 'HLA-DMA', 'HLA-DOB'), ]
range(markers.fib[, 2:6])
markers.fib[markers.fib > 10] <- 10
markers.fib[markers.fib < -10] <- -10
annot.fib <- data.frame(row.names = colnames(markers.fib),
                        Cluster = colnames(markers.fib))
annot.fib.row <- data.frame(row.names = rownames(markers.fib),
                            Ident = c(rep('MHC I', 3), rep('MHC II', 12)))
pheatmap(markers.fib, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = annot.fib, show_colnames = FALSE, 
         color = colorRampPalette(c('#104E8B', 'white', 'firebrick3'))(200), gaps_row = 3, border_color = 'darkgrey',
         annotation_colors = list(Cluster = c('fib' = 'navy', color.fib[1:5])), annotation_row = annot.fib.row)
		 
# GSVA comparison for pre and post CAF
original_gmt_cancer_gsva <- readLines('/data3/R/files/h.all.v7.4.entrez.gmt')
database_list_cancer_gsva <- lapply(original_gmt_cancer_gsva, strsplit_no_name)
for (i in 1:length(database_list_cancer_gsva)) {
  names(database_list_cancer_gsva)[i] <- database_list_cancer_gsva[i][[1]][1]
  database_list_cancer_gsva[i][[1]] <- database_list_cancer_gsva[i][[1]][-1]
}

fib.small <- fib[, fib@meta.data$Pathology %in% c('pre-treatment', 'post-treatment')]
Idents(fib.small) <- 'Cluster'
caf.sample <- fib.small[, sample(colnames(fib.small), size = ncol(fib.small)/2)]
Idents(caf.sample) <- 'Pathology'
counts_caf <- caf.sample@assays$RNA@counts 
counts_caf <- counts_caf + 1

ENTREZ <- bitr(rownames(counts_caf), fromType = 'SYMBOL', toType = c('ENTREZID'), OrgDb = org.Hs.eg.db)
colnames(ENTREZ) <- c('gene.name', 'ENTREZID')
counts_caf <- data.frame(counts_caf)
counts_caf$gene.name <- rownames(counts_caf)
counts_caf <- inner_join(ENTREZ, counts_caf, by = 'gene.name')
counts_caf[1:10, 1:10]
rownames(counts_caf) <- counts_caf$ENTREZID
dim(counts_caf)
colnames(counts_caf) <- gsub(colnames(counts_caf), pattern = '\\.', replacement = '_')
counts_caf <- counts_caf[, c(-1, -2)]

es.caf <- gsva(as.matrix(counts_caf), database_list_cancer_gsva, mx.diff = FALSE, kcdf = 'Poisson', parallel.sz = 1)
colnames(es.caf) <- colnames(caf.sample)
 
rownames(es.caf) <- substr(rownames(es.caf), 10, 100)
t.caf <- data.frame(row.names = rownames(es.caf),
                    Genesets = rownames(es.caf))
for (i in 1:nrow(es.caf)) {
  t.test <- t.test(es.caf[i, colnames(caf.sample)[caf.sample@meta.data$Pathology %in% 'post-treatment']],
                   es.caf[i, colnames(caf.sample)[caf.sample@meta.data$Pathology %in% 'pre-treatment']])
  t.caf$t_value[i] <- t.test$statistic
  t.caf$p_value[i] <- t.test$p.value
}
t.caf$Genesets <- factor(t.caf$Genesets, levels = rownames(t.caf[order(t.caf$t_value, decreasing = FALSE), ]))
t.caf[order(t.caf$t_value, decreasing = FALSE), ]
t.caf$Annot[t.caf$t_value > 0] <- 'Post'
t.caf$Annot[t.caf$t_value < 0] <- 'Pre'
t.caf$Annot[t.caf$p_value >= 0.05] <- 'Unsignificant'
ggplot(data = t.caf, aes(x = t_value, y = Genesets)) + geom_bar(aes(fill = Annot), width =0.8, stat = 'identity') + 
  scale_fill_manual(values = c(Pre = 'firebrick3', Post = 'forestgreen', Unisignificant = 'grey')) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, size = 1.5))

### for EC ---------------------------------------------------------------------------------------
endo <- cc[, Idents(cc) %in% c('Endo')]

DefaultAssay(endo) <- 'integrated'
endo <- ScaleData(endo) 
endo <- RunPCA(endo)
ElbowPlot(endo, ndims = 50)
endo <- FindNeighbors(endo, dims = 1:10, do.plot = TRUE)
endo <- FindClusters(endo, resolution = 0.25)
endo <- RunUMAP(endo, dims = 1:10, n.neighbors = 30, min.dist = 0.3)

new.ident.ec <- c(paste0(rep('EC', 5), 1:5), 'EC5')
names(new.ident.ec) <- levels(Idents(endo))
endo <- RenameIdents(endo, new.ident.ec)
DimPlot(endo, group.by = 'Cluster', cols = color.ec, pt.size = 1) + 
  DimPlot(endo, cells = colnames(endo)[endo@meta.data$Pathology %in% 'pre-treatment'], cols = color.ec, pt.size = 1) +
  DimPlot(endo, cells = colnames(endo)[endo@meta.data$Pathology %in% 'post-treatment'], cols = color.ec, pt.size = 1)

DefaultAssay(endo) <- 'RNA'
endo <- ScaleData(endo, features = rownames(endo))
markers.ec <- FindAllMarkers(endo, min.pct = 0.25)
markers.ec.50 <- markers.ec %>% group_by(cluster) %>% top_n(50, wt = avg_logFC)

gene_list_0 <- markers.ec.50[markers.ec.50$cluster %in% '0',]$gene
gene_0 <- bitr(gene_list_0, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db, drop = TRUE)
gene_0 <- gene_0$ENTREZID
ego_0 <- enrichGO(gene = gene_0, keyType = 'ENTREZID', OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)
ego_0_simp <- simplify(ego_0, cutoff =0.7, select_fun = min)
barplot(ego_0_simp, showCategory = 10) + aes(y = log10(pvalue), x = Description) + scale_y_reverse()

# proportion change in EC
mat.bar.ratio.ec <- data.frame(row.names = levels(endo),
                                Cluster = levels(endo))
for (i in 1:length(levels(endo))) {
  mat.bar.ratio.ec[i, 'pre'] <- nrow(endo@meta.data[endo@meta.data$Cluster %in% levels(endo)[i] & 
                                                           endo@meta.data$Pathology %in% 'pre-treatment', ])/nrow(endo@meta.data[endo@meta.data$Pathology %in% 'pre-treatment', ])
  mat.bar.ratio.ec[i, 'post'] <- nrow(endo@meta.data[endo@meta.data$Cluster %in% levels(endo)[i] & 
                                                            endo@meta.data$Pathology %in% 'post-treatment', ])/nrow(endo@meta.data[endo@meta.data$Pathology %in% 'post-treatment', ])
}
mat.bar.ratio.ec.melt <- melt(mat.bar.ratio.ec, id.vars = 'Cluster')
mat.bar.ratio.ec.melt$Cluster <- factor(mat.bar.ratio.ec.melt$Cluster, levels = levels(endo))
x.end.ec <- c()
y.end.ec <- c()
for (i in 1:(length(levels(endo))-1)) {
  x.end.ec[i] <- sum(mat.bar.ratio.ec$pre[(6-i):5])
  y.end.ec[i] <- sum(mat.bar.ratio.ec$post[(6-i):5])
}
ggplot(mat.bar.ratio.ec.melt, aes(x = variable, y = value)) + 
  geom_bar(aes(fill = Cluster), position = 'fill', stat = 'identity', width = 0.6) +
  scale_fill_manual(values = color.ec) +
  geom_segment(data = data.frame(), aes(x = rep(1.3, 4), y = x.end.ec, xend = rep(1.7, 4), yend = y.end.ec)) +
  theme(panel.background = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(x = '', y = 'Proportion', title = 'Proportion change in EC')
  
# GO term enrichment for EC
gene.list <- markers.ec %>% group_by(cluster) %>% top_n(50, wt = avg_logFC) %>% subset(cluster %in% 'EC5') %>% pull(gene)
entrez <- bitr(gene.list, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db, drop = TRUE)
genes <- entrez$ENTREZID
ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)
ego.simp <- simplify(ego, cutoff = 0.7, select_fun = min)
barplot(ego.simp, showCategory = 10) + aes(y = log10(pvalue), x = Description) + scale_y_reverse() +
  theme(panel.background = element_blank(), panel.grid = element_blank())
  
# R(o/e) for EC subpopulations
pro.ec.data <- dcast(as.data.frame(table(endo@meta.data[, c('Patient', 'Cluster')])), Patient~Cluster, value.var = 'Freq')
pro.ec.data <- pro.ec.data[4:13, ]
rownames(pro.ec.data) <- levels(endo@meta.data$Patient)[4:13]
pro.ec.data$Patient <- as.character(pro.ec.data$Patient)
pro.ec.data$Patient <- c(rep('cc4', 2), rep('cc5', 2), rep('cc6', 2), rep('cc7', 2), rep('cc8', 2))
pro.ec.data$total_num <- apply(pro.ec.data[, 2:6], 1, sum)
pro.ec.data <- trans2ROE(pro.ec.data)
pro.ec.data.melt <- melt(pro.ec.data, id.vars = c('Patient', 'Pathology'), measure.vars = levels(endo@meta.data$Cluster))
ggplot(data = pro.ec.data.melt, aes(x = Pathology, y = value)) + 
  geom_point(aes(color = Patient), size = 3) + 
  geom_line(aes(group = Patient), color = 'grey') +
  facet_wrap(~variable, ncol = 5, scales = 'free') + 
  theme(axis.text.x = element_blank(), panel.background = element_blank(), axis.line = element_line(), axis.ticks.x = element_blank()) + 
  geom_boxplot(aes(y = value, colour = Pathology), fill = NA, size = 0.5, outlier.shape = NA) + 
  scale_color_manual(values = c(color.patient.2, 'pre-treatment' = 'darkred', 'post-treatment' = 'darkgreen')) + 
  ylab(label = 'R(O/E)') +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), size = 0.5, vjust = 2)

# GSVA comparison for pre and post EC
endo.small <- endo[, endo@meta.data$Pathology %in% c('pre-treatment', 'post-treatment')]
Idents(endo.small) <- 'Pathology'
counts_endo <- endo.small@assays$RNA@counts 
counts_endo <- counts_endo + 1

ENTREZ <- bitr(rownames(counts_endo), fromType = 'SYMBOL', toType = c('ENTREZID'), OrgDb = org.Hs.eg.db)
colnames(ENTREZ) <- c('gene.name', 'ENTREZID')
counts_endo <- data.frame(counts_endo)
counts_endo$gene.name <- rownames(counts_endo)
counts_endo <- inner_join(ENTREZ, counts_endo, by = 'gene.name')
counts_endo[1:10, 1:10]
rownames(counts_endo) <- counts_endo$ENTREZID
dim(counts_endo)
colnames(counts_endo) <- gsub(colnames(counts_endo), pattern = '\\.', replacement = '_')
counts_endo <- counts_endo[, c(-1, -2)]

es.endo <- gsva(as.matrix(counts_endo), database_list_cancer_gsva, mx.diff = FALSE, kcdf = 'Poisson', parallel.sz = 1)
colnames(es.endo) <- colnames(endo.small)
save(endo.small, es.endo, file = '/data3/analysis/cc/es_endo.RData')
load('/data3/analysis/cc/es_endo.RData')

rownames(es.endo) <- substr(rownames(es.endo), 10, 100)
t.endo <- data.frame(row.names = rownames(es.endo),
                    Genesets = rownames(es.endo))
for (i in 1:nrow(es.endo)) {
  t.test <- t.test(es.endo[i, colnames(endo.small)[endo.small@meta.data$Pathology %in% 'post-treatment']],
                   es.endo[i, colnames(endo.small)[endo.small@meta.data$Pathology %in% 'pre-treatment']])
  t.endo$t_value[i] <- t.test$statistic
  t.endo$p_value[i] <- t.test$p.value
}
t.endo$Genesets <- factor(t.endo$Genesets, levels = rownames(t.endo[order(t.endo$t_value, decreasing = FALSE), ]))
t.endo[order(t.endo$t_value, decreasing = FALSE), ]
t.endo$Annot[t.endo$t_value > 0] <- 'Post'
t.endo$Annot[t.endo$t_value < 0] <- 'Pre'
t.endo$Annot[t.endo$p_value >= 0.05] <- 'Unsignificant'
ggplot(data = t.endo, aes(x = t_value, y = Genesets)) + geom_bar(aes(fill = Annot), width =0.8, stat = 'identity') + 
  scale_fill_manual(values = c(Pre = 'firebrick3', Post = 'forestgreen', Unisignificant = 'grey')) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, size = 1.5)) 
mean(es.endo['COMPLEMENT', colnames(endo.small)[endo.small@meta.data$Pathology %in% 'post-treatment']])

# cluster contribution by EC subclusters
epi.pro.genes <- strsplit(go.filtered.ec$geneID[go.filtered.ec$Description %in% 'epithelial cell proliferation'], split = '/')[[1]]
ec.dev.genes <-strsplit(go.filtered.ec$geneID[go.filtered.ec$Description %in% 'endothelial cell development'], split = '/')[[1]] 
blood.ves.genes <- strsplit(go.filtered.ec$geneID[go.filtered.ec$Description %in% 'blood vessel remodeling'], split = '/')[[1]]
ec.pro.genes <- strsplit(go.filtered.ec$geneID[go.filtered.ec$Description %in% 'endothelial cell proliferation'], split = '/')[[1]]
hif1.sig.genes <- strsplit(go.filtered.ec$geneID[go.filtered.ec$Description %in% 'HIF-1 signaling pathway'], split = '/')[[1]]
focal.adh.genes <- strsplit(go.filtered.ec$geneID[go.filtered.ec$Description %in% 'Focal adhesion'], split = '/')[[1]]

endo.small <- endo[focal.adh.genes, endo@meta.data$Pathology %in% c('pre-treatment', 'post-treatment')]
mat.endo <- data.frame(row.names = colnames(endo.small),
                       Expression = apply(t(endo.small@assays$RNA@data), 1, mean),
                       Pathology = endo.small@meta.data$Pathology,
                       Cluster = endo.small@meta.data$Cluster)
mat.endo.melt <- melt(mat.endo, id.vars = c('Cluster', 'Pathology'))
ggplot(mat.endo.melt, aes(x = Pathology, y = value)) + geom_violin(aes(fill = Pathology)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) + scale_fill_manual(values = color.pathol[2:3]) +
  facet_wrap(~Cluster, ncol = 5, scales = 'free') +
  theme(panel.background = element_blank(), axis.line = element_line(), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c('pre-treatment', 'post-treatment')), vjust = 2) +
  labs(title = 'Focal Adhesion', y = 'Expression')
  
### CIBERSORTx for hazard ration -------------------------------------------------------------------------
clinical.data <- read.csv('E:/BioInfor learning/cc/CIBERSORT/clinical_data.csv', header=T, row.names=1)
ciber.result <- read.table('E:/BioInfor learning/cc/CIBERSORT/CIBERSORT-Results.txt', header=T, row.names=1)
clinical.data$OS_STATUS <- ifelse(clinical.data$OS_STATUS %in% '0:LIVING', 0, 1)

dim(ciber.result); dim(clinical.data)
ciber.result <- apply(ciber.result, 2, function(x) ifelse(x>median(x), 1, 0))
clinical.data$OS_STATUS <- as.numeric(clinical.data$OS_STATUS)

clinical.data <- cbind(clinical.data[, c('OS_MONTHS', 'OS_STATUS')], ciber.result[, 1:37])
res.cox <- coxph(Surv(OS_MONTHS, OS_STATUS)~TNFRSF9low_Treg+CCL20_Mac+APOE_Mac+apCAF1+FCN1_M.MDSC+CD1C_DC+FCGR3A_Mac+GZMK_CD8_T+CD16_NK+IFIT2_M.MDSC+B+vCAF2+CD4_naive+JUN_T+
                   CCL5_CD8_T+vCAF1+Plasma+Epi2+pDC+CD56_NK+TNFRSF9high_Treg+Epi1+mCAF+exhausted_CD8_T+IFIT_T+EC1+EC4+Th17+proliferating_T+EC5+Epi5+Epi4+apCAF2+EC3+EC2+Epi3+iCAF,
                 data = clinical.data)
ggforest(res.cox)

### cell interaction -------------------------------------------------------------------------------------
library(cellcall)

load(file = '/data3/analysis/cc/cc_integrate.RData')
load('/data3/analysis/cc/colors.RData')
load('/data3/analysis/cc/markers.RData')
DimPlot(cc, group.by = 'Cluster')

cc@meta.data$Cluster_inter <- cc@meta.data$Cluster
cc@meta.data$Cluster_inter <- as.character(cc@meta.data$Cluster_inter)
cc@meta.data$Cluster_inter[cc@meta.data$Subcluster %in% c('FCN1_M-MDSC', 'IFIT2_M-MDSC')] <- 'MDSC'
cc@meta.data$Cluster_inter[cc@meta.data$Subcluster %in% c('APOC_Mac', 'CCL20_Mac', 'FCGR3A_Mac')] <- 'Mac'
cc@meta.data$Cluster_inter[cc@meta.data$Subcluster %in% c('vCAF', 'mCAF')] <- 'Fibro'
cc@meta.data$Cluster_inter[cc@meta.data$Subcluster %in% c('CD16_NK', 'CD56_NK')] <- 'NK'
cc@meta.data$Cluster_inter <- as.factor(cc@meta.data$Cluster_inter)
cc.pre <- cc[, cc@meta.data$Pathology %in% 'pre-treatment']
cc.post <- cc[, cc@meta.data$Pathology %in% 'post-treatment']

levels(cc@meta.data$Cluster_inter)
cc.epi.others.pre <- cc.pre[, cc.pre@meta.data$Cluster_inter %in% c('Epi', 'T', 'NK', 'MDSC', 'Mac', 'DC')]
cc.epi.others.post <- cc.post[, cc.post@meta.data$Cluster_inter %in% c('Epi', 'T', 'NK', 'MDSC', 'Mac', 'DC')]

test.epi.pre <- CreateObj1ect_fromSeurat(Seurat.object= cc.epi.others.pre, 
                                slot="counts", 
                                cell_type="Cluster_inter",
                                data_source="UMI",
                                scale.factor = 10^6, 
                                Org = "Homo sapiens")
test.epi.pre <- TransCommuProfile(object = test.epi.pre,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')

test.epi.post <- CreateObject_fromSeurat(Seurat.object= cc.epi.others.post, 
                                         slot="counts", 
                                         cell_type="Cluster_inter",
                                         data_source="UMI",
                                         scale.factor = 10^6, 
                                         Org = "Homo sapiens")
test.epi.post <- TransCommuProfile(object = test.epi.post,
                                   pValueCor = 0.05,
                                   CorValue = 0.1,
                                   topTargetCor=1,
                                   p.adjust = 0.05,
                                   use.type="median",
                                   probs = 0.9,
                                   method="weighted",
                                   IS_core = TRUE,
                                   Org = 'Homo sapiens')


pair.epi.pre <- test.epi.pre@data$expr_l_r_log2_scale
pair.epi.post <- test.epi.post@data$expr_l_r_log2_scale

pathway.hyper.list.epi.pre <- lapply(colnames(pair.epi.pre), function(i){
  print(i)
  tmp <- getHyperPathway(data = pair.epi.pre, object = test.epi.pre, cella_cellb = i, Org="Homo sapiens")
  return(tmp)
})

pathway.hyper.list.epi.post <- lapply(colnames(pair.epi.post), function(i){
  print(i)
  tmp <- getHyperPathway(data = pair.epi.post, object = test.epi.post, cella_cellb = i, Org="Homo sapiens")
  return(tmp)
})

save(test.epi.pre, test.epi.post, pathway.hyper.list.epi.pre, pathway.hyper.list.epi.post,
     file = '/data3/analysis/cc/cellcall_cc_epi.RData')

colnames(pair.epi.pre)
pheatmap(pair.epi.pre, cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0)
chemo.pairs.epi.pre <- c(rownames(pair.epi.pre)[grep(substr(rownames(pair.epi.pre), 1, 3), pattern = c('CXC'))],
                 rownames(pair.epi.pre)[grep(substr(rownames(pair.epi.pre), 1, 3), pattern = c('CCL'))])
pheatmap(pair.epi.pre[chemo.pairs, ], cluster_rows = TRUE, cluster_cols = FALSE)

pheatmap(pair.epi.post[, colnames(pair.epi.pre)], cluster_cols = FALSE, cluster_rows = TRUE, treeheight_row = 0)
chemo.pairs.epi.post <- c(rownames(pair.epi.post)[grep(substr(rownames(pair.epi.post), 1, 3), pattern = c('CXC'))],
                           rownames(pair.epi.post)[grep(substr(rownames(pair.epi.post), 1, 3), pattern = c('CCL'))])
pheatmap(pair.epi.post[chemo.pairs.epi.post, colnames(pair.epi.pre)], cluster_rows = TRUE, cluster_cols = FALSE)

unique(cc.epi.others.pre@meta.data$Cluster_inter)
color.int <- data.frame(color = c(color.main[7], brewer.pal('Set2', n = 5)),
                        row.names = unique(cc.epi.others.pre@meta.data$Cluster_inter))

ViewInterCircos(object = test.epi.pre, 
                font = 2, cellColor = color.int, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 0.5, #细胞类型多的话设置小点，不然图太大画不出来
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)

start <- character()
end <- character()
for (i in 1:36) {
  x <- strsplit(colnames(pair.epi.pre), split = '-')[[i]][1]
  y <- strsplit(colnames(pair.epi.pre), split = '-')[[i]][2]
  start <- c(start, x)
  end <- c(end, y)
}
mat.circo.epi.pre <- data.frame(start = start, end = end)
mat.circo.epi.pre$value <- colSums(pair.epi.pre > 0)
chordDiagram(mat.circo.epi.pre)
circos.clear()
grid.col <- brewer.pal('Dark2', n = 6)
names(grid.col) <- unique(cc.epi.others.pre@meta.data$Cluster_inter)
mat.col <- matrix(nrow = 36, ncol = 1, data = c(grid.col, rep('#00000000', 30)))
chordDiagram(mat.circo.epi.pre,
             scale = TRUE, link.target.prop = FALSE, 
             grid.col = grid.col, col = mat.col, 
             annotationTrack = c("name", 'grid'), 
             #col = hcl.colors(15), 
             transparency = 0.5,   
             directional = 1,      
             #link.lwd = 2,         
             #link.lty = 2,         
             #link.border = 0
             )  
mat.circo.epi.post <- data.frame(start = start, end = end)
mat.circo.epi.post$value <- colSums(pair.epi.post > 0)
chordDiagram(mat.circo.epi.post)
circos.clear()
grid.col <- brewer.pal('Dark2', n = 6)
names(grid.col) <- unique(cc.epi.others.pre@meta.data$Cluster_inter)
mat.col <- matrix(nrow = 36, ncol = 1, data = c(grid.col, rep('#00000000', 30)))
chordDiagram(mat.circo.epi.post,
             scale = TRUE, link.target.prop = FALSE, 
             grid.col = grid.col, col = mat.col, 
             annotationTrack = c("name", 'grid'), 
             #col = hcl.colors(15),
             transparency = 0.5,   
             directional = 1,      
             #link.lwd = 2,         
             #link.lty = 2,         
             #link.border = 0
)  
circos.clear()

