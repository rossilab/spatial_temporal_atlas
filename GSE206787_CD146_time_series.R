################################# Loading R packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


################################# Remove doublets using hashing data in D3 dataset
# Loading gene/HTO count matrix
mcam_day3PI <- Read10X(data.dir = "INSERT CORRESPONDING DIRECTORY") # Reading data into expression matrix
colnames(x = mcam_day3PI$`Gene Expression`) <- paste('MCAM_3dpi', colnames(x = mcam_day3PI$`Gene Expression`), sep = '_')
colnames(x = mcam_day3PI$`Antibody Capture`) <- paste('MCAM_3dpi', colnames(x = mcam_day3PI$`Antibody Capture`), sep = '_')
exprData <- CreateSeuratObject(counts = mcam_day3PI$`Gene Expression`, min.cells = 3, min.features = 200, project = "MCAM_3dpi")
sharedCells <- intersect(colnames(exprData), colnames(mcam_day3PI$`Antibody Capture`))
HTOData <- mcam_day3PI$`Antibody Capture`[, sharedCells]

# Processing HTO data
exprData[["HTO"]] <- CreateAssayObject(counts = HTOData)
exprData <- NormalizeData(exprData, assay = "HTO", normalization.method = "CLR")
exprData <- HTODemux(exprData, assay = "HTO", positive.quantile = 0.99)
# Summary of annotations
table(exprData$HTO_classification.global)

# Remove doublets
Idents(exprData) <- "HTO_classification.global"
exprData <- subset(exprData, idents = c("Singlet", "Negative"))
exprData@active.assay <- "RNA"
Idents(exprData) <- "orig.ident"

# Saving processed Seurat object
saveRDS(exprData, file = "HTO_processed_object.rds")


################################# Loading gene count matrices
mcam_undamaged <- Read10X(data.dir = "INSERT CORRESPONDING DIRECTORY")
colnames(x = mcam_undamaged) <- paste('mcam_undamaged', colnames(x = mcam_undamaged), sep = '_')
mcam_undamaged <- CreateSeuratObject(counts = mcam_undamaged, min.cells = 3, min.features = 200, project = "mcam_undamaged")
mcam_undamaged@meta.data$group <- "undamaged"

mcam_day3PI <- readRDS(file = file.choose())
mcam_day3PI@meta.data$group <- "3d PI"

mcam_day7PI <- Read10X(data.dir = "INSERT CORRESPONDING DIRECTORY")
colnames(x = mcam_day7PI) <- paste('mcam_day7PI', colnames(x = mcam_day7PI), sep = '_')
mcam_day7PI <- CreateSeuratObject(counts = mcam_day7PI, min.cells = 3, min.features = 200, project = "mcam_day7PI")
mcam_day7PI@meta.data$group <- "7d PI"

mcam_day42PI <- Read10X(data.dir = "INSERT CORRESPONDING DIRECTORY")
colnames(x = mcam_day42PI) <- paste('mcam_day42PI', colnames(x = mcam_day42PI), sep = '_')
mcam_day42PI <- CreateSeuratObject(counts = mcam_day42PI, min.cells = 3, min.features = 200, project = "mcam_day42PI")
mcam_day42PI@meta.data$group <- "42d PI"


################################# Quality control and Normalization
mcam_undamaged[["percent.mt"]] <- PercentageFeatureSet(mcam_undamaged, pattern = "^mt-")
mcam_day3PI[["percent.mt"]] <- PercentageFeatureSet(mcam_day3PI, pattern = "^mt-")
mcam_day7PI[["percent.mt"]] <- PercentageFeatureSet(mcam_day7PI, pattern = "^mt-")
mcam_day42PI[["percent.mt"]] <- PercentageFeatureSet(mcam_day42PI, pattern = "^mt-")

VlnPlot(mcam_undamaged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mcam_day3PI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mcam_day7PI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mcam_day42PI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mcam_undamaged <- subset(mcam_undamaged, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA < 8000)
mcam_day3PI <- subset(mcam_day3PI, subset = nFeature_RNA > 200 & percent.mt < 15 & nCount_RNA < 30000)
mcam_day7PI <- subset(mcam_day7PI, subset = nFeature_RNA > 200 & percent.mt < 15 & nCount_RNA < 10000)
mcam_day42PI <- subset(mcam_day42PI, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA < 7500)

mcam_undamaged <- NormalizeData(mcam_undamaged, normalization.method = "LogNormalize", scale.factor = 10000)
mcam_day3PI <- NormalizeData(mcam_day3PI, normalization.method = "LogNormalize", scale.factor = 10000)
mcam_day7PI <- NormalizeData(mcam_day7PI, normalization.method = "LogNormalize", scale.factor = 10000)
mcam_day42PI <- NormalizeData(mcam_day42PI, normalization.method = "LogNormalize", scale.factor = 10000)

mcam_undamaged <- FindVariableFeatures(mcam_undamaged, selection.method = "vst", nfeatures = 2000)
mcam_day3PI <- FindVariableFeatures(mcam_day3PI, selection.method = "vst", nfeatures = 2000)
mcam_day7PI <- FindVariableFeatures(mcam_day7PI, selection.method = "vst", nfeatures = 2000)
mcam_day42PI <- FindVariableFeatures(mcam_day42PI, selection.method = "vst", nfeatures = 2000)


################################# Data integration
anchors <- FindIntegrationAnchors(object.list = allData, dims = 1:30)
combinedData <- IntegrateData(anchorset = anchors, dims = 1:30)
combinedData <- ScaleData(combinedData, features = rownames(combinedData), verbose = T, vars.to.regress = c("percent.mt"))


################################# PCA, clustering, and UMAP projection
combinedData <- RunPCA(combinedData, features = VariableFeatures(object = combinedData), npcs = 30)
# ElbowPlot(combinedData, 30)
combinedData <- FindNeighbors(combinedData, dims = 1:20)
combinedData <- FindClusters(combinedData, resolution = 0.45)
combinedData <- RunUMAP(combinedData, dims = 1:20)


################################# Saving processed Seurat object & aesthetics
saveRDS(combinedData, file = "seurat_only.rds")
cols_annot <- c("#FD8D3C", "#41B6C4", "#7FCDBB", "#081D58", "#225EA8", "#E31A1C", "black")
cols_peri <- c("#081D58", "#225EA8", "#E31A1C")
cols_timepoints <- c("#081D58", "#225EA8", "#E31A1C", "black")
cols_clusters <- brewer.pal(12,"Paired")


################################# Cluster annotations
new.cluster.ids <- c("PER-1", "Ven. VSMC", "Ven. VSMC", "Ven. VSMC", "Ven. VSMC", "Art. VSMC", "Act. PER", "PER-2", "Schwann",
                     "Ven. VSMC", "Fibroblast")
names(new.cluster.ids) <- levels(combinedData)
combinedData <- RenameIdents(combinedData, new.cluster.ids)
combinedData@active.ident <- factor(combinedData@active.ident, levels = c("Fibroblast", "PER-1", "PER-2", "Act. PER", "Ven. VSMC", "Art. VSMC", "Schwann"))
combinedData@meta.data$celltype <- combinedData@active.ident


################################# Visualizations
# UMAP - cell annotations
Idents(combinedData) <- combinedData@meta.data$time
combinedData@meta.data$celltype <- factor(combinedData@meta.data$celltype, levels = c("Fibroblast", "Ven. VSMC", "Art. VSMC", "PER-1", "PER-2", "Act. PER", "Schwann"))
subset_SS <- subset(combinedData, idents = "SS")
subset_D3 <- subset(combinedData, idents = "D3")
subset_D7 <- subset(combinedData, idents = "D7")
subset_D42 <- subset(combinedData, idents = "D42")
Idents(subset_SS) <- subset_SS@meta.data$celltype
Idents(subset_D3) <- subset_D3@meta.data$celltype
Idents(subset_D7) <- subset_D7@meta.data$celltype
Idents(subset_D42) <- subset_D42@meta.data$celltype
UMAP_annot <- DimPlot(subset_SS, reduction = "umap", label = F, pt.size = 0.5, cols = cols_annot) +
  theme(text = element_text(face = "bold", size=16), axis.text=element_text(size=6), axis.title.x=element_blank(),
        axis.title.y=element_blank(), legend.position = "left") + NoAxes() # change object accordingly

# UMAP - clusters
Idents(combinedData) <- combinedData@meta.data$time
subset_SS <- subset(combinedData, idents = "SS")
subset_D3 <- subset(combinedData, idents = "D3")
subset_D7 <- subset(combinedData, idents = "D7")
subset_D42 <- subset(combinedData, idents = "D42")
Idents(subset_SS) <- subset_SS@meta.data$seurat_clusters
Idents(subset_D3) <- subset_D3@meta.data$seurat_clusters
Idents(subset_D7) <- subset_D7@meta.data$seurat_clusters
Idents(subset_D42) <- subset_D42@meta.data$seurat_clusters
UMAP_annot <- DimPlot(subset_D42, reduction = "umap", label = F, pt.size = 0.5, cols = cols_clusters) +
  theme(text = element_text(face = "bold", size=16), axis.text=element_text(size=6), axis.title.x=element_blank(),
        axis.title.y=element_blank(), legend.position = "left") + NoAxes() # change object accordingly

# Violin plots
VlnPlot(object = combinedData, features = "GENE", same.y.lims = TRUE, pt.size = 0.1, sort = F, cols = cols_clusters) +
  theme(text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(vjust=-1, face = "italic")) + NoLegend()

# Violin plots - PER only
Idents(combinedData) <- combinedData@meta.data$celltype
subsetData <- subset(combinedData, idents = c("PER-1", "PER-2", "Act. PER"))
subsetData@active.ident <- factor(subsetData@active.ident, levels = c("PER-1", "PER-2", "Act. PER"))
VlnPlot(object = subsetData, features = "GENE", same.y.lims = TRUE, pt.size = 0.1, sort = F, cols = cols_peri) +
  theme(text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        plot.title=element_text(vjust=-1, face = "italic")) + NoLegend()


################################# Pearson's correlation analysis
Idents(combinedData) <- combinedData@meta.data$seurat_clusters
cluster.averages <- AverageExpression(object = combinedData, return.seurat = TRUE, verbose = T, assays = "RNA")
mat <- as.matrix(cluster.averages@assays$RNA@data)
cormat <- as.data.frame(cor(mat), na.rm = T)
colors <- colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(255)
p <- pheatmap(cormat,
              main = NA,
              treeheight_row = 20, treeheight_col = 20,
              color = colors,
              angle_col = 45,
              show_colnames = T,
              show_rownames  = F,
              fontsize = 12)


################################# Temporal analysis - PER only
Idents(combinedData) <- combinedData@meta.data$celltype
subsettedData <- SubsetData(object = combinedData, ident.use = c("PER-1", "PER-2", "Act. PER"))
Idents(subsettedData) <- subsettedData@meta.data$time

# Find all DEGs relative to D0
clusterMarkers_D3 <- FindMarkers(subsettedData, ident.1 = c("D3"), ident.2 = c("SS"), min.pct = 0.25, only.pos = F)
clusterMarkers_D3 <- clusterMarkers_D3[which(clusterMarkers_D3$p_val_adj < 0.05),]
clusterMarkers_D7 <- FindMarkers(subsettedData, ident.1 = c("D7"), ident.2 = c("SS"), min.pct = 0.25, only.pos = F)
clusterMarkers_D7 <- clusterMarkers_D7[which(clusterMarkers_D7$p_val_adj < 0.05),]
clusterMarkers_D42 <- FindMarkers(subsettedData, ident.1 = c("D42"), ident.2 = c("SS"), min.pct = 0.25, only.pos = F)
clusterMarkers_D42 <- clusterMarkers_D42[which(clusterMarkers_D42$p_val_adj < 0.05),]
DEGs_D3 <- rownames(clusterMarkers_D3)
DEGs_D7 <- rownames(clusterMarkers_D7)
DEGs_D42 <- rownames(clusterMarkers_D42)

# Load pathways of interest
fibrosis_genes <- read.delim(file.choose(), header = F)
fibrosis_genes <- as.character(fibrosis_genes$V1)
angiogenesis_genes <- read.delim(file.choose(), header = F)
angiogenesis_genes <- as.character(angiogenesis_genes$V1)
inflam_genes <- read.delim(file.choose(), header = F)
inflam_genes <- as.character(inflam_genes$V1)

# Intersect all temporally regulated genes with specific pathways
all_markers <- unique(c(DEGs_D3, DEGs_D7, DEGs_D42))
fibro_genes <- intersect(fibrosis_genes,all_markers)
angio_genes <- intersect(angiogenesis_genes,all_markers)
inflammatory_genes <- intersect(inflam_genes, all_markers)

# Scale all genes and compute average scaled expression
subsettedData <- ScaleData(subsettedData, features = rownames(subsettedData), verbose = T)
cluster.averages <- AverageExpression(object = subsettedData, return.seurat = F, verbose = T, assays = "RNA", use.scale = T)
mat <- cluster.averages$RNA

# Subset scaled matrices of temporally regulated genes associated with each pathway
fibro_mat <- mat[fibro_genes,]
angio_mat <- mat[angio_genes,]
inflam_mat <- mat[inflammatory_genes,]
fibro_mat$genes <- rownames(fibro_mat)
angio_mat$genes <- rownames(angio_mat)
inflam_mat$genes <- rownames(inflam_mat)
fibro_mat_long <- melt(fibro_mat, id.vars = c("genes"))
angio_mat_long <- melt(angio_mat, id.vars = c("genes"))
inflam_mat_long <- melt(inflam_mat, id.vars = c("genes"))
colnames(fibro_mat_long) <- c("Gene", "Time", "Expression")
colnames(angio_mat_long) <- c("Gene", "Time", "Expression")
colnames(inflam_mat_long) <- c("Gene", "Time", "Expression")

# Plotting trends
g <- ggplot(data = fibro_mat_long, 
            aes(x = Time, y = Expression)) +
  geom_point(size = 3, color = "grey") +
  stat_summary(fun.y=mean, colour="red", geom="line", size = 1.1, aes(group = 1), linetype = "dashed") +
  ylab("Relative expression") +
  xlab("") +
  theme_bw() +
  theme(text=element_text(size=20), legend.position = "none")
g + stat_summary(fun="mean", geom="point", size=5, shape = 18,
                 aes(group = 1), fill="black")

# Heatmap of selected genes
subsettedData@meta.data$time <- factor(subsettedData@meta.data$time, levels = c("SS", "D3", "D7", "D42"))
DE_PER_ECM_final <- c("Col1a1", "Col3a1", "Eln", "Fn1", "Adamts4", "Mmp2", "Mmp11", "Loxl1", "Loxl2", "Ccn1", "Ccn2", "Postn", "Sparc", "Pdgfa")
DE_PER_vascular_final <- c("Ets1", "Meox2", "Prrx1", "Tbx20", "Apold1", "Thy1", "Actg1", "Myh9", "Mylk", "Angpt2", "Jag1", "Notch1", "Thbs2",  "Ednra", "Ednrb")
DE_PER_inflam_final <- c("Casp12", "Cfh", "Cxcl1", "Hmgb2", "Il1r1", "Il6", "Nfkbia", "Nfkbib", "Tnfrsf1a")
DE_VSMC_final <- c("Adamts4", "Adamts9", "Ccn2", "Apold1", "Actg1", "Myh9", "Mylk", "Thbs1", "Vegfb", "Egln1", "Hbegf", "Jag1", "Cxcl1",
                   "Hmgb2", "Nfkbia", "Nfkbib", "Tnfrsf1a")
cluster.averages <- AverageExpression(object = subsettedData, return.seurat = TRUE, verbose = T, assays = "RNA")
mat <- as.matrix(cluster.averages@assays$RNA@data)
cluster.averages@meta.data$annot <- c("SS", "D3", "D7", "D42")
colData <- data.frame(cluster.averages@meta.data[,4])
rownames(colData) <- rownames(cluster.averages@meta.data)
colnames(colData) <- c("annot")
cols_annot <- list(annot = c("SS" = "#80CDC1", "D3" = "#CB181D", "D7" = "#1B7837",
                             "D42" = "#542788"))
df <- mat[DE_PER_ECM_final,]
df <- mat[DE_PER_vascular_final,]
df <- mat[DE_PER_inflam_final,]
df <- mat[DE_VSMC_final,]
breaksList = seq(-2, 2, length.out = 255)
colors <- colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(255)
p <- pheatmap(df,
              main = NA,
              treeheight_row = 30, treeheight_col = 15,
              color = colors,
              angle_col = 45,
              show_colnames = T,
              show_rownames  = T,
              clustering_method = "complete",
              scale = "row",
              annotation_col = colData,
              cluster_cols = F,
              cluster_rows = F,
              annotation_names_col = F,
              annotation_colors = cols_annot,
              breaks = breaksList,
              fontsize = 36,
              annotation_legend = F,
              legend = F)
