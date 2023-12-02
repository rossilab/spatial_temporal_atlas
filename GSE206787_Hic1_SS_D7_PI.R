################################# Loading R packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


################################# Loading gene count matrices
hic1_undamaged <- readRDS(file.choose()) # Load doublet-removed object
hic1_undamaged@meta.data$group <- "undamaged"

hic1_day7PI <- Read10X(data.dir = "INSERT CORRESPONDING DIRECTORY")
colnames(x = hic1_day7PI) <- paste('hic1_day7PI', colnames(x = hic1_day7PI), sep = '_')
hic1_day7PI <- CreateSeuratObject(counts = hic1_day7PI, min.cells = 3, min.features = 200, project = "hic1_day7PI")
hic1_day7PI@meta.data$group <- "7d PI"


################################# Quality control and Normalization
hic1_undamaged[["percent.mt"]] <- PercentageFeatureSet(hic1_undamaged, pattern = "^mt-")
hic1_day7PI[["percent.mt"]] <- PercentageFeatureSet(hic1_day7PI, pattern = "^mt-")
# VlnPlot(hic1_undamaged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(hic1_day7PI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hic1_undamaged <- subset(hic1_undamaged, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA < 15000)
hic1_day7PI <- subset(hic1_day7PI, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA < 10000)
hic1_undamaged <- NormalizeData(hic1_undamaged, normalization.method = "LogNormalize", scale.factor = 10000)
hic1_day7PI <- NormalizeData(hic1_day7PI, normalization.method = "LogNormalize", scale.factor = 10000)
hic1_undamaged <- FindVariableFeatures(hic1_undamaged, selection.method = "vst", nfeatures = 2000)
hic1_day7PI <- FindVariableFeatures(hic1_day7PI, selection.method = "vst", nfeatures = 2000)
allData <- list(hic1_undamaged, hic1_day7PI)


################################# Data integration
anchors <- FindIntegrationAnchors(object.list = allData, dims = 1:30)
combinedData <- IntegrateData(anchorset = anchors, dims = 1:30)
combinedData <- ScaleData(combinedData, features = rownames(combinedData), verbose = T, vars.to.regress = c("percent.mt"))


################################# PCA, clustering, and UMAP projection
combinedData <- RunPCA(combinedData, features = VariableFeatures(object = combinedData), npcs = 30)
# ElbowPlot(combinedData, 30)
combinedData <- FindNeighbors(combinedData, dims = 1:20)
combinedData <- FindClusters(combinedData, resolution = 0.6)
combinedData <- RunUMAP(combinedData, dims = 1:20)


################################# Saving processed Seurat object & aesthetics
saveRDS(combinedData, file = "seurat_obj.rds")
cols_timepoints <- c("#081D58", "#E31A1C")


################################# Cluster annotations
new.cluster.ids <- c("Fibroblast", "Mural cell", "Fibroblast", "Mural cell", "Mural cell", "Fibroblast", "Mural cell", "Fibroblast",
                     "Mural cell", "Mural cell", "Mural cell")
names(new.cluster.ids) <- levels(combinedData)
combinedData <- RenameIdents(combinedData, new.cluster.ids)
combinedData@active.ident <- factor(combinedData@active.ident, levels = c("Fibroblast", "Mural cell"))
combinedData@meta.data$cellgroup <- Idents(combinedData)


################################# Visualizations
# UMAP
Idents(combinedData) <- combinedData@meta.data$timepoint
SS_barcodes <- rownames(combinedData@meta.data[which(combinedData@meta.data$timepoint == "SS"),])
D7_barcodes <- rownames(combinedData@meta.data[which(combinedData@meta.data$timepoint == "7d PI"),])
UMAP_annot_SS <- DimPlot(combinedData, reduction = "umap", label = F, pt.size = 0.5, cols.highlight = "#081D58",
                         cells.highlight = SS_barcodes) + theme(text = element_text(face = "bold", size=16), axis.text=element_text(size=6), axis.title.x=element_blank(),
                                                                axis.title.y=element_blank(), legend.position = "left") + NoAxes()
UMAP_annot_D7 <- DimPlot(combinedData, reduction = "umap", label = F, pt.size = 0.5, cols.highlight = "#E31A1C",
                         cells.highlight = D7_barcodes) + theme(text = element_text(face = "bold", size=16), axis.text=element_text(size=6), axis.title.x=element_blank(),
                                                                axis.title.y=element_blank(), legend.position = "left") + NoAxes()

# Violin plots
combinedData@meta.data$timepoint <- "SS"
combinedData@meta.data$timepoint[which(combinedData@meta.data$group == "7d PI")] <- "D7"
combinedData@meta.data$timepoint <- factor(combinedData@meta.data$timepoint, levels = c("SS","D7"))
Idents(combinedData) <- combinedData@meta.data$cellgroup
VlnPlot(object = combinedData, features = "GENE", same.y.lims = TRUE, pt.size = 0.1, sort = F, cols = cols_timepoints, split.by = "timepoint", assay = "RNA") +
  theme(text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(vjust=-1, face = "italic")) + NoLegend()
