################################# Loading R packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


################################# Loading gene count matrices
cspg4_mcam_day7PI <- Read10X(data.dir = "INSERT CORRESPONDING DIRECTORY")
cspg4_mcam_day7PI <- CreateSeuratObject(counts = cspg4_mcam_day7PI, min.cells = 3, min.features = 200, project = "cspg4_mcam_day7PI")


################################# Quality control, Normalization, and Scaling
cspg4_mcam_day7PI[["percent.mt"]] <- PercentageFeatureSet(cspg4_mcam_day7PI, pattern = "^mt-")
# VlnPlot(cspg4_mcam_day7PI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cspg4_mcam_day7PI <- subset(cspg4_mcam_day7PI, subset = nFeature_RNA > 200 & percent.mt < 15 & nCount_RNA < 20000 & nCount_RNA > 2000)
cspg4_mcam_day7PI <- NormalizeData(cspg4_mcam_day7PI, normalization.method = "LogNormalize", scale.factor = 10000)
cspg4_mcam_day7PI <- FindVariableFeatures(cspg4_mcam_day7PI, selection.method = "vst", nfeatures = 2000)
cspg4_mcam_day7PI <- ScaleData(cspg4_mcam_day7PI, vars.to.regress = "percent.mt", features = rownames(cspg4_mcam_day7PI))


################################# PCA, clustering, and UMAP projection
cspg4_mcam_day7PI <- RunPCA(cspg4_mcam_day7PI, features = VariableFeatures(object = cspg4_mcam_day7PI), npcs = 50)
# ElbowPlot(cspg4_mcam_day7PI, ndims = 50)
cspg4_mcam_day7PI <- FindNeighbors(cspg4_mcam_day7PI, dims = 1:20)
cspg4_mcam_day7PI <- FindClusters(cspg4_mcam_day7PI, resolution = 1.2)
cspg4_mcam_day7PI <- RunUMAP(cspg4_mcam_day7PI, dims = 1:20)


################################# Saving processed Seurat object & aesthetics
saveRDS(cspg4_mcam_day7PI, file = "seurat_obj.rds")
cspg4_mcam_day7PI <- readRDS(file = file.choose())
cols_annot <- c("#FD8D3C", "#41B6C4", "#7FCDBB", "#081D58", "#225EA8", "#E31A1C", "black")


################################# Cluster annotations
new.cluster.ids <- c("PER-1", "Ven. VSMC", "Fibroblast", "Ven. VSMC", "PER-1", "Art. VSMC", "Ven. VSMC", "Act. PER",
                     "Art. VSMC", "Art. VSMC", "PER-1", "Schwann", "Art. VSMC", "Ven. VSMC", "PER-2", "Fibroblast")
names(new.cluster.ids) <- levels(cspg4_mcam_day7PI)
cspg4_mcam_day7PI <- RenameIdents(cspg4_mcam_day7PI, new.cluster.ids)
cspg4_mcam_day7PI@active.ident <- factor(cspg4_mcam_day7PI@active.ident, levels = c("Fibroblast", "Ven. VSMC", "Art. VSMC", "PER-1", "PER-2",
                                                                                    "Act. PER",  "Schwann"))
cspg4_mcam_day7PI@meta.data$celltype <- Idents(cspg4_mcam_day7PI)


################################# Visualizations
# UMAP
UMAP_annot <- DimPlot(cspg4_mcam_day7PI, reduction = "umap", label = F, pt.size = 0.5, cols = cols_annot) + theme(text = element_text(face = "bold", size=16), axis.text=element_text(size=6), axis.title.x=element_blank(),
                                                                                                                  axis.title.y=element_blank(), legend.position = "left") + NoAxes()

# Violin plots
VlnPlot(object = cspg4_mcam_day7PI, features = "GENE", same.y.lims = TRUE, sort = F, cols = cols_annot, pt.size = 0.1) +
  theme(text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(vjust=-1, face = "italic"))
