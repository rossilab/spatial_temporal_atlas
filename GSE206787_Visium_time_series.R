################################# Loading R packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(CellChat)
library(ComplexHeatmap)

################################# Loading spot-count matrices and associated images
data_dir <- "INSERT CORRESPONDING DIRECTORY"
list.files(data_dir)
Heart_D0 <- Load10X_Spatial(data.dir = data_dir, assay="Spatial", filename = "filtered_feature_bc_matrix.h5", slice = "Heart_D0", filter.matrix = T)
Heart_D0@meta.data$group <- "D0"
Heart_D0[["percent.mt_Spatial"]] <- PercentageFeatureSet(Heart_D0, pattern = "^mt-")
VlnPlot(Heart_D0, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2, combine = T, pt.size = 0.05, cols = c("grey"))

data_dir <- "INSERT CORRESPONDING DIRECTORY"
list.files(data_dir)
Heart_D3<- Load10X_Spatial(data.dir = data_dir, assay="Spatial", filename = "filtered_feature_bc_matrix.h5", slice = "Heart_D3",
                           filter.matrix = T)
Heart_D3@meta.data$group <- "D3"
Heart_D3[["percent.mt_Spatial"]] <- PercentageFeatureSet(Heart_D3, pattern = "^mt-")
VlnPlot(Heart_D3, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2, combine = T, pt.size = 0.05, cols = c("grey"))

data_dir <- "INSERT CORRESPONDING DIRECTORY"
list.files(data_dir)
Heart_D7<- Load10X_Spatial(data.dir = data_dir, assay="Spatial", filename = "filtered_feature_bc_matrix.h5", slice = "Heart_D7",
                           filter.matrix = T)
Heart_D7@meta.data$group <- "D7"
Heart_D7[["percent.mt_Spatial"]] <- PercentageFeatureSet(Heart_D7, pattern = "^mt-")
VlnPlot(Heart_D7, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2, combine = T, pt.size = 0.05, cols = c("grey"))

data_dir <- "INSERT CORRESPONDING DIRECTORY"
list.files(data_dir)
Heart_D14<- Load10X_Spatial(data.dir = data_dir, assay="Spatial", filename = "filtered_feature_bc_matrix.h5", slice = "Heart_D14",
                            filter.matrix = T)
Heart_D14@meta.data$group <- "D14"
Heart_D14[["percent.mt_Spatial"]] <- PercentageFeatureSet(Heart_D14, pattern = "^mt-")
VlnPlot(Heart_D14, features = c("nFeature_Spatial", "nCount_Spatial"), ncol = 2, combine = T, pt.size = 0.05, cols = c("grey"))


################################# Data integration
data_list = list(D0 = Heart_D0, D3 = Heart_D3, D7 = Heart_D7, D14 = Heart_D14)
data_list = lapply(data_list, SCTransform, assay = "Spatial", return.only.var.genes = F)
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB
var_features = SelectIntegrationFeatures(data_list, nfeatures = 3000, verbose = FALSE)
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = var_features,
                                verbose = T)
anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT",
                                  verbose = T, anchor.features = var_features, dims = 1:30)
visium_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                                   verbose = T, dims = 1:30)

################################# PCA, clustering, and UMAP projection
visium_integrated <- RunPCA(visium_integrated, features = VariableFeatures(object = visium_integrated), npcs = 30)
# ElbowPlot(visium_integrated, 30)
visium_integrated <- FindNeighbors(visium_integrated, dims = 1:20)
visium_integrated <- FindClusters(visium_integrated, resolution = 0.4)
visium_integrated <- RunUMAP(visium_integrated, dims = 1:20)


################################# Modifying metadata
visium_integrated@images$Heart_D0@coordinates$tissue <- as.integer(visium_integrated@images$Heart_D0@coordinates$tissue)
visium_integrated@images$Heart_D0@coordinates$row <- as.integer(visium_integrated@images$Heart_D0@coordinates$row)
visium_integrated@images$Heart_D0@coordinates$col <- as.integer(visium_integrated@images$Heart_D0@coordinates$col)
visium_integrated@images$Heart_D0@coordinates$imagerow <- as.integer(visium_integrated@images$Heart_D0@coordinates$imagerow)
visium_integrated@images$Heart_D0@coordinates$imagecol <- as.integer(visium_integrated@images$Heart_D0@coordinates$imagecol)

visium_integrated@images$Heart_D3@coordinates$tissue <- as.integer(visium_integrated@images$Heart_D3@coordinates$tissue)
visium_integrated@images$Heart_D3@coordinates$row <- as.integer(visium_integrated@images$Heart_D3@coordinates$row)
visium_integrated@images$Heart_D3@coordinates$col <- as.integer(visium_integrated@images$Heart_D3@coordinates$col)
visium_integrated@images$Heart_D3@coordinates$imagerow <- as.integer(visium_integrated@images$Heart_D3@coordinates$imagerow)
visium_integrated@images$Heart_D3@coordinates$imagecol <- as.integer(visium_integrated@images$Heart_D3@coordinates$imagecol)

visium_integrated@images$Heart_D7@coordinates$tissue <- as.integer(visium_integrated@images$Heart_D7@coordinates$tissue)
visium_integrated@images$Heart_D7@coordinates$row <- as.integer(visium_integrated@images$Heart_D7@coordinates$row)
visium_integrated@images$Heart_D7@coordinates$col <- as.integer(visium_integrated@images$Heart_D7@coordinates$col)
visium_integrated@images$Heart_D7@coordinates$imagerow <- as.integer(visium_integrated@images$Heart_D7@coordinates$imagerow)
visium_integrated@images$Heart_D7@coordinates$imagecol <- as.integer(visium_integrated@images$Heart_D7@coordinates$imagecol)

visium_integrated@images$Heart_D14@coordinates$tissue <- as.integer(visium_integrated@images$Heart_D14@coordinates$tissue)
visium_integrated@images$Heart_D14@coordinates$row <- as.integer(visium_integrated@images$Heart_D14@coordinates$row)
visium_integrated@images$Heart_D14@coordinates$col <- as.integer(visium_integrated@images$Heart_D14@coordinates$col)
visium_integrated@images$Heart_D14@coordinates$imagerow <- as.integer(visium_integrated@images$Heart_D14@coordinates$imagerow)
visium_integrated@images$Heart_D14@coordinates$imagecol <- as.integer(visium_integrated@images$Heart_D14@coordinates$imagecol)


################################# Saving processed Seurat object & aesthetics
visium_integrated@active.assay <- "SCT"
visium_integrated@meta.data$group <- factor(visium_integrated@meta.data$group, levels = c("D0", "D3", "D7", "D14"))
saveRDS(visium_integrated, file = "integrated_visium.rds")

cols_clusters <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#6A3D9A")
cols_region <- c("#1F78B4", "#FB9A99", "#E31A1C")
cols_timepoints <- c("#80CDC1", "#CB181D", "#1B7837", "#542788")
colors <- brewer.pal(n=9,name="Reds")


################################# Cluster annotations
new.cluster.ids <- c("Remote", "Remote", "Remote", "Infarct-1", "Infarct-2", "Remote", "Remote")
names(new.cluster.ids) <- levels(visium_integrated)
visium_integrated <- RenameIdents(visium_integrated, new.cluster.ids)
visium_integrated@active.ident <- factor(visium_integrated@active.ident, levels = c("Remote", "Infarct-1", "Infarct-2"))
visium_integrated@meta.data$annot <- as.character(visium_integrated@active.ident)
# Manually setting all regions of D0 to remote
visium_integrated@meta.data$annot[which(visium_integrated@meta.data$group == "D0")] <- "Remote"
Idents(visium_integrated) <- visium_integrated@meta.data$annot
# Binarize to remote or infarct regions
visium_integrated@meta.data$region <- as.character(visium_integrated@meta.data$annot)
visium_integrated@meta.data$region[which(visium_integrated@meta.data$region %in% c("Infarct-1", "Infarct-2"))] <- "Infarct"


################################# Visualization
# Spatial dim plot
Idents(visium_integrated) <- visium_integrated@meta.data$group
Heart_D0 <- subset(visium_integrated, idents = c("D0"))
Heart_D3 <- subset(visium_integrated, idents = c("D3"))
Heart_D7 <- subset(visium_integrated, idents = c("D7"))
Heart_D14 <- subset(visium_integrated, idents = c("D14"))

Idents(Heart_D0) <- Heart_D0@meta.data$seurat_clusters
Idents(Heart_D3) <- Heart_D3@meta.data$seurat_clusters
Idents(Heart_D7) <- Heart_D7@meta.data$seurat_clusters
Idents(Heart_D14) <- Heart_D14@meta.data$seurat_clusters

p <- SpatialDimPlot(Heart_D0, label = F, images = "Heart_D0", cols = cols_clusters, pt.size.factor = 1.65, crop = TRUE, alpha = c(0.1,1.5), stroke=1, repel = T) +
  theme(text = element_text(size = 12)) + NoLegend()
p <- SpatialDimPlot(Heart_D3, label = F, images = "Heart_D3", cols = cols_clusters, pt.size.factor = 1.5, crop = TRUE, alpha = c(0.1,1.5), stroke=1, repel = T) +
  theme(text = element_text(size = 12)) + NoLegend()
p <- SpatialDimPlot(Heart_D7, label = F, images = "Heart_D7", cols = cols_clusters, pt.size.factor = 1.7, crop = TRUE, alpha = c(0.1,1.5), stroke=1, repel = T) +
  theme(text = element_text(size = 12)) + NoLegend()
p <- SpatialDimPlot(Heart_D14, label = F, images = "Heart_D14", cols = cols_clusters, pt.size.factor = 1.4, crop = TRUE, alpha = c(0.1,1.5), stroke=1, repel = T) +
  theme(text = element_text(size = 12)) + NoLegend()

# Dot plot
visium_integrated@active.assay <- "Spatial"
visium_integrated <- NormalizeData(visium_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
visium_integrated@meta.data$time_region <- visium_integrated@meta.data$annot
visium_integrated@meta.data$time_region[which(visium_integrated@meta.data$group == "D0")] <- "SS"
Idents(visium_integrated) <- factor(visium_integrated@meta.data$time_region,
                                    levels = c("Infarct-2", "Infarct-1", "Remote", "SS"))
visium_integrated@active.assay <- "Spatial"
d <- DotPlot(visium_integrated, features = "GENES", dot.scale = 10, scale.by = "size") + xlab("") + ylab("")+
  RotatedAxis() + theme_bw() + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size = 16)) +
  geom_vline(xintercept=c(4.5,8.5,12.5), linetype = "dashed", colour = "grey60") + scale_color_gradient2(low = colors[1], mid = colors[5], high = colors[9])

################################# Spatiotemporal analysis
# Identify temporally regulated genes in the infarct zone
visium_integrated@meta.data$region[which(visium_integrated@meta.data$group == "D0")] <- "Infarct"
Idents(visium_integrated) <- visium_integrated@meta.data$region
subsetData <- subset(visium_integrated, idents = "Infarct")
Idents(subsetData) <- subsetData@meta.data$group
subsetData@active.assay <- "SCT"
clusterMarkers_D3 <- FindMarkers(subsetData, ident.1 = c("D3"), ident.2 = c("D0"), min.pct = 0.25, only.pos = F)
clusterMarkers_D3 <- clusterMarkers_D3[which(clusterMarkers_D3$p_val_adj < 0.05),]
clusterMarkers_D7 <- FindMarkers(subsetData, ident.1 = c("D7"), ident.2 = c("D0"), min.pct = 0.25, only.pos = F)
clusterMarkers_D7 <- clusterMarkers_D7[which(clusterMarkers_D7$p_val_adj < 0.05),]
clusterMarkers_D14 <- FindMarkers(subsetData, ident.1 = c("D14"), ident.2 = c("D0"), min.pct = 0.25, only.pos = F)
clusterMarkers_D14 <- clusterMarkers_D14[which(clusterMarkers_D14$p_val_adj < 0.05),]
DEGs_D3 <- rownames(clusterMarkers_D3)
DEGs_D7 <- rownames(clusterMarkers_D7)
DEGs_D14 <- rownames(clusterMarkers_D14)
infarct_all_markers <- unique(c(DEGs_D3, DEGs_D7, DEGs_D14))

# Identify temporally regulated genes in the remote zone
visium_integrated@meta.data$region[which(visium_integrated@meta.data$group == "D0")] <- "Remote"
Idents(visium_integrated) <- visium_integrated@meta.data$region
subsetData <- subset(visium_integrated, idents = "Remote")
Idents(subsetData) <- subsetData@meta.data$group
subsetData@active.assay <- "SCT"
clusterMarkers_D3 <- FindMarkers(subsetData, ident.1 = c("D3"), ident.2 = c("D0"), min.pct = 0.25, only.pos = F)
clusterMarkers_D3 <- clusterMarkers_D3[which(clusterMarkers_D3$p_val_adj < 0.05),]
clusterMarkers_D7 <- FindMarkers(subsetData, ident.1 = c("D7"), ident.2 = c("D0"), min.pct = 0.25, only.pos = F)
clusterMarkers_D7 <- clusterMarkers_D7[which(clusterMarkers_D7$p_val_adj < 0.05),]
clusterMarkers_D14 <- FindMarkers(subsetData, ident.1 = c("D14"), ident.2 = c("D0"), min.pct = 0.25, only.pos = F)
clusterMarkers_D14 <- clusterMarkers_D14[which(clusterMarkers_D14$p_val_adj < 0.05),]
DEGs_D3 <- rownames(clusterMarkers_D3)
DEGs_D7 <- rownames(clusterMarkers_D7)
DEGs_D14 <- rownames(clusterMarkers_D14)
remote_all_markers <- unique(c(DEGs_D3, DEGs_D7, DEGs_D14))

# K-means clustering of spatio-temporally regulated programs
visium_integrated@meta.data$region[which(visium_integrated@meta.data$group == "D0")] <- "Remote"
visium_integrated@active.assay <- "Spatial"
visium_integrated <- NormalizeData(visium_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
visium_integrated@meta.data$time_region <- paste(visium_integrated@meta.data$group, visium_integrated@meta.data$region)
visium_integrated@meta.data$time_region[which(visium_integrated@meta.data$time_region == "D0 Remote")] <- "SS"
visium_integrated@meta.data$time_region <- factor(visium_integrated@meta.data$time_region, levels = c("SS", "D3 Remote", "D3 Infarct", "D7 Remote", "D7 Infarct", "D14 Remote", "D14 Infarct"))
Idents(visium_integrated) <- visium_integrated@meta.data$time_region
cluster.averages <- AverageExpression(object = visium_integrated, return.seurat = TRUE, verbose = T, assays = "Spatial")
mat <- as.matrix(cluster.averages@assays$Spatial@data)
cluster.averages@meta.data$annot <- c("SS", "D3 Remote", "D3 Infarct", "D7 Remote", "D7 Infarct", "D14 Remote", "D14 Infarct")
colData <- data.frame(cluster.averages@meta.data[,4])
rownames(colData) <- rownames(cluster.averages@meta.data)
colnames(colData) <- c("annot")
cols_annot <- list(annot = c("SS" = "#80CDC1", "D3 Remote" = "#CB181D", "D3 Infarct" = "#CB181D", "D7 Remote" = "#1B7837",
                             "D7 Infarct" = "#1B7837", "D14 Remote" = "#542788", "D14 Infarct" = "#542788"))
all_markers <- unique(c(remote_all_markers, infarct_all_markers))
df <- mat[all_markers,]
elbow <- fviz_nbclust(df, kmeans, method = "wss")
elbow_df <- elbow$data
g <- ggplot(data=elbow_df, aes(x=clusters, y=y, group =1)) +
  geom_line(size=1.1) +
  geom_point(size = 3, color = "grey") +
  xlab("Number of clusters") +
  ylab("") +
  theme_bw() +
  theme(text=element_text(size=20), legend.position = "none")
breaksList = seq(-2, 2, length.out = 255)
colors <- colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(255)
p <- pheatmap(df,
              main = "",
              treeheight_row = 30, treeheight_col = 15,
              color = colors,
              angle_col = 45,
              show_colnames = F,
              show_rownames  = F,
              clustering_method = "complete",
              scale = "row",
              annotation_col = colData,
              cluster_cols = F,
              cluster_rows = T,
              annotation_names_col = F,
              annotation_colors = cols_annot,
              breaks = breaksList,
              fontsize = 16,
              kmeans_k = 4,
              annotation_legend = F,
              legend = T)
clus_assign <- p$kmeans$cluster
clus_1 <- names(clus_assign)[which(clus_assign == "1")]
clus_2 <- names(clus_assign)[which(clus_assign == "4")]
clus_3 <- names(clus_assign)[which(clus_assign == "2")]
clus_4 <- names(clus_assign)[which(clus_assign == "3")]

# heatmap of selected genes across regions
vascular_final <- c("Pdgfb",
                    "Hif1a",
                    "Flt1", "Kdr", "Nrp1", "Sulf1", "Sulf2", "Vegfa", "Vegfb", "Vegfc", "Vegfd",
                    "Hes1", "Jag1", "Notch1", "Notch3",
                    "Pxn", "Robo1","Robo4", "Slit3",
                    "Mob1b", "Taz", "Tead1", "Yap1",
                    "Angpt1", "Angptl2", "Angptl4", "Tek",
                    "Thbs1", "Thbs2", "Thbs3", "Thbs4",
                    "Aggf1", "Ang", "Anxa1", "Anxa2", "Apln", "Emc10", "Ets2", "Fgf1", "Ptn",
                    "Agt", "Nppa")
ECM_final <- c("Pdgfa", "Col1a1", "Fn1", "Col15a1", "Lamb1", "Eln", "Comp", "Lum",
               "Adam8", "Adamts1", "Adamts2", "Bmp1", 'Mmp2', "Mmp10", "Timp1",
               "Lox", "Loxl2", 'Nid1',"Plod3",
               "Ccn2", "Ccn5", "Postn", "Sparc","Spp1",
               "Acvrl1", "Smad1", "Smad5", "Smad3", "Tgfb1", 'Tgfb2', "Tgfb3", "Tgfbr2",
               "Apc", "Ctnnb1", "Fzd1", "Fzd2", "Gsk3b","Lrp5", "Lrp6", "Sfrp1", "Sfrp2", "Wnt5a")
inflam_final <- c("Ager", "Hmgb1", "Myd88", "Tirap","Tlr9",
                  "Nfkb1", "Nfkb2", "Rela", "Relb", 
                  "Ifitm3", "Irf3", "Irf5", "Irf7", 
                  "C1qa", "C1qb", "C1ra", "C3", "C3ar1",
                  "Csf1", "Csf1r", "Il3ra", "Il17ra", 
                  "Tnf", "Tnfrsf1a", "Tnfrsf1b", "Tradd",
                  "Il33", "Il1rn", "Il4ra", "Il10ra",
                  "Ccl3", "Ccl5", "Ccl7", "Cxcl1", "Cxcl10", "Cxcl12", "Icam1", "Selp",
                  "Gpx1", "Prdx4", "Nos2", "Selenos", "Sod3", 
                  "Apoe", "Clta", "Cltb", "Ctss", "Ctsc", "Lrp1", "Trem2")
breaksList = seq(-2, 2, length.out = 255)
colors <- colorRampPalette(brewer.pal(n=7,name="YlGnBu"))(255)
selected_genes <- mat[inflam_final,] # change gene list accordingly
p <- pheatmap(selected_genes,
              main = NA,
              treeheight_row = 30, treeheight_col = 15,
              color = colors,
              angle_col = 45,
              show_colnames = F,
              show_rownames  = T,
              clustering_method = "complete",
              scale = "row",
              annotation_col = colData,
              cluster_cols = F,
              cluster_rows = F,
              annotation_names_col = F,
              annotation_colors = cols_annot,
              breaks = breaksList,
              fontsize = 24,
              annotation_legend = F,
              legend = F)


################################# Mapping of cells from CD146-sorted cells onto Visium data (D3)
# SCTransform CD146 D3 data
mcam_day3PI_sct <- readRDS(file = file.choose())
mcam_day3PI_sct[["percent.mt"]] <- PercentageFeatureSet(mcam_day3PI_sct, pattern = "^mt-")

# VlnPlot(mcam_day3PI_sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mcam_day3PI_sct <- subset(mcam_day3PI_sct, subset = nFeature_RNA > 200 & percent.mt < 15 & nCount_RNA < 30000)
mcam_day3PI_sct <- SCTransform(mcam_day3PI_sct, vars.to.regress = "percent.mt", verbose = T)
mcam_day3PI_sct <- RunPCA(mcam_day3PI_sct, assay = "SCT", verbose = FALSE)
mcam_day3PI_sct <- FindNeighbors(mcam_day3PI_sct, reduction = "pca", dims = 1:20)
mcam_day3PI_sct <- FindClusters(mcam_day3PI_sct, verbose = T, resolution = 0.6)
mcam_day3PI_sct <- RunUMAP(mcam_day3PI_sct, reduction = "pca", dims = 1:20, verbose = T)

# Grabbing metadata from timeseries object
mcam_anchored <- readRDS(file = file.choose())
Idents(mcam_anchored) <- mcam_anchored@meta.data$group
mcam_day3PI <- subset(mcam_anchored, idents = "3d PI")
Idents(mcam_day3PI) <- mcam_day3PI@meta.data$celltype
mcam_day3PI_sct@meta.data$celltype <- mcam_day3PI@meta.data$celltype
Idents(mcam_day3PI_sct) <- mcam_day3PI_sct@meta.data$celltype
saveRDS(mcam_day3PI_sct, file = "seurat_D3_sct.rds")

# Label transfer and prediction of localization
anchors <- FindTransferAnchors(reference = mcam_day3PI, query = Heart_D0, normalization.method = "SCT", dims = 1:20,
                               reference.assay = "SCT", query.assay = "SCT")
mcam_day3PI@meta.data$celltype <- as.character(Idents(mcam_day3PI))
predictions.assay <- TransferData(anchorset = anchors, refdata = mcam_day3PI@meta.data$celltype, prediction.assay = TRUE, 
                                  weight.reduction = Heart_D0[["pca"]], dims = 1:20)
Heart_D0[["predictions"]] <- predictions.assay
Heart_D0@active.assay <- "predictions"
SpatialFeaturePlot(Spatial_Heart, features = "ANNOTATION", pt.size.factor = 2, crop = TRUE, stroke=1, max.cutoff = "q95", min.cutoff = 0) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme(text = element_text(size = 16), legend.key.size = unit(0.9, "cm"), legend.position = "right", legend.title = element_blank())


################################# CellChat - analysis
# Prepare input data for CellChat analysis
Idents(Heart_D3) <- Heart_D3@meta.data$annot
data_input_D3 <- GetAssayData(Heart_D3, slot = "data", assay = "SCT")
meta_D3 = data.frame(labels = Idents(Heart_D3), row.names = names(Idents(Heart_D3)))
Heart_D3@images$Heart_D0 <- NULL
Heart_D3@images$Heart_D7 <- NULL
Heart_D3@images$Heart_D14 <- NULL

Idents(Heart_D7) <- Heart_D7@meta.data$annot
data_input_D7 <- GetAssayData(Heart_D7, slot = "data", assay = "SCT")
meta_D7 = data.frame(labels = Idents(Heart_D7), row.names = names(Idents(Heart_D7)))
Heart_D7@images$Heart_D0 <- NULL
Heart_D7@images$Heart_D3 <- NULL
Heart_D7@images$Heart_D14 <- NULL

Idents(Heart_D14) <- Heart_D14@meta.data$annot
data_input_D14 <- GetAssayData(Heart_D14, slot = "data", assay = "SCT")
meta_D14 = data.frame(labels = Idents(Heart_D14), row.names = names(Idents(Heart_D14)))
Heart_D14@images$Heart_D0 <- NULL
Heart_D14@images$Heart_D3 <- NULL
Heart_D14@images$Heart_D7 <- NULL

# Load spatial imaging information
spatial_locs_D3 = GetTissueCoordinates(Heart_D3, scale = NULL, cols = c("imagerow", "imagecol")) 
scale_factors_D3 = jsonlite::fromJSON(txt = file.path("C:/Users/David/Desktop/Henry/Visium/Data/D3_Visium/spatial", 'scalefactors_json.json'))
scale_factors_D3 = list(spot.diameter = 55, spot = scale_factors_D3$spot_diameter_fullres, # these two information are required
                        fiducial = scale_factors_D3$fiducial_diameter_fullres, hires = scale_factors_D3$tissue_hires_scalef, lowres = scale_factors_D3$tissue_lowres_scalef
)

spatial_locs_D7 = GetTissueCoordinates(Heart_D7, scale = NULL, cols = c("imagerow", "imagecol")) 
scale_factors_D7 = jsonlite::fromJSON(txt = file.path("C:/Users/David/Desktop/Henry/Visium/Data/D7_Visium/spatial", 'scalefactors_json.json'))
scale_factors_D7 = list(spot.diameter = 55, spot = scale_factors_D7$spot_diameter_fullres, # these two information are required
                        fiducial = scale_factors_D7$fiducial_diameter_fullres, hires = scale_factors_D7$tissue_hires_scalef, lowres = scale_factors_D7$tissue_lowres_scalef
)

spatial_locs_D14 = GetTissueCoordinates(Heart_D14, scale = NULL, cols = c("imagerow", "imagecol")) 
scale_factors_D14 = jsonlite::fromJSON(txt = file.path("C:/Users/David/Desktop/Henry/Visium/Data/D14_Visium/spatial", 'scalefactors_json.json'))
scale_factors_D14 = list(spot.diameter = 55, spot = scale_factors_D14$spot_diameter_fullres, # these two information are required
                         fiducial = scale_factors_D14$fiducial_diameter_fullres, hires = scale_factors_D14$tissue_hires_scalef, lowres = scale_factors_D14$tissue_lowres_scalef
)

# Create CellChat objects
cellchat_D3 <- createCellChat(object = data_input_D3, meta = meta_D3, group.by = "labels",
                              datatype = "spatial", coordinates = spatial_locs_D3, scale.factors = scale_factors_D3)

cellchat_D7 <- createCellChat(object = data_input_D7, meta = meta_D7, group.by = "labels",
                              datatype = "spatial", coordinates = spatial_locs_D7, scale.factors = scale_factors_D7)

cellchat_D14 <- createCellChat(object = data_input_D14, meta = meta_D14, group.by = "labels",
                               datatype = "spatial", coordinates = spatial_locs_D14, scale.factors = scale_factors_D14)

# Select CellChat database and category of interactions
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- CellChatDB

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
secreted_signaling <- unique(CellChatDB.use$interaction$pathway_name)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use cell-cell contact
cell_contact <- unique(CellChatDB.use$interaction$pathway_name)
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use ECM receptor
ECM_receptor <- unique(CellChatDB.use$interaction$pathway_name)

cellchat_D3@DB <- CellChatDB.use
cellchat_D7@DB <- CellChatDB.use
cellchat_D14@DB <- CellChatDB.use

cellchat_D3 <- subsetData(cellchat_D3) # This step is necessary even if using the whole database
cellchat_D7 <- subsetData(cellchat_D7) # This step is necessary even if using the whole database
cellchat_D14 <- subsetData(cellchat_D14) # This step is necessary even if using the whole database

# Compute over-expressed ligands or receptors in each annotated group
cellchat_D3 <- identifyOverExpressedGenes(cellchat_D3, thresh.p = 0.1)
cellchat_D3 <- identifyOverExpressedInteractions(cellchat_D3)
cellchat_D3 <- projectData(cellchat_D3, PPI.mouse)

cellchat_D7 <- identifyOverExpressedGenes(cellchat_D7, thresh.p = 0.1)
cellchat_D7 <- identifyOverExpressedInteractions(cellchat_D7)
cellchat_D7 <- projectData(cellchat_D7, PPI.mouse)

cellchat_D14 <- identifyOverExpressedGenes(cellchat_D14, thresh.p = 0.1)
cellchat_D14 <- identifyOverExpressedInteractions(cellchat_D14)
cellchat_D14 <- projectData(cellchat_D14, PPI.mouse)

# Compute communication probability and infer cell communication network
cellchat_D3 <- computeCommunProb(cellchat_D3, type = "truncatedMean", trim = 0.1, 
                                 distance.use = T, interaction.length = 275, scale.distance = 0.01,
                                 raw.use = T, population.size = T)
cellchat_D7 <- computeCommunProb(cellchat_D7, type = "truncatedMean", trim = 0.1, 
                                 distance.use = T, interaction.length = 275, scale.distance = 0.01,
                                 raw.use = T, population.size = T)
cellchat_D14 <- computeCommunProb(cellchat_D14, type = "truncatedMean", trim = 0.1, 
                                  distance.use = T, interaction.length = 275, scale.distance = 0.01,
                                  raw.use = T, population.size = T)

# Remove cell-cell communication if an annotated group has less than 10 cells
cellchat_D3 <- filterCommunication(cellchat_D3, min.cells = 10)
cellchat_D7 <- filterCommunication(cellchat_D7, min.cells = 10)
cellchat_D14 <- filterCommunication(cellchat_D14, min.cells = 10)

# Compute pathway-specific communication probabilities by summarizing all ligand/receptor interactions associated with a pathway
cellchat_D3 <- computeCommunProbPathway(cellchat_D3)
cellchat_D7 <- computeCommunProbPathway(cellchat_D7)
cellchat_D14 <- computeCommunProbPathway(cellchat_D14)

# Tabulate communication strength and links
cellchat_D3 <- aggregateNet(cellchat_D3)
cellchat_D7 <- aggregateNet(cellchat_D7)
cellchat_D14 <- aggregateNet(cellchat_D14)

# Compute the network centrality scores
cellchat_D3 <- netAnalysis_computeCentrality(cellchat_D3, slot.name = "netP") 
cellchat_D7 <- netAnalysis_computeCentrality(cellchat_D7, slot.name = "netP") 
cellchat_D14 <- netAnalysis_computeCentrality(cellchat_D14, slot.name = "netP") 

# Saving processed Cellchat object & aesthetics
cellchat_list <- list(D3 = cellchat_D3, D7 = cellchat_D7, D14 = cellchat_D14)
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
saveRDS(cellchat_list, "cellchat_list_allDB_275um.rds")
saveRDS(cellchat_merged, "cellchat_merged_allDB_275um.rds")
cols_timepoints <- c("#CB181D", "#1B7837", "#542788")
cols_region <- c("#1F78B4", "#FB9A99", "#E31A1C")

# Identify active signaling pathways in at least one time point
pathway_union <- Reduce(union, list(cellchat_list[[1]]@netP$pathways, cellchat_list[[2]]@netP$pathways, cellchat_list[[3]]@netP$pathways))

# Selected pathways involving secreted ligands
ECM_paths_heatmap <- c("ACTIVIN", "BMP", "GDF", "PERIOSTIN", "SPP1", "TENASCIN", "TGFb", "WNT", "ncWNT")
vasc_paths_heatmap <- c("ANGPT", "ANGPTL", "EGF", "FGF", "IGF", "PTN", "THBS", "VEGF")
immune_paths_heatmap <- c("CCL", "COMPLEMENT", "CSF", "GALECTIN", "IL1", "IL6", "IL16", "MIF", "OSM", "TNF", "VISFATIN")


################################# CellChat - visualizations
# Display # of interactions across regions
n <- netVisual_circle(cellchat_list[["INSERT TIMEPOINT"]]@net$count, weight.scale = T,
                      label.edge = T, edge.weight.max = weight.max[2],
                      edge.width.max = 12, edge.label.cex = 1.8, color.use = cols_region,
                      arrow.size = 0.5, arrow.width = 1, vertex.label.cex = NA,
                      margin = 0.0000000001)
n[[1]][[34]][[2]][[3]] <- NA

# Graphing of information flow (i.e. sum of communication probability among cell groups)
g <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = F, comparison = c(1,2,3),
             color.use = cols_timepoints, signaling = pathway_union, measure = "count", do.flip = F,
             font.size = 16) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend() # num of interactions
g <- rankNet(cellchat_merged, mode = "comparison", stacked = T, do.stat = F, comparison = c(1,2,3),
             color.use = cols_timepoints, signaling = pathway_union, measure = "weight", do.flip = T,
             font.size = 16) +
  theme(axis.title.x = element_blank()) + NoLegend() # weight of interactions

# Graphing of incoming or outgoing signals
all_secreted <- c(immune_paths_heatmap, ECM_paths_heatmap, vasc_paths_heatmap)
ht1 = netAnalysis_signalingRole_heatmap(cellchat_list[["INSERT TIMEPOINT"]], pattern = "outgoing", signaling = all_secreted,
                                        title = names(cellchat_list)[1], width = 8, height = 22,
                                        color.use = cols_region, font.size = 12, color.heatmap = "Reds")
ht1@column_names_param$show <- F
ht1@column_title <- NULL
ht1 = netAnalysis_signalingRole_heatmap(cellchat_list[["INSERT TIMEPOINT"]], pattern = "incoming", signaling = all_secreted,
                                        title = names(cellchat_list)[1], width = 8, height = 22,
                                        color.use = cols_region, font.size = 12, color.heatmap = "GnBu")
ht1@column_names_param$show <- F
ht1@column_title <- NULL


################################# CellChat - Graphing of dot plots corresponding to the cumulative activity of each ligand across regions
ECM_paths_secreted <- c("TGFb", "BMP", "GDF", "ncWNT", "PERIOSTIN", "TENASCIN")
ECM_paths_secreted_problem <- c("ACTIVIN", "WNT", "SPP1")
vasc_paths_secreted <- c("VEGF", "ANGPT", "PTN", "FGF", "EGF", "IGF", "ANGPTL", "THBS")
immune_paths_secreted <- c("CCL", "OSM", "IL1", "CSF", "VISFATIN", "COMPLEMENT", "GALECTIN")
immune_paths_secreted_problem <- c("IL6", "MIF", "TNF", "IL16")
structural <- c("COLLAGEN", "FN1", "LAMININ", "HSPG")
contact <- c("CD34", "EPHB", "ICAM", "NOTCH", "SEMA4", "SEMA6", "SEMA7", "VCAM")
structural_contact <- c(structural, contact)
problem_contact <- c("EPHA", "MHC-II", "SELE")

dt <- netVisual_bubble(cellchat_merged, sources.use = 1:3, targets.use = 1:3,
                       signaling = structural_contact, comparison = c(1,2,3), angle.x = 45, return.data = T)

# Data wrangling/processing for plotting
# Identify significant LR interactions that were missed by netVisual_bubble
# If a signaling pathway (all LR pairs) is missing in one of the timepoints, netVisual_bubble
# is unable to capture the pathway
problem_df <- data.frame()
for (signal in ECM_paths_secreted_problem) {
  i <- 1
  while (i <= 3) {
    LR <- try(netVisual_bubble(cellchat_list[[i]], sources.use = 1:3, targets.use = 1:3,
                               signaling = signal,  angle.x = 45, return.data = T))
    if("try-error" %in% class(LR)){
      print(paste(signal, "not found in", i, sep = " "))
      i <- i + 1
      break
    }
    else{
      LR <- LR$communication
      LR$group.names <- LR$source.target
      if(i == 1){
        LR$dataset <- "D3"
      }
      if(i == 2){
        LR$dataset <- "D7"
      }
      if(i == 3){
        LR$dataset <- "D14"
      }
      LR$source.target <- paste(LR$source.target, LR$dataset, sep = " ")
      problem_df <- rbind(problem_df, LR)
      i <- i + 1
    }
  }
}
for (signal in vasc_paths_secreted_problem) {
  i <- 1
  while (i <= 3) {
    LR <- try(netVisual_bubble(cellchat_list[[i]], sources.use = 1:3, targets.use = 1:3,
                               signaling = signal,  angle.x = 45, return.data = T))
    if("try-error" %in% class(LR)){
      print(paste(signal, "not found in", i, sep = " "))
      i <- i + 1
      break
    }
    else{
      LR <- LR$communication
      LR$group.names <- LR$source.target
      if(i == 1){
        LR$dataset <- "D3"
      }
      if(i == 2){
        LR$dataset <- "D7"
      }
      if(i == 3){
        LR$dataset <- "D14"
      }
      LR$source.target <- paste(LR$source.target, LR$dataset, sep = " ")
      problem_df <- rbind(problem_df, LR)
      i <- i + 1
    }
  }
}
for (signal in immune_paths_secreted_problem) {
  i <- 1
  while (i <= 3) {
    LR <- try(netVisual_bubble(cellchat_list[[i]], sources.use = 1:3, targets.use = 1:3,
                               signaling = signal,  angle.x = 45, return.data = T))
    if("try-error" %in% class(LR)){
      print(paste(signal, "not found in", i, sep = " "))
      i <- i + 1
      break
    }
    else{
      LR <- LR$communication
      LR$group.names <- LR$source.target
      if(i == 1){
        LR$dataset <- "D3"
      }
      if(i == 2){
        LR$dataset <- "D7"
      }
      if(i == 3){
        LR$dataset <- "D14"
      }
      LR$source.target <- paste(LR$source.target, LR$dataset, sep = " ")
      problem_df <- rbind(problem_df, LR)
      i <- i + 1
    }
  }
}
for (signal in problem_contact) {
  i <- 1
  while (i <= 3) {
    LR <- try(netVisual_bubble(cellchat_list[[i]], sources.use = 1:3, targets.use = 1:3,
                               signaling = signal,  angle.x = 45, return.data = T))
    if("try-error" %in% class(LR)){
      print(paste(signal, "not found in", i, sep = " "))
      i <- i + 1
      break
    }
    else{
      LR <- LR$communication
      LR$group.names <- LR$source.target
      if(i == 1){
        LR$dataset <- "D3"
      }
      if(i == 2){
        LR$dataset <- "D7"
      }
      if(i == 3){
        LR$dataset <- "D14"
      }
      LR$source.target <- paste(LR$source.target, LR$dataset, sep = " ")
      problem_df <- rbind(problem_df, LR)
      i <- i + 1
    }
  }
}

# Append "missing" LR interactions to "identified" LR interactions
pre_df <- dt$communication
df <- rbind(pre_df, problem_df) # Run if there are "missing" LR interactions
df$ligand[which(df$pathway_name == "ACTIVIN")] <- "Inhbb"

# [Optional] subset for "autocrine" signals
same_region <- c("Infarct-2 -> Infarct-2", "Infarct-1 -> Infarct-1", "Remote -> Remote")
same_region_signal <- df[which(df$group.names %in% same_region),]
df <- same_region_signal

# Collapsing LR interactions based on ligands
# Pval is ditched as every interaction is already significant
# The number of unique, significant interactions is counted instead
# Communication probability is summed across interactions with the same ligand
unique_ligands <- as.character(na.omit(unique(df$ligand)))
source_targets <- unique(as.character(df$source.target))
new_df <- data.frame()
for (ligand in unique_ligands) {
  ligand_df <- df[which(df$ligand == ligand),]
  for (source_target in source_targets) {
    source_df <- ligand_df[which(ligand_df$source.target == source_target),]
    source_df <- source_df[which(source_df$pval > 1),]
    if(nrow(source_df) > 0){
      prob <- sum(source_df$prob)
      prob_orig <- sum(source_df$prob.original)
      source_df$prob <- prob
      source_df$prob.original <- prob_orig
      source_df$pval <- length(unique(source_df$receptor))
      new_df <- rbind(new_df, source_df)
    }
  }
}

# MISC. modifications
new_df$interaction_name_2 <- new_df$ligand
colnames(new_df)[6] <- "Count"

# Modify x-axis labels
new_df$source <- as.character(new_df$source)
new_df$target <- as.character(new_df$target)
i <- 1
while (i <= nrow(new_df)) {
  source <- new_df$source[i]
  source <- paste(substr(source,1,1),substr(source,nchar(source),nchar(source)), sep = "")
  target <- new_df$target[i]
  target <- paste(substr(target,1,1),substr(target,nchar(target),nchar(target)), sep = "")
  new_df$source[i] <- source
  new_df$target[i] <- target
  i <- i + 1
}
new_df$source[which(new_df$source == "Re")] <- "R"
new_df$target[which(new_df$target == "Re")] <- "R"
new_df$source <- factor(new_df$source, levels = c("R", "I1", "I2"))
new_df$target <- factor(new_df$target, levels = c("R", "I1", "I2"))
new_df$group.names2 <- paste(new_df$source, new_df$target, sep = " -> ")
new_df$group.names2 <- paste(new_df$group.names2, " (", new_df$dataset, ")", sep = "")
new_df$group.names2 <- factor(new_df$group.names2, levels = c("R -> R (D3)", "R -> R (D7)", "R -> R (D14)",
                                                              "R -> I1 (D3)", "R -> I1 (D7)", "R -> I1 (D14)",
                                                              "R -> I2 (D3)", "R -> I2 (D7)", "R -> I2 (D14)",
                                                              "I1 -> R (D3)", "I1 -> R (D7)", "I1 -> R (D14)",
                                                              "I1 -> I1 (D3)", "I1 -> I1 (D7)", "I1 -> I1 (D14)",
                                                              "I1 -> I2 (D3)", "I1 -> I2 (D7)", "I1 -> I2 (D14)",
                                                              "I2 -> R (D3)", "I2 -> R (D7)", "I2 -> R (D14)",
                                                              "I2 -> I1 (D3)", "I2 -> I1 (D7)", "I2 -> I1 (D14)",
                                                              "I2 -> I2 (D3)", "I2 -> I2 (D7)", "I2 -> I2 (D14)"))
vasc_ligands_keep <- c("Angpt1", "Angpt2", "Angptl2", "Angptl4", "Fgf1", "Fgf2",
                       "Igf1", "Igf2", "Ptn", "Thbs1", "Thbs2", "Thbs3", 'Thbs4',
                       "Vegfa", "Vegfb", "Vegfc", "Vegfd")
ecm_ligands_keep <- c("Bmp4", "Bmp6", "Gdf6", "Gdf11", "Gdf15", 'Inhbb', "Postn", "Spp1", "Tgfb1", "Tgfb2", "Tgfb3", "Tnc", "Tnxb",
                      "Wnt4", "Wnt5a", "Wnt5b", "Wnt9a", "Wnt11")
immune_ligands_keep <- c("C3", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Ccl11", "Ccl12", "Ccl22",
                         "Csf1", "Il1b", "Il11", "Il34", "Mif", "Tnf", "Osm")
structural_contact_keep <- c("Cd34", "Col1a1", "Col1a2", "Col4a1", "Col4a2", "Col6a1", "Col6a2", "Col9a2",
                             "Fn1", "Lama2", "Lama4", "Lama5", "Lamb1", "Lamb2", "Lamc1", "Hspg2", "Icam1", "Icam2",
                             "Dll1", "Dll4", "Jag1", "Jag2", "Efna2", "Efna4", "Efna5", "Efnb1", "Efnb2", "Sele")
ligands_to_keep <- structural_contact_keep
new_df <- new_df[which(new_df$ligand %in% ligands_to_keep),]
new_df$interaction_name_2 <- factor(new_df$interaction_name_2, levels = rev(ligands_to_keep))

df <- new_df
df$binarized <- 1

#Modification of plot parameters
color.use <- brewer.pal(n = 9, "Reds")
# Plot prototype
g <- ggplot(df, aes(x = group.names2, y = interaction_name_2, color = prob, size = binarized, label = Count)) +
  geom_point(pch = 16) +
  geom_text(hjust=-0.95, vjust=1.25, colour = "grey50", size = 3) +
  theme_linedraw() + theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom", drop=F)

g <- ggplot(df, aes(x = group.names2, y = interaction_name_2, color = prob, size = binarized, label = Count)) +
  geom_point(pch = 16) +
  geom_text(hjust=-0.95, vjust=1.25, colour = "grey50", size = 3) +
  theme_linedraw() + theme(panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom", drop=T)

# g <- ggplot(df, aes(x = group.names2, y = interaction_name_2, color = prob, label = Count)) +
#   geom_text(hjust=1, vjust=1) +
#   theme_linedraw() + theme(panel.grid.major = element_blank()) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   scale_x_discrete(position = "bottom", drop=F)

# Rescale color for communication probability
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
  g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                  breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
} else {
  g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), na.value = "white") +
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
}

# Change text size
g <- g + theme(axis.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16))

# Add grids
g <- g + geom_vline(xintercept=seq(1.5, length(levels(df$group.names2))),lwd=0.1,colour="grey60")
g <- g + geom_vline(xintercept=seq(1.5, 9),lwd=0.1,colour="grey60")

g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour="grey60")         

# Add dashed lines between sets of regional interactions
x_int <- seq(3.5,
             24.5,
             by = 3)
x_int <- seq(3.5,
             6.5,
             by = 3)
g <- g + geom_vline(xintercept=x_int, linetype="dashed", color = "grey50", size = 1)

# Color text based on timepoint
cols_timepoints <- c("#EF3B2C", "#5AAE61", "#8073AC")
g <- g + theme(axis.text.x = element_text(colour = rep(cols_timepoints,9)))

