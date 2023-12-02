################################# Loading R packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(SCENIC)
library(RcisTarget)
library(AUCell)
library(GENIE3)


################################# Remove doublets using hashing data
# Loading gene/HTO count matrix
hic1_undamaged <- Read10X(data.dir = "INSERT CORRESPONDING DIRECTORY") # Reading data into expression matrix
colnames(x = hic1_undamaged$`Gene Expression`) <- paste('hic1_undamaged', colnames(x = hic1_undamaged$`Gene Expression`), sep = '_')
colnames(x = hic1_undamaged$`Antibody Capture`) <- paste('hic1_undamaged', colnames(x = hic1_undamaged$`Antibody Capture`), sep = '_')
exprData <- CreateSeuratObject(counts = hic1_undamaged$`Gene Expression`, min.cells = 3, min.features = 200, project = "hic1_undamaged")
sharedCells <- intersect(colnames(exprData), colnames(hic1_undamaged$`Antibody Capture`))
HTOData <- hic1_undamaged$`Antibody Capture`[, sharedCells]
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
hic1_undamaged <- exprData
saveRDS(hic1_undamaged, "HTO_processed_object.rds")


################################# Quality control, Normalization, and Scaling
# hic1_undamaged <- readRDS(file.choose())
hic1_undamaged[["percent.mt"]] <- PercentageFeatureSet(hic1_undamaged, pattern = "^mt-")
# VlnPlot(hic1_undamaged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hic1_undamaged <- subset(hic1_undamaged, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA < 15000)
hic1_undamaged <- NormalizeData(hic1_undamaged, normalization.method = "LogNormalize", scale.factor = 10000)
hic1_undamaged <- FindVariableFeatures(hic1_undamaged, selection.method = "vst", nfeatures = 2000)
hic1_undamaged <- ScaleData(hic1_undamaged, vars.to.regress = "percent.mt", features = rownames(hic1_undamaged))


################################# PCA, clustering, and UMAP projection
hic1_undamaged <- RunPCA(hic1_undamaged, features = VariableFeatures(object = hic1_undamaged), npcs = 50)
# ElbowPlot(hic1_undamaged, ndims = 50)
hic1_undamaged <- FindNeighbors(hic1_undamaged, dims = 1:20)
hic1_undamaged <- FindClusters(hic1_undamaged, resolution = 0.6)
hic1_undamaged <- RunUMAP(hic1_undamaged, dims = 1:20)


################################# Saving processed Seurat object & aesthetics
saveRDS(hic1_undamaged, file = "seurat_obj.rds")
cols_clusters <- brewer.pal(12,"Paired")
cols_annot <- c("#FD8D3C", "#E31A1C", "#41B6C4", "#081D58", "#225EA8")


################################# Cluster annotations
new.cluster.ids <- c("Fibroblast", "VSMC", "VSMC", "VSMC", "PER-1", "PER-1", "Fibroblast", "PER-1",
                     "PER-1", "FAP", "PER-2", "Fibroblast")
names(new.cluster.ids) <- levels(hic1_undamaged)
hic1_undamaged <- RenameIdents(hic1_undamaged, new.cluster.ids)
hic1_undamaged@active.ident <- factor(hic1_undamaged@active.ident, levels = c("FAP", "Fibroblast", "VSMC", "PER-1", "PER-2"))
hic1_undamaged@meta.data$annotation <- Idents(hic1_undamaged)


################################# Visualizations
# UMAP
Idents(hic1_undamaged) <- hic1_undamaged@meta.data$annotation
hic1_undamaged@active.ident <- factor(hic1_undamaged@active.ident, levels = c("FAP", "Fibroblast", "VSMC", "PER-1", "PER-2"))
UMAP_annot <- DimPlot(hic1_undamaged, reduction = "umap", label = F, pt.size = 0.5, cols = cols_annot) + theme(text = element_text(face = "bold", size=16), axis.text=element_text(size=6), axis.title.x=element_blank(),
                                                                                                               axis.title.y=element_blank(), legend.position = "left") + NoAxes()
# Feature plot
FeaturePlot(hic1_undamaged, features = "X", pt.size = 0.1, reduction = "umap",
            min.cutoff = "0", max.cutoff = "q95", cols = c("lightgrey", "blue"), order = T) + theme(text = element_text(size=12), axis.text=element_text(size=12),
                                                                                                    axis.title.x=element_blank(),
                                                                                                    axis.title.y=element_blank(),
                                                                                                    plot.title=element_text(vjust=-1, face = "italic"),
                                                                                                    legend.position=c(0.05, 0.15),
                                                                                                    legend.direction = "horizontal",
                                                                                                    legend.key.size=unit(0.5, "cm"),
                                                                                                    legend.text=element_text(size = 8)) + NoAxes()
# Violin plot
VlnPlot(object = hic1_undamaged, features = "X", same.y.lims = TRUE, pt.size = 0.1, sort = F, cols = cols_annot) +
  theme(text = element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        plot.title=element_text(vjust=-1, face = "italic"))
        
# Heatmap
heatmap <- DoHeatmap(subset(hic1_undamaged, downsample = 100), features = "GENES", assay = "RNA", size = 1, group.colors = cols_annot, label = F)


################################# Peason's correlation analysis
cluster.averages <- AverageExpression(object = hic1_undamaged, return.seurat = TRUE, verbose = T, assays = "RNA")
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


################################# Differential gene expression
clusterMarkers <- FindMarkers(hic1_undamaged, ident.1 = c("X"), ident.2 = c("Y"), min.pct = 0.25, only.pos = T) # X vs. Y
allMarkers <- FindAllMarkers(hic1_undamaged, only.pos = TRUE, min.pct = 0.25) # One vs. all


################################# SCENIC - pre-processing
x = subset(hic1_undamaged)
save(x, file = "x.Robj")
dir.create("int") # Creating a folder to store intermediate objects: cellInfo and colVars

# Creating cellInfo, a dataframe annotating each cell by cluster annotation
cellInfo <- data.frame(celltype = as.vector(x@meta.data$celltype),
                       row.names = colnames(x@assays$RNA@data)) # colouring by what metadata
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Creating colVars, a dataframe assigning colors to each cluster 
cols <- brewer.pal(n = 12, name = "Paired")
colVars <- list(celltype=setNames(cols, unique(cellInfo$celltype)))
saveRDS(colVars, file="int/colVars.Rds")

# Creating an unnormalized, filtered expression matrix
x <- as.matrix(x@assays$RNA@counts)
exprMat <- x
dim(exprMat)

# Creating scenicOptions, storing information on cistarget databases and metadata
org= "mgi" # "hgnc" for human and "mgi" for mouse
dbDir= "/brcwork/rossi_lab/Henry/bioinf/cisTarget_databases" # set directory of the cisTarget database
scenicOptions <- initializeScenic(org= org, dbDir=dbDir, nCores=24)
scenicOptions@inputDatasetInfo$datasetTitle <- "title"
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Soft fitlering
# Counts per gene per cell
nCountsPerGene <- apply(exprMat, 1, sum) # calculating # counts per gene
minReads <- 3*.01*ncol(exprMat) # thresholding at 1% of cells expressing the gene (where 3 counts = cell expressed this gene)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)] # Keeping genes that are "expressed" in 1% of cells
length(genesLeft_minReads)

# Cells per gene per cell
nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0)) # calculating # cells expressing a gene
minSamples <- ncol(exprMat)*.01 # 1% of cells
genesLeft_minCells <- names(nCellsPerGene)[which(nCellsPerGene > minSamples)]
length(genesLeft_minCells)

# Genes kept after filtering
genesKept <- unique(genesLeft_minReads,genesLeft_minCells)

# Remove unconsidered genes in the database
motifRankings <- importRankings("/brcwork/rossi_lab/Henry/bioinf/cisTarget_databases/mm9-tss-centered-10kb-7species.mc9nr.feather") # loading database
genesInDatabase <- colnames(getRanking(motifRankings)) # calling genes in the database
genesKept <- genesKept[which(genesKept %in% genesInDatabase)] # find genes shared between the database and filtered data
length(genesKept)
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept")) # storing filtered genes in scenicOptions

# Performing filters
exprMat_filtered <- exprMat[genesKept, ]

# Computing Spearman correlation coefficient amongst genes kept
corrMat <- cor(t(exprMat_filtered), method="spearman")
allTFs <- getDbTfs(scenicOptions) # Grab all TFs in the database

# Since we are only interested in specifically correlations between TFs and downstream targets, we subset for TFs in the database
corrMat <- corrMat[which(rownames(corrMat) %in% allTFs),] # subsetting for TFs in the database
dim(corrMat)

# Saving correlation matrix to scenicOptions
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


################################# SCENIC - analysis
# Retrieving prepped data 
scenicOptions <- readRDS("int/scenicOptions.Rds")
load("x.Robj")
exprMat <- x@assays$RNA@counts
exprMat <- as.matrix(exprMat)
genesKept <- loadInt(scenicOptions, "genesKept")
exprMat_filtered <- exprMat[genesKept,]
exprMat_filtered <- log2(exprMat_filtered+1)

# Running SCENIC 
runGenie3(exprMat_filtered, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)


################################# SCENIC - post-analysis
# Heatmap of Binarized AUC
scenicOptions <- readRDS(file.choose())
cellInfo <- readRDS(file.choose())
minPerc <- .2 # filter out regulons that are not active (20% of cells per cluster) in any cell type

# Grabbing the binarized matrix
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")

# Grabbing cell identifications
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]

# Split into % per cell type
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

pheatmap::pheatmap(binaryActPerc_subset,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA, show_rownames = T, angle_col = 45)
