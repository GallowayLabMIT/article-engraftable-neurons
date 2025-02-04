###############################################################################
# Before starting, place a txt file in the same folder this .R file is in called
# 'datadir.txt' which points to the folder which contains the aligned/processed
# files (GSE: )
# e.g. In 'C:\Users\username\Documents\GitHub\article-engraftable-neurons\r'
# Place 'datadir.txt' containing:
# 'C:\Users\username\Documents\article-engraftable-neurons_data'
###############################################################################


# Load packages and everything else ---------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(modules)
  library(sctransform)
  library(ggplot2)
})

# Load tools from modules
tools <- modules::use("src/modules/gene_analysis_tools.R")
dir_tools <- modules::use("src/modules/dir_tools.R")
original_wd <- dir_tools$set_datadir() # Saves original wd and goes to
                                       # specified path in datadir.txt

###############################################################################


#' Loads filtered_feature_bc_matrix from 10X -> CellRanger as Seurat obj
#' @param data_10X_wd A string representing the directory where the 10X data is.
#' @param dirname A string representing the directory where the
#'                outs/filtered_feature_bc_matrix folder is
#' @param condition A string representing the condition.
#' @param dpi A numeric representing dpi cells come from
#' @return A Seurat object with the percent mitochondrial reads
load_10X_data <- function(data_10X_wd, dirname, condition, dpi) {
  
  print(paste("Loading" , dirname, "..."))
  
  # Load count matrix file
  data_dir <- paste(data_10X_wd, dirname, sep="/")
  dataCounts <- Read10X(data.dir=data_dir)
  
  # Create Seurat Object
  countObj = CreateSeuratObject(counts=dataCounts, project=dirname)
  
  # Add metadata
  countObj <- AddMetaData(countObj, condition, "condition")
  countObj <- AddMetaData(countObj, dpi, "dpi")
  countObj <- AddMetaData(countObj, paste0(condition,'-',dpi), "condDpi")
  countObj <- AddMetaData(countObj, paste0(condition,'-',dpi,'-','10X'), "condDpiData")
  countObj <- AddMetaData(countObj, '10X', "datatype") 
  countObj[["percent.mt"]] <- PercentageFeatureSet(countObj, pattern = "^mt-")
  
  print("done!")
  return(countObj)
}


#-------------Load processed data ----------------

# Save current directory and get 10X and data
data_10X_wd <- getwd()

# Directory to store outputs
output_dir <- '2023.10.16_embMN-filter' 
all_output_dir <-  paste0(original_wd, '/../r-outputs/', output_dir, "/")

# Directory for image outputs
img_dir <- paste0(all_output_dir, 'img_output/')
# Make it if it doesn't exist
if (!file.exists(img_dir)) {
  if (!file.exists(all_output_dir)) {
    dir.create(all_output_dir)
  }
  dir.create(img_dir)
}


#-------------Get data ----------------


# Saves Seurat object in compressed Rds for easier loading in future
seuratObj_filename <- paste0(all_output_dir, output_dir, '.Rds')

# If already saved as Rds, load
if (file.exists(seuratObj_filename)) {
  seuratObj <- readRDS(seuratObj_filename)
# Else, load data and convert to Seurat objects
} else {
  
  # Load count matrices
  all_obj <- c(
    load_10X_data(data_10X_wd, 'EmbMN-1', 'EmbMN', 14), # (data_10X_wd, dirname, condition, dpi)
    load_10X_data(data_10X_wd, 'EmbMN-2', 'EmbMN', 14),
    load_10X_data(data_10X_wd, 'MEFs', 'MEF', 0)
  )
  
  # Get sample names for cell.ids
  sample_ids <- character() # Initialize char vector
  for (i in 1:length(all_obj)) {
    # eg. "EmbMN-1" (orig.ident)
    sample_id <- paste0(unique(all_obj[[i]]@meta.data$orig.ident), "-", unique(all_obj[[i]]@meta.data$datatype))
    sample_ids <- append(sample_ids, sample_id)
  }
  
  # Merge into single Seurat object with a project name
  seuratObj <- merge(all_obj[[1]], y=all_obj[-1], add.cell.ids=sample_ids)
  
  # Clear unused vars
  rm(sample_ids, sample_id, all_obj)
  gc()
  saveRDS(seuratObj, seuratObj_filename)
  
}

#------------- QC steps and filtering ----------------

# plot the QC metrics in a violin plot
VlnPlot(seuratObj, group.by = 'condDpiData',
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, ) # &
  # theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

ggsave( paste0(img_dir, 'QC1_prefilt_1.png'), width=12, height=5, dpi=300)

# plot scatter QC 
plot1 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='condDpiData', jitter=TRUE, pt.size=0.8)
plot2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='condDpiData', jitter=TRUE, pt.size=0.8)
plot1 + plot2
ggsave(paste0(img_dir, 'QC1_prefilt_2.png'), width=13, height=5, dpi=300)

# How many cells are recovered before filtering
sprintf("%i cells before filtering", ncol(seuratObj))
# Filter cells with <8% mitochondria reads, captured genes between [200, 10,000]
seuratObj_filt <- subset(seuratObj, subset=(nFeature_RNA>200 & nFeature_RNA<10e3 & percent.mt<8))
# How many cells are recovered after filtering
sprintf("%i cells after filtering", ncol(seuratObj_filt))
sprintf("%i cells (%.0f%% of all cells) were filtered out", ncol(seuratObj)-ncol(seuratObj_filt), (ncol(seuratObj)-ncol(seuratObj_filt))/ncol(seuratObj)*100)

# How many cells per condition
table(seuratObj@meta.data$condDpiData, seuratObj@meta.data$orig.ident)
table(seuratObj_filt@meta.data$condDpiData, seuratObj_filt@meta.data$orig.ident)

# Replot the QC metrics after filtering
VlnPlot(seuratObj_filt, group.by = 'condDpiData',
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
ggsave(paste0(img_dir, 'QC1_postfilt_1.png'), width=13, height=5, dpi=300)
plot1 <- FeatureScatter(seuratObj_filt, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='condDpiData', jitter=TRUE, pt.size=0.8)
plot2 <- FeatureScatter(seuratObj_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='condDpiData', jitter=TRUE, pt.size=0.8)
plot1 + plot2
ggsave(paste0(img_dir, 'QC1_postfilt_2.png'), width=13, height=5, dpi=300)


# ------------- Normalize data -----------------

# Saves normalized Seurat object in compressed Rds for easier loading in future
seuratObjNorm_filename <- paste0(all_output_dir, output_dir, '_NORM.Rds')
rm(seuratObj)
gc()

# If already saved as Rds, load
if (file.exists(seuratObjNorm_filename)) {
  obj <- readRDS(seuratObjNorm_filename)
  # Else, load data and convert to Seurat objects
} else {

  # Normalize using SCTransform
  options(future.globals.maxSize = 8000 * 1024^2) # may need this
  obj  <- SCTransform(seuratObj_filt, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

  # Save
  saveRDS(obj, seuratObjNorm_filename)

}

rm(seuratObj_filt)
gc()

#------------- Linear dimensional reduction! ----------------

obj <- RunPCA(obj, verbose=FALSE)

# Top 5 PCA dim and features
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA dimension loadings
VizDimLoadings(obj, dims = 1:2, reduction = "pca")
ggsave(paste0(img_dir, 'QC2_PC1-2feature.png'), width=13, height=5, dpi=300)

# Visualize PCA on graph by orig ident
DimPlot(obj, reduction="pca") + xlim(-125, 125) + ylim(-75, 100)
ggsave(paste0(img_dir, 'QC2_PC1v2_ident.png'), width=13, height=5, dpi=300)

# Do by condition and dpi
DimPlot(obj, reduction="pca", group.by='condDpi') + xlim(-125, 125) + ylim(-75, 100)
ggsave(paste0(img_dir, 'QC2_PC1v2_condDpi.png'), width=13, height=5, dpi=300)

# Visualize PCA via heatmap
DimHeatmap(obj, dims=1, cells = 500, balanced = TRUE)
DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE) # Save as 1000x2000

#-------------  Determine dimensionality of data set ----------------

# Elbow plot
plot1 <- ElbowPlot(obj, ndims=50)
plot1
# Looks like 40 PCs is good?
# Note: Explored only 40 PCs
numPC <- 40
plot2 <- geom_vline(xintercept = numPC, linetype = "dashed", color = 'red', linewidth=1)
plot1 + plot2
ggsave(paste0(img_dir, 'QC2_PCA_elbow.png'), width=6, height=5, dpi=300)

#------------- Cluster cells and run non-lin dedim ----------------

# Use same PCs as above
obj <- FindNeighbors(obj, dims = 1:numPC)

# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells.
# Optimal resolution often increases for larger datasets.
# Use resolution above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
obj <- FindClusters(obj, resolution=0.2)

# Perform UMAP
obj <- RunUMAP(obj, dims = 1:numPC)

# Plot by cluster
DimPlot(obj, reduction = "umap", label=TRUE) + NoLegend()
ggsave(paste0(img_dir, '1-UMAP_byCluster.png'), width=4, height=4, dpi=300)

# Plot by condition
DimPlot(obj, reduction = "umap", group.by='condDpi', label=FALSE)
ggsave(paste0(img_dir, '1-UMAP_byCondDpi.png'), width=6, height=4, dpi=300)
DimPlot(obj, reduction = "umap", group.by='condDpi', label=FALSE) + NoLegend()
ggsave(paste0(img_dir, '1-UMAP_byCondDpi-nolabel.png'), width=4, height=4, dpi=300)

DimPlot(obj, reduction = "umap", group.by='orig.ident', label=FALSE)
ggsave(paste0(img_dir, '1-UMAP_byOrigIdent.png'), width=6, height=5, dpi=300)

# How many cells per condition
table(obj@active.ident)
table(obj@active.ident, obj@meta.data$orig.ident)

# ----- Look at general gene expression across clusters ---------------------

# Directory to store feature plots
feature_dir <- paste0(img_dir, 'features/')
# Make it if it doesn't exist
if (!file.exists(feature_dir)) {
  dir.create(feature_dir)
}

# Plot all genes of interest
genes.of.interest <- c(
  'Isl1', 'Lhx3', 'Mnx1', 'v.EGFP'
)
# Feature plot
for (gene in genes.of.interest) {
  FeaturePlot(obj, features=gene) & 
    scale_color_viridis_c()
  ggsave(paste0(feature_dir, gene, '-feature.png'), width=3, height=2, dpi=300)
}
# Violin plot
for (gene in genes.of.interest) {
  VlnPlot(obj, features=gene)
  ggsave(paste0(feature_dir, gene, '-vln.png'), width=10, height=5, dpi=300)
}

DotPlot(obj, features=genes.of.interest) & 
  scale_color_viridis_c()
ggsave(paste0(feature_dir, 'dotplot-features.png'), width=20, height=6, dpi=300)


# ------------------ Extract MMCs and save ----------------------

# Saves MMCs in compressed Rds for easier loading in future
objMMC_filename <- paste0(all_output_dir, output_dir, '_NORM_embMN-MMC.Rds')

# If not saved as Rds, save
if (!file.exists(objMMC_filename)) {
  # MMC is cluster 10 (Isl1, Lhx3, Mnx1+)
  MMCobj <- subset(x=obj, ident=10)
  table(MMCobj@active.ident, MMCobj@meta.data$orig.ident)
  saveRDS(MMCobj, objMMC_filename)
}









































