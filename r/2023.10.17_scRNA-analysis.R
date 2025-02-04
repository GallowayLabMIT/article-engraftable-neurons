###############################################################################
# Before starting, place a txt file in the same folder this .R file is in called
# 'datadir.txt' which points to the folder which contains the aligned/processed
# files (GSE: )
# e.g. In 'C:\Users\username\Documents\GitHub\article-engraftable-neurons\r'
# Place 'datadir.txt' containing:
# 'C:\Users\username\Documents\article-engraftable-neurons_data'
#
#
#
# ALso run 2023.10.16_embMN-filter.R before this to extract MMC motor neurons
#
###############################################################################


# Load packages and everything else ---------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(speckle)
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
output_dir <- '2023.10.17_scRNA-analysis' 
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
    load_10X_data(data_10X_wd, '6FDDRR-14dpi-1-iMN', '6FDDRR', 14),
    load_10X_data(data_10X_wd, '6FDDRR-14dpi-2-iMN', '6FDDRR', 14),
    load_10X_data(data_10X_wd, 'LNI-DDRR-iMN', 'LNI-DDRR', 14),
    load_10X_data(data_10X_wd, 'LNI-DDRR-nosupplements', 'LNI-DDRR-noNT', 14),
    load_10X_data(data_10X_wd, 'LNI-DDRR-4dpi', 'LNI-DDRR', 4),
    load_10X_data(data_10X_wd, 'MEFs', 'MEF', 0)
  )
  
  # Get sample names for cell.ids
  sample_ids <- character() # Initialize char vector
  for (i in 1:length(all_obj)) {
    # eg. "EmbMN-1" (orig.ident) + "10X" (datatype)
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


# ------------ Load in filtered embMNs (MMCs only) ---------------

# Load in MMC file
MMC_dir <- '2023.10.16_embMN-filter'
objMMC_filename <- paste0(original_wd, '/../r-outputs/', MMC_dir, "/", MMC_dir, "_NORM_embMN-MMC.Rds")

# If saved as Rds, read
if (file.exists(objMMC_filename)) {
  MMCobj <- readRDS(objMMC_filename)
  # Rename metadata to make easier to read
  MMCobj$condition <- 'EmbMN'
  MMCobj$dpi <- 14
  MMCobj$condDpi <- 'EmbMN-14'
  MMCobj$condDpiData <- 'EmbMN-14-10X'
  Idents(object=MMCobj) <- 'EmbMN-MMC'
}
table(MMCobj@active.ident, MMCobj@meta.data$orig.ident)

# Add MMCs to seuratObj
seuratObj_filt <- merge(MMCobj, seuratObj_filt)
table(seuratObj_filt@active.ident, seuratObj_filt@meta.data$orig.ident)

rm(MMCobj)
gc()

# Replot the QC metrics after merging
VlnPlot(seuratObj_filt, group.by = 'condDpiData',
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
ggsave(paste0(img_dir, 'QC2_postfilt_1.png'), width=13, height=5, dpi=300)
plot1 <- FeatureScatter(seuratObj_filt, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by='condDpiData', jitter=TRUE, pt.size=0.8)
plot2 <- FeatureScatter(seuratObj_filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by='condDpiData', jitter=TRUE, pt.size=0.8)
plot1 + plot2
ggsave(paste0(img_dir, 'QC2_postfilt_2.png'), width=13, height=5, dpi=300)


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
# seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj), npcs=200)

# Top 5 PCA dim and features
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA dimension loadings
VizDimLoadings(obj, dims = 1:2, reduction = "pca")
ggsave(paste0(img_dir, 'QC3_PC1-2feature.png'), width=13, height=5, dpi=300)

# Visualize PCA on graph by orig ident
DimPlot(obj, reduction="pca") + xlim(-125, 125) + ylim(-75, 100)
ggsave(paste0(img_dir, 'QC3_PC1v2_ident.png'), width=13, height=5, dpi=300)
# Do by condition and dpi
DimPlot(obj, reduction="pca", group.by='condDpi') + xlim(-125, 125) + ylim(-75, 100)
ggsave(paste0(img_dir, 'QC3_PC1v2_condDpi.png'), width=13, height=5, dpi=300)

# Separate some out so easier to see
smallObj <- subset(x=obj, condition==c('MEF', 'EmbMN'))
DimPlot(smallObj, reduction="pca") + xlim(-125, 125) + ylim(-75, 100)
ggsave(paste0(img_dir, 'QC3_PC1v2_ident_1.png'), width=13, height=5, dpi=300)

smallObj <- subset(x=obj, condition==c('6FDDRR', 'LNIx3HARIDD', 'LNIx3HARIDD-noNT'))
DimPlot(smallObj, reduction="pca") + xlim(-125, 125) + ylim(-75, 100)
ggsave(paste0(img_dir, 'QC3_PC1v2_ident_2.png'), width=13, height=5, dpi=300)
DimPlot(smallObj, reduction="pca", group.by='condDpi') + xlim(-125, 125) + ylim(-75, 100)
ggsave(paste0(img_dir, 'QC3_PC1v2_condDpi_2.png'), width=13, height=5, dpi=300)

rm(smallObj)
gc()

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
ggsave(paste0(img_dir, 'QC3_PCA_elbow.png'), width=6, height=5, dpi=300)

#------------- Cluster cells ----------------

# Use same PCs as above
obj <- FindNeighbors(obj, dims = 1:numPC)

# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells.
# Optimal resolution often increases for larger datasets.
# Use resolution above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
# NOTE: Found 0.4 is okay, could probably go lower
obj <- FindClusters(obj, resolution=0.2)



#------------- Run non-linear dimensional reduction ----------------

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

#------------- Investigate satellite MEF cluster ----------------
# A small satellite cluster from the MEF sample appears in cluster 8

smallObj_MEFs <- subset(x=obj, idents=c(4, 8))
# plot the QC metrics in a violin plot

VlnPlot(smallObj_MEFs,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, )
ggsave(paste0(img_dir, '2-MEF-clusters_QC.png'), width=10, height=8, dpi=300)

# What separates out the MEF clusters?
cluster4v8.markers <- FindMarkers(smallObj_MEFs, ident.1=4, ident.2=8)
head(cluster4v8.markers, n=20)

# Looks like some type of immune cell rather than MEF
VlnPlot(smallObj_MEFs,
            features = c("Col1a1", "Spi1", "Laptm5", "Csf1r", 'Lcp1'))
ggsave(paste0(img_dir, '2-MEF-clusters_vlnplt-1.png'), width=10, height=8, dpi=300)
rm(smallObj_MEFs)

# Look at just MEFs
obj.MEFs <- subset(x=obj, condDpi=='MEF-0')
obj.MEFs <- FindNeighbors(obj.MEFs, dims = 1:numPC)

# Adj resolution to be low to look at pretty separate things
obj.MEFs <- FindClusters(obj.MEFs, resolution=0.15)

# Perform UMAP
obj.MEFs <- RunUMAP(obj.MEFs, dims = 1:numPC)

# Plot by cluster
DimPlot(obj.MEFs, reduction = "umap", label=TRUE) + NoLegend()
ggsave(paste0(img_dir, '2-MEF-clusters_UMAP.png'), width=4, height=4, dpi=300)

# What is small cluster 2?
# Again see immune and hematopoietic markers
cluster2.markers <- FindMarkers(obj.MEFs, ident.1=2)
head(cluster2.markers, n=20)
VlnPlot(obj.MEFs,
        features = c("Col1a1", "Rac2", "C5ar1", "Spi1", "Ncf1", "Laptm5", "Csf1r"))
ggsave(paste0(img_dir, '2-MEF-clusters_vlnplt_2.png'), width=10, height=8, dpi=300)

# Eliminate cluster 2
# assume they came in during MEF isolation/dissection
DimPlot(obj.MEFs, reduction = "umap", label=TRUE) + NoLegend()
ggsave(paste0(img_dir, '2-MEF-clusters_vlnplt_beforeElim.png'), width=4, height=4, dpi=300)

obj.MEFs <- subset(x=obj.MEFs, idents=2, invert=TRUE)
DimPlot(obj.MEFs, reduction = "umap", label=TRUE) + NoLegend()
ggsave(paste0(img_dir, '2-MEF-clusters_vlnplt_afterElim.png'), width=4, height=4, dpi=300)


#----- Re-run without satellite MEF cluster ------
otherData <- subset(x = obj,  subset = condition != 'MEF')
obj <- merge(obj.MEFs, otherData)
rm(obj.MEFs, otherData)
gc()

# Renormalize
seuratObjNorm_filename_MEFadj <- paste0(all_output_dir, output_dir, '_NORM_MEFadj.Rds')

# If already saved as Rds, load
if (file.exists(seuratObjNorm_filename_MEFadj)) {
  obj <- readRDS(seuratObjNorm_filename_MEFadj)
  # Else, load data and convert to Seurat objects
} else {
  
  # Normalize using SCTransform
  obj  <- SCTransform(obj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
  
  # Cluster
  obj <- RunPCA(obj, verbose=FALSE)
  obj <- FindNeighbors(obj, dims = 1:numPC)
  obj <- FindClusters(obj, resolution=0.2)
  obj <- RunUMAP(obj, dims = 1:numPC)
  
  # Reorder conditions
  obj$condDpi <- factor(obj$condDpi,
                        levels = c("MEF-0",
                                   "LNI-DDRR-4",
                                   "LNI-DDRR-noNT-14",
                                   "LNI-DDRR-14",
                                   "6FDDRR-14",
                                   "EmbMN-14"))
  
  # Save
  saveRDS(obj, seuratObjNorm_filename_MEFadj)

}

# Plot by cluster
DimPlot(obj, reduction = "umap", label=TRUE)
ggsave(paste0(img_dir, '3-UMAP_byCluster.png'), width=6, height=4, dpi=300)
ggsave(paste0(img_dir, '3-UMAP_byCluster-nolabel.png'), width=4, height=4, dpi=300)

# Plot by condition
DimPlot(obj, reduction = "umap", group.by='condDpi', label=FALSE)
ggsave(paste0(img_dir, '3-UMAP_byCondDpi.png'), width=6, height=4, dpi=300)
DimPlot(obj, reduction = "umap", group.by='condDpi', label=FALSE) + NoLegend()
ggsave(paste0(img_dir, '3-UMAP_byCondDpi-nolabel.png'), width=4, height=4, dpi=300)

# Plot by original identity
DimPlot(obj, reduction = "umap", group.by='orig.ident', label=FALSE)
ggsave(paste0(img_dir, '3-UMAP_byOrigIdent.png'), width=8, height=4, dpi=300)


#------------- Finding cluster biomarkers ----------------

# Directory to store biomarker plots
biomarker_dir <- paste0(img_dir, 'biomarkers/')
# Make it if it doesn't exist
if (!file.exists(biomarker_dir)) {
  dir.create(biomarker_dir)
}

# How many cells per condition
table(obj@active.ident, obj@meta.data$condition)


# -------------- What are the 2 clusters for LNI-DDRR-4dpi ---------------------

# What markers separate cluster 5 from 0 in LNI-DDRR-4
cluster5.markers <- FindMarkers(obj, ident.1=5, ident.2=0)
# Plot top 20 genes from cluster 5
# Just grab clusters 5, 0 and MEFs (4)/EmbMN (7) for ref
smallobj <- subset(x=obj, idents=c(0, 4, 5, 7))
smallobj@active.ident <- factor(smallobj@active.ident,
                                levels = c(4, 0, 5, 7))

# Sort by FC value
custom_clus5v0genes <- head(cluster5.markers, n=20)
custom_clus5v0genes <- rbind(custom_clus5v0genes, cluster5.markers['Isl1', ])
custom_clus5v0genes <- custom_clus5v0genes[order(custom_clus5v0genes$avg_log2FC, decreasing=FALSE), ]
DotPlot(smallobj,  features=rownames(custom_clus5v0genes)) & 
  scale_color_viridis_c() &
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(biomarker_dir, '0v5-dotplot-features_rev.png'), width=10, height=4, dpi=300)

# -------------- What is central cluster 6? ---------------------

# What markers make up cluster 6, the proliferative one
cluster6.markers <- FindMarkers(obj, ident.1=6)
head(cluster6.markers, n=20)
# Cluster 1: 6F + DDRR
cluster6v1.markers <- FindMarkers(obj, ident.1=6, ident.2=1)
head(cluster6v1.markers, n=20)
# Cluster 2: LNI + RIDD + no NT
cluster6v2.markers <- FindMarkers(obj, ident.1=6, ident.2=2)
head(cluster6v2.markers, n=20)
# Cluster 3: LNI + RIDD
cluster6v3.markers <- FindMarkers(obj, ident.1=6, ident.2=3)
head(cluster6v3.markers, n=20)

for (gene in rownames(head(cluster6.markers, n=20))) {
  VlnPlot(obj, features=gene, idents=c(1,2,3, 6) # , cols=c('#E4E1E3', '#F6222E', '#FE00FA', '#1CFFCE')
                   )
  ggsave(paste0(biomarker_dir, 'cluster6-', gene, '-vln.png'), width=4, height=3, dpi=300)
}

# Plot how things of interest change between clusters 6 and other MNs
genes.of.interest <- c(
  'Mki67', 'Top2a', 'Top1','Cdkn1a', 'Trp53', # Cell cycle
  'Col1a1', # Fib genes 
  'Scg2', 'Chgb', 'Snca', 'Sncg', # Neurosignaling
  'Isl2', 'Mnx1', # MMC genes (Mnx1, Isl1/2, Lhx3)
  'Foxp1', 'Lhx1', # LMC genes (Foxp1, Isl1/2, Lhx1, Mnx1)
  'Etv1', # HMC (Mnx1, Isl1, Etv1)
  'Nos1', 'Zeb2', # PGC (Isl1, nNos, pSmad, Zeb2)
  'Chat', 'Mapt', 'Rbfox3', # Neuronal markers
  'Tubb3', 'Elavl3', # Delile paper - neuron
  'Pou4f1', # Delile paper - dl3 vs. MN
  'Neurog2', 'Isl1', 'Lhx3', 'Ascl1', 'Pou3f2', 'Myt1l', 'v.EGFP' # 6F
)
# Feature plot
for (gene in genes.of.interest) {
  VlnPlot(obj, features=gene, idents=c(1,2,3, 6) #, 1:6F+DDRR, 2:LNI-DDRR-nosupplements, 3:LNI-DDRR
          )
  ggsave(paste0(biomarker_dir,'cluster6-GOI-', gene, '-feature.png'), width=3, height=2, dpi=300)
}

#------------- Investigate cluster 6 (incomplete iMN) ----------------

# Look at rough overall breakdown
cluster.prop <- prop.table(table(Idents(obj), obj$condDpi), margin = 2)
print(cluster.prop)

# Plot cell type proportions for just 14 dpi conditions
dpi14_cond <- c('LNI-DDRR-noNT-14', 'LNI-DDRR-14', '6FDDRR-14', 'EmbMN-14')
smallObj <- subset(x=obj,
                   condDpi==dpi14_cond)
# Reorder conditions
smallObj$condition <- factor(smallObj$condition,
                      levels = c('LNI-DDRR-noNT',
                                 'LNI-DDRR',
                                 '6FDDRR',
                                 'EmbMN'))

# Plot cluster percents for each condition

plotCellTypeProps(clusters=smallObj$seurat_clusters, sample=smallObj$condition) +
  geom_bar(stat="identity", color="black") +
  theme_classic() +
  scale_x_discrete(labels=c('6FDDRR'='6F\n+ DDRR',
                            'EmbMN'="EmbMN",
                            'LNI-DDRR'="LNI\n+ DDRR",
                            'LNI-DDRR-noNT'="LNI\n+ DDRR\n(no NTR)"))
ggsave(paste0(img_dir, '4-clusters_14dpi.png'), width=4.5, height=3, dpi=300)
cluster.prop

# Plot % cluster 6 (incomplete iMN) for each condition
# Convert into df, then barplot
df <- data.frame(
  name=dpi14_cond,  
  cluster0=cluster.prop["0",dpi14_cond],
  cluster1=cluster.prop["1",dpi14_cond],
  cluster2=cluster.prop["2",dpi14_cond],
  cluster3=cluster.prop["3",dpi14_cond],
  cluster4=cluster.prop["4",dpi14_cond],
  cluster5=cluster.prop["5",dpi14_cond],
  cluster6=cluster.prop["6",dpi14_cond],
  cluster7=cluster.prop["7",dpi14_cond]
)
ggplot(df, aes(x=name, y=cluster6*100)) + 
  geom_bar(stat = "identity", color='black', fill='grey') +
  xlab("") + ylab("(%) Incomplete iMN\n(cluster 6)") +
  theme_classic() +
  scale_x_discrete(labels=c('6FDDRR-14'='6F\n+ DDRR',
                            'EmbMN-14'="EmbMN",
                            'LNIx3HARIDD-14'="LNI\n+ DDRR",
                            'LNIx3HARIDD-noNT-14'="LNI\n+ DDRR\n(no NTR)")
  )
ggsave(paste0(img_dir, '4-cluster7-bar_14dpi.png'), width=3, height=2, dpi=300)

# ----- Look at gene expression without cluster 6 (incomplete iMN) -------------

# Eliminate cluster 6
DimPlot(obj, reduction = "umap", label=TRUE)
ggsave(paste0(img_dir, '5-UMAP-beforeElimCluster6.png'), width=4, height=4, dpi=300)

obj <- subset(x=obj, idents=6, invert=TRUE)
DimPlot(obj, reduction = "umap", label=TRUE)
ggsave(paste0(img_dir, '5-UMAP-afterElimCluster6.png'), width=4, height=4, dpi=300)

# Directory to store feature plots
final_feature_dir <- paste0(img_dir, 'featuresAfterElim6/')
# Make it if it doesn't exist
if (!file.exists(final_feature_dir)) {
  dir.create(final_feature_dir)
}

# Look at some general features (note: Pou3f2 = Brn2)
FeaturePlot(obj,
            features = c("Ascl1", "Pou3f2", "Myt1l", "Neurog2", "Isl1", "Lhx3", "v.EGFP")) & 
  scale_color_viridis_c()
ggsave(paste0(img_dir, '6-UMAP_6F+eGFP_norm.png'), width=10, height=8, dpi=300)

# Plot UMIs
FeaturePlot(obj, features='nCount_RNA') & 
  scale_color_viridis_c()
ggsave(paste0(final_feature_dir, '6-UMAP_UMI.png'), width=3, height=2, dpi=300)


# Plot all genes of interest
genes.of.interest <- c(
  'Mki67', 'Top2a', 'Top1','Cdkn1a', 'Trp53', # Cell cycle
  'Col1a1', # Fib genes 
  'Scg2', 'Chgb', 'Snca', 'Sncg', # Neurosignaling
  'Isl1', 'Mnx1', # MMC genes
  'Chat', 'Mapt','Tubb3', 'Elavl3','Rbfox3', # Neuronal markers
  'Neurog2', 'Isl1', 'Lhx3', 'Ascl1', 'Pou3f2', 'Myt1l', 'v.EGFP' # 6F
)
# Feature plot
for (gene in genes.of.interest) {
  FeaturePlot(obj, features=gene) & 
    scale_color_viridis_c()
  ggsave(paste0(final_feature_dir, gene, '-feature.png'), width=3, height=2, dpi=300)
}
# Violin plot
for (gene in genes.of.interest) {
  VlnPlot(obj, features=gene, group.by='condDpi')
  ggsave(paste0(final_feature_dir, gene, '-vln.png'), width=7, height=5, dpi=300)
}