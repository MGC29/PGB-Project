#########################################################################################################
# Single cell RNA-seq analysis###########################################################################
#########################################################################################################

setwd("")
library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(readr)

FACS_files = list.files("TabulaMuris/FACS/", full.names = TRUE)
# Only considers cells from Brain tissue.
FACS_files <- grep("Brain" ,FACS_files, value = TRUE)
# Create a list of 3 elements, each a table of one of the cell types from both these tissues
raw.data.list <- list()
for (file in FACS_files){
  raw.data <- read.csv(file, row.names = 1)
  raw.data <- Matrix(as.matrix(raw.data), sparse = TRUE) # when no value -> "." instead of "0"
  raw.data.list <- append(raw.data.list, raw.data)
}

# Merge all the three tables using cbind function (bind by columns)
raw.data <- do.call(cbind, raw.data.list)

#PREPARE THE METADATA: Link the identifiers with the count matrix
meta.data <- read.csv("TabulaMuris/metadata_FACS.csv") 
head(meta.data)
# plate.barcode mouse.id  tissue subtissue FACS.selection mouse.sex
head(colnames(raw.data))
#cell.platebarcode.mouseId
#I need to repeat the metadata for each plate barcode times the number of cells from each plate
# Split the column names in the matrix by the "." divider and get the 2nd spot = Plate Barcode
plates <- str_split(colnames(raw.data),"[.]", simplify = TRUE)[,2]
# Put the row names as the names of the plates
# Create row names for each plate as many times as needed to match the matrix columns
rownames(meta.data) <- meta.data$plate.barcode
cell.meta.data <- meta.data[plates,]
rownames(cell.meta.data) <- colnames(raw.data)

#Spiking genes > exeperimentally determine the normalization factor of my samples
# 92 genes from bacteria, we put known RNA concentration of this genes, sequence it
# and see if the resulting number of reads is equal. If one cell shows more of this reads,
# it has to be normalized. (If one cell x2 reads > divide by 2)

#Remove the ERCC genes
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
#Proportion of the reads going to the ERCC elements compared with the total values
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,]

######################################################################################################
# Create a Seurat object #############################################################################
######################################################################################################

tiss <- CreateSeuratObject(counts = raw.data) #23341 features across 10561 samples within 1 assay 
dim(tiss)
colnames(tiss)
head(tiss@meta.data) #nCount_RNA nFeature_RNA
tiss <- AddMetaData(object = tiss, cell.meta.data) # nCount_RNA nFeature_RNA + cols of metadata
tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
# Last cmd is equal to tiss$my_new_column <- percent.ercc

# Number of genes per cell: median around 2500
# However, there are a few cells containing 0 genes
VlnPlot(tiss, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# We subset our data based on the number of genes and counts
# Only keep cells containing >500 genes and > 50000 counts
tiss <- subset(tiss, subset = nFeature_RNA >500 & nCount_RNA > 50000)
VlnPlot(tiss, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Normalize the counts by 1M as an scaling factor
tiss <- NormalizeData(object = tiss, scale.factor = 1e6) #tiss@assays$RNA@layers$data
tiss <- ScaleData(object = tiss) #tiss@assays$RNA@layers$scale.data

###########################################################################################
# We will calculate PCA but only on the top variable genes: ###############################
###########################################################################################
tiss <- FindVariableFeatures(object = tiss)
VariableFeatures(tiss)
top10 <- head(VariableFeatures(tiss), 10)
plot1 <- VariableFeaturePlot(tiss)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Choose a number of PCs based on the elbow plot.
tiss <- RunPCA(tiss, features = VariableFeatures(object = tiss))
ElbowPlot(tiss, ndims = 50)

#Find clusters
tiss <- FindNeighbors(tiss, dims = 1:43)
tiss <- FindClusters(tiss, resolution = 0.5)
# Number of communities: 1
table(tiss$seurat_clusters)
table(tiss$RNA_snn_res.0.5)

#Run UMAP to visualize them:
tiss <- RunUMAP(tiss, dims = 1:43)
DimPlot(tiss, reduction = "umap", group.by = 'tissue')
DimPlot(tiss, reduction = "umap", group.by = 'subtissue')
DimPlot(tiss, reduction = "umap", label = TRUE)

#Find marker genes
# Look at the level of expression as well as the number of cells expressing a gene
tiss.markers <- FindAllMarkers(tiss, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Rows = Genes Cols = p_val  avg_log2FC  pct.1   pct.2   p_val_adj
FeaturePlot(tiss, features = "Neurod2")
RidgePlot(tiss, features = "Neurod2", group.by = "seurat_clusters")
VlnPlot(tiss, features = "Neurod2", group.by = "seurat_clusters")
#Enrichr > we copy the more expressed genes for a cluster (marker genes) and get a type of cell for which those genes are markers

#Order by fold change and plot their expression using the function FeaturePlot
#Find out real annotation
anno <- read_csv("TabulaMuris/annotations_FACS.csv")
tiss@meta.data$cell <- rownames(tiss@meta.data)
meta2 <- tiss@meta.data %>% left_join(anno[,c(1,3)], by='cell')
tiss <- AddMetaData(object = tiss, meta2$cell_ontology_class, col.name = "cell_ontology_class")
tiss@meta.data$cell_ontology_class[is.na(tiss@meta.data$cell_ontology_class)] <- "unknown"
tiss$cell_ontology_class <-as.factor(tiss$cell_ontology_class)

neuron.markers <- FindMarkers(tiss, ident.1 = "neuron", min.pct = 0.25) #, test.use = "") default is wilcoxon

DimPlot(tiss, reduction = "umap", group.by = 'cell_ontology_class')

FeaturePlot(tiss, features = "Neurod2")
RidgePlot(tiss, features = "Neurod2", group.by = "cell_ontology_class")
VlnPlot(tiss, features = "Neurod2", group.by = "cell_ontology_class")
