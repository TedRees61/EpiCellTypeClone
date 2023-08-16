#
#   Code broadly based on https://satijalab.org/seurat/articles/atacseq_integration_vignette.html
#
#   Section 1 - reference data used for annotation transfer
#
#   Section 2 - Benchmark RNA part and labelling for Ground Truth
#
#   Section 3 - Benchmark ATAC part and Labelling for test results
#

# Section 1 Code to read in Reference object and process
# Installs when required
#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
#remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
#remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
#remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
#4remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
#library(Rsamtools) not sure this is needed
library(Seurat)
library(SeuratObject)

library(tidyverse)
#
# read in Reference data file
#
Cereb <- readRDS("/Users/tedrees/ProjectTest/Linardsson/CB Cerebellar Vermis.rds")
head(Cereb)
library(ggplot2)
library(cowplot)
library(patchwork)
#

#
# View QC Metrics
#
range(Cereb$total_genes)
range(Cereb$total_UMIs)
range(Cereb$fraction_mitochondrial)

#
# Populate Counts - since this is not in the original file in the right place..:-(
#
Cereb@assays$RNA@counts <- Cereb@assays$RNA@data

#
# Run normalisation & Find Variable Features Step
#

Cereb <- NormalizeData(Cereb, verbose = TRUE)

# Run the standard workflow for visualization and clustering
Cereb <- FindVariableFeatures(Cereb, selection.method = "vst", nfeatures = 2000)
Cereb <- ScaleData(Cereb, verbose = TRUE)

Cereb <- RunPCA(Cereb, npcs = 30, verbose = TRUE)
#
# Warning! next step took ~10 hours
#
Cereb <- RunUMAP(Cereb, reduction = "pca", dims = 1:30, verbose = TRUE)
p1 <- DimPlot(Cereb, reduction = "umap", group.by = "supercluster_term")
p1
p2 <- DimPlot(Cereb, reduction = "umap", group.by = "cell_type") 
p2 
# Checkpoint Save & Read for restart
#
saveRDS(Cereb,file = "~/ProjectTest/10X Filtered PBMC/CerebPostUMAP.rds")
#
# Restart Step if needed
#
Cereb <- readRDS("/Users/tedrees/ProjectTest/10X Filtered PBMC/CerebPostUMAP.rds")
Cereb = UpdateSeuratObject(object = Cereb)
#
# Section 2 - Do the same for query file (multiome) This time read in from cellranger output
#
# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("~/ProjectTest/10X Filtered PBMC/human_brain_3k_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

#
# Convert Rownames to ensemble format
#
GeneLookup<-read.csv("/Users/tedrees/ProjectTest/GeneLookup.csv")
Cereb10XGeneNames <- rownames(rna_counts)
#
# Convert Gene Names to Ensembl format
# Lookup file created from source with ~55k genes
#
Cereb10XEnsemblNames <- GeneLookup$EnsemID[match(unlist(Cereb10XGeneNames),GeneLookup$GeneID)]
#
# Get rid of Na rows (should be < 5%) - these didn't get converted
#
na_rows <-is.na(Cereb10XEnsemblNames)
rownames(rna_counts) <- Cereb10XEnsemblNames
dim(rna_counts)
rna_counts <- rna_counts[!na_rows, ]
dim(rna_counts)



# Create Seurat object
Cereb10X <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
#
# Need to do QC before changing rownames!
#
#Cereb10X[["percent.mt"]] <- PercentageFeatureSet(Cereb10X, pattern = "^MT-")

# Section 3 Now add in the ATAC-seq data
#
# we'll only use peaks in standard chromosomes
#library(Seurat)
#library(Signac)
#library(EnsDb.Hsapiens.v86)
#library(EnsDb.Hsapiens.v75)
library(dplyr)
library(ggplot2)
#
# Section for old v75 / hg19 (may need to remove)
#
#grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
#grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
#atac_counts <- atac_counts[as.vector(grange.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "hg19"

#frag.file <- "~/ProjectTest/10X Filtered PBMC/human_brain_3k_atac_fragments.tsv.gz"
#chrom_assay <- CreateChromatinAssay(
#  counts = atac_counts,
#  sep = c(":", "-"),
#  genome = 'hg19',
#  fragments = frag.file,
#  min.cells = 10,
#  annotation = annotations
#)
#
# New Section for v86 / hg38 - taken from 
#
atac_counts <- inputdata.10x$Peaks
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotationv86 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# ignore 24 warnings
seqlevels(annotationv86) <- paste0('chr', seqlevels(annotationv86))
seqlevelsStyle(annotationv86) <- 'UCSC'
genome(annotationv86) <- "hg38"

frag.file <- "~/ProjectTest/10X Filtered PBMC/human_brain_3k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotationv86
)




Cereb10X[["peaks"]] <- chrom_assay
Cereb10X
#
# Save file checkpoint
#

saveRDS(Cereb10X,file = "~/ProjectTest/Cereb10XincPeaks.rds")
#
# Restart
Cereb10X<-readRDS("~/ProjectTest/Cereb10XincPeaks.rds")
UpdateSeuratObject(Cereb10X)
Assays(Cereb10X)

#
# Check Default Assay and set to RNA
#
DefaultAssay(Cereb10X) <- "RNA"
DefaultAssay(Cereb10X)
#


# Visualize QC metrics as a violin plot
#VlnPlot(Cereb10X, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Cereb10X <- NormalizeData(Cereb10X, verbose = TRUE)
Cereb10X <- ScaleData(Cereb10X, verbose = TRUE)
Cereb10X <- FindVariableFeatures(Cereb10X, selection.method = "vst", nfeatures = 2000)
Cereb10X <- RunPCA(Cereb10X, npcs = 30, verbose = TRUE)
Cereb10X <- RunUMAP(Cereb10X, reduction = "pca", dims = 1:30, verbose = TRUE)
#
#
# Section to transfer anchors from Reference (Cereb) to query (Cereb10X) based on RNA seq data to create benchmark
#
Cereb.anchors <- FindTransferAnchors(reference = Cereb, query = Cereb10X,
                                        dims = 1:30, reference.reduction = "pca",verbose=TRUE)
#
# Could transfer different cluster info if required
#predictions <- TransferData(anchorset = Spc.anchors, refdata = Spc$supercluster_term,
#                            dims = 1:30)
predictions <- TransferData(anchorset = Cereb.anchors, refdata = Cereb$cell_type,
                            dims = 1:30)

Cereb10X <- AddMetaData(Cereb10X, metadata = predictions)
DimPlot(Cereb10X, reduction = "umap", group.by = "predicted.id")
saveRDS(Cereb10X,file = "~/ProjectTest/Cereb10XGroundTruth.rds")
write.csv(Cereb10X@meta.data,file="Cereb10Xmeta.csv")
Cereb10X<-readRDS(file="~/ProjectTest/Cereb10XGroundTruth.rds")
#
# ATAC code - see ATAC Roxy code
#
# We exclude the  first dimension as this is typically correlated with sequencing depth
DefaultAssay(Cereb10X)<-"peaks"
DefaultAssay(Cereb10X)
Cereb10Xpeaks<-DietSeurat(Cereb10X,assays = "peaks")
#
#
# optional QA step here
#
saveRDS(Cereb10X,file = "~/ProjectTest/Cereb10XPeaks.rds")
Cereb10Xpeaks<-readRDS(file="~/ProjectTest/Cereb10XPeaks.rds")
Cereb10Xpeaks<-UpdateSeuratObject(Cereb10Xpeaks)
Cereb10Xpeaks <- RunTFIDF(Cereb10Xpeaks)
Cereb10Xpeaks <- FindTopFeatures(Cereb10Xpeaks, min.cutoff = "q0")

#
Cereb10Xpeaks <- RunSVD(Cereb10Xpeaks)
Cereb10Xpeaks <- RunUMAP(Cereb10Xpeaks, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
Cereb10Xpeaks <- FindNeighbors(object = Cereb10Xpeaks, reduction = 'lsi', dims = 2:30)
Cereb10Xpeaks <- FindClusters(object = Cereb10Xpeaks, verbose = TRUE, algorithm = 3)

p2<-DimPlot(Cereb10Xpeaks,group.by = "seurat_clusters")
p2

p3<-DimPlot(Cereb10Xpeaks,group.by = "predicted.id")
p3
#
# Gene Activity (quite slow) maybe tune a bit and use gene-id....
#
gene.activities <- GeneActivity(Cereb10Xpeaks)
#
#
# Convert Rownames to ensemble format
#
Cereb10XGeneNames <- rownames(gene.activities)
#
# Convert Gene Names to Ensembl format
# Lookup file created from source with ~55k genes
#
Cereb10XEnsemblNames <- GeneLookup$EnsemID[match(unlist(Cereb10XGeneNames),GeneLookup$GeneID)]
#
# Get rid of Na rows (should be < 5%)
#
na_rows <-is.na(Cereb10XEnsemblNames)
rownames(gene.activities) <- Cereb10XEnsemblNames
dim(gene.activities)
gene.activities <-gene.activities[!na_rows, ]
dim(gene.activities)
# add the gene activity matrix to the Seurat object as a new synthetic RNA assay and normalize it
Cereb10Xpeaks[['RNA']] <- CreateAssayObject(counts = gene.activities)
Cereb10Xpeaks <- NormalizeData(
  object = Cereb10Xpeaks,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(Cereb10Xpeaks$nCount_RNA)
)
# Save the object as RDS
saveRDS(Cereb10Xpeaks,file = "~/ProjectTest/10X Filtered PBMC/10XPeaksPostAC2.rds")

#
# Now we determine cell types from the synthetic RNA assay
#
DefaultAssay(Cereb10Xpeaks) <- "RNA"
DefaultAssay(Cereb10Xpeaks)
#Cereb10Xpeaks <- NormalizeData(Cereb10Xpeaks, verbose = TRUE)
Cereb10Xpeaks <- ScaleData(Cereb10Xpeaks, verbose = TRUE)
Cereb10Xpeaks <- FindVariableFeatures(Cereb10Xpeaks, selection.method = "vst", nfeatures = 2000)
Cereb10Xpeaks <- RunPCA(Cereb10Xpeaks, npcs = 30, verbose = TRUE)
Cereb10Xpeaks <- RunUMAP(Cereb10Xpeaks, reduction = "pca", dims = 1:30, verbose = TRUE)
#
# Section to transfer anchors from Reference (Cereb) to query (Cereb10Xpeaks) based on synthetic RNA seq data to create results
#
#Cereb.anchors <- FindTransferAnchors(reference = Cereb, query = Cereb10Xpeaks,
#                                     dims = 1:30, reference.reduction = "pca",verbose=TRUE)
# different version taken from https://github.com/AprilYuge/ATAC-annotation-benchmark/blob/main/method_running/seurat3.r
# takes 10+minutes running CCA and long time filtering anchors
transfer.anchors <- FindTransferAnchors(reference = Cereb, query = Cereb10Xpeaks, 
                                        features = Cereb@assays$RNA@var.features,
                                         reduction = "cca",verbose=TRUE)
#
# Could transfer different cluster info if required. here we use cell_type

#
# New Version based on benchmarking paper
#


celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = Cereb$cell_type,
                                     weight.reduction = Cereb10Xpeaks[["lsi"]], dims = 2:30)
#
# Old Version below
#
#predictions <- TransferData(anchorset = Cereb.anchors, refdata = Cereb$cell_type,
#                            dims = 1:30)
# modified as data was sparse for V75 version
#predictions <- TransferData(anchorset = Cereb.anchors, refdata = Cereb$cell_type, k.weight=38,
#                            +                             dims = 1:30)
Cereb10Xpeaks <- AddMetaData(Cereb10Xpeaks, metadata = celltype.predictions)
DimPlot(Cereb10Xpeaks, reduction = "umap", group.by = "predicted.id")
saveRDS(Cereb10X,file = "~/ProjectTest/Cereb10XpeaksPredictions.rds")
write.csv(Cereb10Xpeaks@meta.data,file="Cereb10XpeaksPred.csv")
