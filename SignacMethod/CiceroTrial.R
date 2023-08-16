#
# Code to try out Cicero
#
# Install Cicero
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
library(cicero)
library(monocle3)
# convert to CellDataSet format and make the cicero object
ATAC.cds <- as.CellDataSet(x = ATACObj)
ATAC.cicero <- make_cicero_cds(ATAC.cds, reduced_coordinates = reducedDims(ATAC.cds)$UMAP)
#
# A different attempt
#
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
#
# Now create Seurat/Signac chromatin assay
#
CICchrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotationv86
)
CICATACObj <- CreateSeuratObject(counts = CICchrom_assay, assay = "ATAC") 
CICATACObj
#
# A 2nd different attempt
#
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg38)
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c("_", "_"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotationv86 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# ignore 24 warnings
seqlevels(annotationv86) <- paste0('chr', seqlevels(annotationv86))
seqlevelsStyle(annotationv86) <- 'UCSC'
genome(annotationv86) <- "hg38"
#
# Now create Seurat/Signac chromatin assay
#
CICchrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotationv86
)
ATACObj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks") 
ATACObj
#
# From GPT4 suggestions
#
# Assuming your data is in a directory called 'cellranger_output_directory'

# Read in the fragments file
fragments <- read.table(paste0("~/ProjectTest/", "10X Filtered PBMC/human_brain_3k_atac_fragments.tsv.gz"))

# Read in the peak-barcode matrix
peaks <- readMM(paste0("~/ProjectTest/", "10X Filtered PBMC/human_brain_3k_filtered_feature_bc_matrix.h5"))
rownames(peaks) <- readLines(paste0(cellranger_output_directory, "/filtered_peak_bc_matrix/peaks.bed"))
colnames(peaks) <- readLines(paste0(cellranger_output_directory, "/filtered_peak_bc_matrix/barcodes.tsv.gz"))
#
# From Cicero Website https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#loading-10x-scatac-seq-data
#
# read in matrix data using the Matrix package
indata <- Matrix::readMM("/Users/tedrees/ProjectTest/10X Filtered PBMC/matrix.mtx") 
# binarize the matrix
indata@x[indata@x > 0] <- 1

# format cell info
cellinfo <- read.csv("/Users/tedrees/ProjectTest/10X Filtered PBMC/human_brain_3k_per_barcode_metrics.csv")
row.names(cellinfo) <- cellinfo$atac_barcode
names(cellinfo) <- "cells"

# format peak info
peakinfo <- read.table("/Users/tedrees/ProjectTest/10X Filtered PBMC/human_brain_3k_atac_peaks.bed")
#
# may need to deal with dodgy data in here
#
peakinfo <- peakinfo[substr(peakinfo$chr, 1, 3) == "chr", ]
names(peakinfo) <- c("chr", "bp1", "bp2")

peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 