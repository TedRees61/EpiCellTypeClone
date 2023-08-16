#
#     R program to run benchmarking for cell typing based on scATAC-seq data
#     ======================================================================
#
#     Required Input Data
#
#     RefData:      A set of RNA-seq data in SeuratObject format
#                   This data should be labelled by a reputatable source with Cell type data
#                   This data will be used to transfer labels
#
#     BenchmarkRNA: The scRNA-seq assay part of a 10X MultiOmic data set
#                   This will be used to transfer cell type labels 
#                   to create the ground truth cell type per barcode          
#
#     BenchmarkAtac:  The scATAC-seq part of a 10X MultiOmic data set
#                     This will be used to generate results data for the chosen method
#                     to create benchmark prediction to compare with the ground truth 
#                     cell type per barcode 
#
#
#     please note:  The ATAC seq fragment files and indices are also required!
#
#
#     Outputs:
#
#     groundTruth.csv   per barcode celltypes based on BenchmarkRNA with probabilities
#
#     predictions.csv   per barcode celltypes based on BenchmarkATAC with probabilities
#
#     Installations if required (optional)
#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
#remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
#remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
#remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
#remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
# install below seems to work better with some signac functions
#setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
#install.packages("Signac")#install.packages("tidyverse")
library(Seurat)
library(SeuratObject)
library(tidyverse)
#
library(ggplot2)
library(cowplot)
library(patchwork)
library(Signac)
#
#
# Functions Listed Below
#
source("ECTProcessRefData.R")
source("ECTProcessRNABenchmark.R")
source("ECTProcessATAC.R")
#
#   Part 1 - Process Reference File
#   -------------------------------
#
#   Read in Reference file (A Seurat Object)
#
#Cereb <- readRDS("/Users/tedrees/ProjectTest/Linardsson/CB Cerebellar Vermis.rds")
#Cereb<- readRDS("/Users/tedrees/ProjectTest/10X Filtered PBMC/Cereb2CbDN.rds")
Cereb<- readRDS("/Users/tedrees/ProjectTest/10X Filtered PBMC/CerebLatHem.rds")
#       Check its a Seurat object
Cereb
#       Carry out pre-processing
#
Cereb <- ECTProcessRefData(Cereb)
#
# Checkpoint Save & Read for restart
#
#saveRDS(Cereb,file = "~/ProjectTest/10X Filtered PBMC/CerebPostUMAP.rds")
saveRDS(Cereb,file = "~/ProjectTest/10X Filtered PBMC/Cereb3PostUMAP.rds")
#
# Restart point if needed to save time
#
Cereb <- readRDS("/Users/tedrees/ProjectTest/10X Filtered PBMC/Cereb2PostUMAP.rds")
Cereb <- UpdateSeuratObject(object = Cereb)
#
#   Part 2 - Read in MultiOmics Data for Benchmark 
#   ----------------------------------------------
#
# Read data (both RNA and ATAC) in from 10X genomics multiomics cellranger output
#
# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("~/ProjectTest/10X Filtered PBMC/human_brain_3k_filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

#
# location of ATAC-seq fragments
frag.file <- "~/ProjectTest/10X Filtered PBMC/human_brain_3k_atac_fragments.tsv.gz"

#
#   Part 3 - Process the RNA part to derive ground truth celltype labels per barcode
#
#   Input :   rna_counts matrix, labelled reference object
#   output:   seurat object with "RNA" assay including predicted cell ids from supplied refererence
RNAObj<-ECTProcessRNABenchmark(rna_counts,Cereb)
DimPlot(RNAObj, reduction = "umap", group.by = "predicted.id")
table(RNAObj$predicted.id)
saveRDS(RNAObj,file = "~/ProjectTest/RNAObjGroundTruth2.rds")
write.csv(RNAObj@meta.data,file="RNAObjmeta3.csv")  
#
#   Part 4 - Process the ATAC part to derive benchmark results celltype labels per barcode. 
#   This includes generation of synthetic RNA Assay
#
ATACObj<-ECTProcessATAC(atac_counts,Cereb,frag.file)
DimPlot(ATACObj, reduction = "umap", group.by = "predicted.id")
table(ATACObj$predicted.id)
saveRDS(ATACObj,file = "~/ProjectTest/ATACObjPredictions3.rds")
write.csv(ATACObj@meta.data,file="ATACObjPred3.csv")
# in case we need to build docker etc later
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

  