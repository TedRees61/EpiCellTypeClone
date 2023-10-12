#
#   Code to read in raw data from Azimuth motor cortex reference and create Seurat object
#
library(Seurat)
library(readr)
library(SeuratObject)
AziMatrix<-read.csv("~/ProjectTest/10X Filtered PBMC/AziMCmatrix.csv")
AziMeta<-read.csv("~/ProjectTest/10X Filtered PBMC/AziMCmeta.csv")
head(AziMatrix)
head(AziMeta)
AziSmall<-read.csv("~/ProjectTest/10X Filtered PBMC/AziSmall.csv")
dim(AziSmall)
AziRows<-AziSmall$sample_name
rownames(AziSmall)<-AziRows
AziRows[1:10]
max(AziRows)
rownames(AziSmall)
TAziSmall<-t(AziSmall)
AziSmeta<-read.csv("~/ProjectTest/10X Filtered PBMC/AziSmeta.csv")
data_numeric <- as.data.frame(lapply(AziSmall, as.numeric))
rownames(data_numeric)<-AziRows

> rownames(data_numeric)<-AziRows
> dim(data_numeric)
[1]   299 50282
> newdata<-data_numeric[,2:50282]
> View(newdata)
> tnewdata<-t(newdata)
CreateSeuratObject(counts = tnewdata, project = "ted_single_cell", min.cells = 3, min.features = 200)
