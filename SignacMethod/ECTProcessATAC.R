#
#   Function to carry out pre-processing on ATAC seq data object
#
#
#   Input - RefObj = Seurat object containing reference data (inc labels)
#           atac_counts = array of ATAC counts per cell/gene
#           frag.file = reference to file location containing the fragment files and idx
#   Output - ATACObj = Seurat Object with ATAC assay and generated ATAC seq object 
#           including predictions of cell labels and some QC metrics


ECTProcessATAC<-function(atac_counts,RefObj,frag.file) {
  
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
  chrom_assay <- CreateChromatinAssay(
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
  #
  #
  ATACObj <- RunTFIDF(ATACObj)
  ATACObj <- FindTopFeatures(ATACObj, min.cutoff = "q0")
  
  #
  ATACObj <- RunSVD(ATACObj)
  ATACObj <- RunUMAP(ATACObj, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  ATACObj <- FindNeighbors(object = ATACObj, reduction = 'lsi', dims = 2:30)
  ATACObj <- FindClusters(object = ATACObj, verbose = TRUE, algorithm = 3)
  #
  # Gene Activity (quite slow) maybe tune a bit and use gene-id....
  #
  print("Gene Activity - Can be slow")
  gene.activities <- GeneActivity(ATACObj)
  #
  #
  # Convert Rownames to ensemble format
  #
  GeneLookup<-read.csv("/Users/tedrees/ProjectTest/GeneLookup.csv")
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
  ATACObj[['RNA']] <- CreateAssayObject(counts = gene.activities)
  ATACObj <- NormalizeData(
    object = ATACObj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(ATACObj$nCount_RNA)
  )
  # Save the object as RDS
  #saveRDS(ATACObj,file = "~/ProjectTest/10X Filtered PBMC/QueryObjPostAC2.rds")
  
  #
  # Now we determine cell types from the synthetic RNA assay
  #
  DefaultAssay(ATACObj) <- "RNA"
  DefaultAssay(ATACObj)
  #ATACObj <- NormalizeData(ATACObj, verbose = TRUE)
  ATACObj <- ScaleData(ATACObj, verbose = TRUE)
  ATACObj <- FindVariableFeatures(ATACObj, selection.method = "vst", nfeatures = 2000)
  ATACObj <- RunPCA(ATACObj, npcs = 30, verbose = TRUE)
  ATACObj <- RunUMAP(ATACObj, reduction = "pca", dims = 1:30, verbose = TRUE)
  #
  # Section to transfer anchors from Reference Object to query object based on synthetic RNA seq data to create results
  #
  # based on code from https://github.com/AprilYuge/ATAC-annotation-benchmark/blob/main/method_running/seurat3.r
  # takes 10+minutes running CCA and long time filtering anchors
  transfer.anchors <- FindTransferAnchors(reference = RefObj, query = ATACObj, 
                                          mapping.score.k=100,
                                          features = RefObj@assays$RNA@var.features,
                                          reduction = "cca",verbose=TRUE)
  #
  # Could transfer different cluster info if required. here we use cell_type

  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = RefObj$cell_type,
                                       weight.reduction = ATACObj[["lsi"]], dims = 2:30)
  #  
  ATACObj <- AddMetaData(ATACObj, metadata = celltype.predictions)
  #mapscores<-MappingScore(transfer.anchors,ndim=30)
  #ATACObj <- AddMetaData(ATACObj, metadata = mapscores)
  # may add this later - needs pca data unfortunately
  
}
