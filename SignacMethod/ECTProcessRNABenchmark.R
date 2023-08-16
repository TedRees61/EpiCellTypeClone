#
#   Function to carry out processing for RNA query data object
#
#
#   Input - RefObj = Seurat object containing reference data (inc labels)
#           Rna_counts = array of raw RNA counts per cell/gene
#
#   Output - RNAObj = Seurat Object with RNA assay including predictions of cell labels and some QC metrics


ECTProcessRNABenchmark<-function(rna_counts,RefObj) {
  
  # Optional Code to create percent MT
  #Create Seurat object
  #RNAObj <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
  #
  # Populate QA metrics - pct_mitocondrial
  #  #
  #
  #RNAObj[["percent.mt"]] <- PercentageFeatureSet(RNAObj, pattern = "^MT-")
  #
  # Convert Gene Names to Esnembl ids
  #
  # Convert Rownames to ensemble format
  #
  GeneLookup<-read.csv("/Users/tedrees/ProjectTest/GeneLookup.csv")
  GeneNames <- rownames(rna_counts)
  #
  # Convert Gene Names to Ensembl format
  # Lookup file created from source with ~55k genes
  #
  EnsemblNames <- GeneLookup$EnsemID[match(unlist(GeneNames),GeneLookup$GeneID)]
  #
  # Get rid of Na rows (should be < 5%) - these didn't get converted
  #
  na_rows <-is.na(EnsemblNames)
  rownames(rna_counts) <- EnsemblNames
  dim(rna_counts)
  rna_counts <- rna_counts[!na_rows, ]
  dim(rna_counts)
  
  
    # Create Seurat object
  RNAObj <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
  #
  #
  # Review Quality metrics for possible analysis use
  #
  RNAObj <- NormalizeData(RNAObj, verbose = TRUE)
  RNAObj <- ScaleData(RNAObj, verbose = TRUE)
  RNAObj <- FindVariableFeatures(RNAObj, selection.method = "vst", nfeatures = 2000)
  RNAObj <- RunPCA(RNAObj, npcs = 30, verbose = TRUE)
  RNAObj <- RunUMAP(RNAObj, reduction = "pca", dims = 1:30, verbose = TRUE)
  #
  #
  # Section to transfer anchors from Reference (RefObj) to query (RNAObj) based on RNA seq data to create benchmark
  #
  rna.anchors <- FindTransferAnchors(reference = RefObj, query = RNAObj,
                                       dims = 1:30, reference.reduction = "pca",verbose=TRUE)
  #
  # Could transfer different cluster info if required
 
  predictions <- TransferData(anchorset = rna.anchors, refdata = RefObj$cell_type,
                              dims = 1:30)
  
  RNAObj <- AddMetaData(RNAObj, metadata = predictions)


  
}
