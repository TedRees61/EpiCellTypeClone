#
#   Function to carry out pre-processing on reference data object
#
#
#   Input - RefObj = Seurat object containing reference data (inc labels)

ECTProcessRefData<-function(RefObj) {

  #
  # Populate Counts - since this is not in the original file in the right place..:-(
  #
  if (nrow(RefObj@assays$RNA@counts)==0)
    
  {RefObj@assays$RNA@counts <- RefObj@assays$RNA@data}
  
  #
  # Run normalisation & Find Variable Features Step
  #
  
  RefObj <- NormalizeData(RefObj, verbose = TRUE)
  
  # Run the standard workflow for visualization and clustering
  RefObj <- FindVariableFeatures(RefObj, selection.method = "vst", nfeatures = 2000)
  RefObj <- ScaleData(RefObj, verbose = TRUE)
  
  RefObj <- RunPCA(RefObj, npcs = 30, verbose = TRUE)
  #
  # Warning! next step took ~10 hours
  #
  RefObj <- RunUMAP(RefObj, reduction = "pca", dims = 1:30, verbose = TRUE)
  
  
}