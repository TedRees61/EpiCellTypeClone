counts<-GetAssayData(Cereb, layer="data")
meta<-Cereb@meta.data
NewCereb <- CreateSeuratObject(
  counts,
  project = "SeuratProject",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  names.field = 1,
  names.delim = "_",
  meta.data = meta
)
