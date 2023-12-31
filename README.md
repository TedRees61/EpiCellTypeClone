# benchmark_epig_celltyping

This repository contains the code and details of the data used to investigate the classification of brain cell types from epigenetic data using machine learning. 
This was initially carried out as part of an M.Sc. Genomic Medicine Project at Imperial College London. 

The code provided here provides the materials used to carry out this benchmarking process as shown in the schematic below: -

![](EpiCellTypingSchematic2.png?raw=true)

High level Benchmark Process Outline
* a.	scRNA-seq and scATAC-seq from multi-omic benchmarking data processed separately
* b.	Standard scRNA-seq celltyping method applied to scRNA-seq benchmark component
* c.	Annotated reference data used for label transfer
* d.	Machine learning classifier uses scATAC-seq data to generate annotations for query data
* e.	Results from scRNA-seq and scATAC-seq annotations compared to assess performance


The overall process to generate both sets of results is covered in the code module
EpiCellTypBenchmark1.R in the SignacMethod directory, with different versions (2,3) for each reference dataset. This code is designed to be modified for different reference and benchmarking datasets as required. 
Please note it contains extra code to allow for saving and restarting since some steps are time consuming and may only need to be carried out once (for example pre-processing of reference datasets). 
These parts of the code can be included or excluded using #.
Please note the data paths will also need to be change from local references


For the **benchmark dataset (a)** we used the multi-omic dataset for Flash-frozen human healthy brain tissue (cerebellum) obtained from the 10XGenomics website available here:

https://www.10xgenomics.com/resources/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0

A copy of this data is on the RDS data store in this directory projects/neurogenomics-lab/live/Projects/benchmark_epig_celltyping/Data/10XGenomicsHumanBrainscMulti/

The **standard celltyping method (b)** was based on Seurat V4 described here:

https://satijalab.org/seurat/articles/integration_mapping

This is covered in the function ECTProcessRNABenchmark.R in the Signac Method directory
 
For the **scRNA-seq annotated reference (c)** we used 3 datasets related to the cerebellum from https://cellxgene.cziscience.com/datasets. These were derived from a [recently published paper](https://www.science.org/doi/10.1126/science.add7046) and are available as Seurat Objects on the CEllXGene data repository.

Copies of these have been stored in RDS here - projects/neurogenomics-lab/live/Projects/benchmark_epig_celltyping/Data/SilettiHumanBrain/

The **Machine Learning Classifier (d)** was based on the Signac method documented in the GeneActivity part of [this vignette](https://stuartlab.org/signac/articles/pbmc_vignette) and using the CCA version of the label transfer method described for the standard method for scRNA explained above.

This is contained in the ECTProcessATAC.R  Code 

The **performance analysis (e)** is carried out in EpiCTAnalysis1/2/3.R  for each of the reference datasets as above. This uses the ground truth (RNAObjmeta.csv) and ATACseq based predictions (ATACObjPred.csv) from steps (c) and (d) above. Summary results are generated for each cell type including confusion matrix, AuPR, AUROC and MCC and output in CellTypingResults1/2/3.csv This data can then be presented in tabular form or analysed in the Plotting.R code which generates suitable plots.

The dissertation paper is available in this repository in the dissertation pdf.

