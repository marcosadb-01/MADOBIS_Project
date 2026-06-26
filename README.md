# MADOBIS_Project
Bioinformatic Identification of Novel Open Reading Frames in Differentially Expressed Unassigned Transcripts in Breast Cancer

Master's Final Project — Master's Degree in Omic Data Analysis and Systems Biology (MADOBIS)

Author: Marcos Antonio Dorantes Benitez



## Description

This repository contains the code developed for the Master's Thesis focused on the bioinformatic identification of novel ORFs (Open Reading Frames) in unannotated transcripts that are differentially expressed between healthy and tumour tissue in breast cancer (TCGA-BRCA).

The project combines two major blocks of work:


Transcriptomic and multi-omic analysis (TFM_DorantesBenitez.Rmd): processing of isoform-level expression data, exploratory analysis, differential expression, functional annotation, smORF identification, and regulatory network modelling.
Structural prediction of candidate peptides (AlphaFoldGenerico_DorantesBenitez.ipynb): prediction of the 3D structure of peptides encoded by the smORFs identified as relevant, using AlphaFold.



## Repository structure

MADOBIS_Project/

├── TFM_DorantesBenitez.Rmd                  # Complete pipeline in R (omic analysis)

├── AlphaFoldGenerico_DorantesBenitez.ipynb  # Structural prediction notebook (AlphaFold)

└── README.md                                # Repository description


Note: input data (TCGA-BRCA expression matrices, sample metadata, etc.) are not included in the repository due to their size. See the "Data and Code Availability Statement" section in the manuscript.



## TFM_DorantesBenitez.Rmd

R Markdown script implementing the complete analysis pipeline, divided into the following sections:

### 1.1. Data processing
Loading of Bioconductor/CRAN libraries (NOISeq, biomaRt, clusterProfiler, FactoMineR, mixOmics, MORE, igraph, etc.).
Loading and cleaning of the TCGA-BRCA isoform-level expression matrix (RNA-seq, RSEM normalised) and sample metadata.
Identification of vials per patient from TCGA barcodes and sample pairing (primary tumour vs. adjacent normal tissue from the same patient).

### 1.2. Exploratory analysis (PCA)
PCA of technical replicates (vials 01A/01B and 11A/11B) to rule out batch effects.
Global PCA (Normal / Primary tumour / Metastasis).
PCA of paired Tumour vs. Normal samples — basis for the differential expression analysis.

### 1.3. Differential expression
*Paired Wilcoxon test* (one-sample, H₀: median of the difference = 0) on the differential matrix (Healthy − Tumour) per patient, with filtering of isoforms with excess zeros, calculation of the contrast statistic (Z), and FDR p-value correction.
*NOISeq as a complementary method* (non-parametric, based on log-fold-change and probability of differential expression), including exploratory quality control (RNA composition, count distribution, sensitivity plot).
Comparison and intersection of significant isoforms between both methods (Venn diagrams).

### 1.4. Functional annotation and transcript characterisation
Query to Ensembl (biomaRt) and UCSC Genome Browser (RMySQL) to retrieve coordinates, biotype, and annotation of differentially expressed transcripts.
Identification of unannotated transcripts / lncRNAs with undescribed coding potential.
Query to the sORFs.org REST API for detection of small ORFs (smORFs) within those transcripts.
Functional enrichment analysis (GO / KEGG) with clusterProfiler and enrichplot.

### 1.5. Multi-omic modelling / regulatory networks
Regularised CCA (rCCA) with mixOmics to evaluate the correlation between smORF expression and protein-coding genes.
MORE (PLS1 with variable selection by permutation) to model regulator–target relationships and build the smORF–gene network (networkMORE).
Export of correlation networks to .gml format with igraph for visualisation in Cytoscape.
Intersection of nodes between the MORE and rCCA networks using Venn diagrams.

### Output
The result of this script is a list of candidate smORFs (unannotated transcripts, differentially expressed between tumour and healthy tissue, with coding evidence and regulatory relationships with known genes), which constitute the input for the next block.


## AlphaFoldGenerico_DorantesBenitez.ipynb

Generic notebook for the prediction of the three-dimensional structure of peptides encoded by the candidate smORFs selected in the previous analysis, using AlphaFold.

Takes as input the peptide sequences derived from the smORFs of interest.
Generates 3D structure predictions and associated confidence metrics (e.g. pLDDT).
Enables a first structural characterisation of peptides with no known homology, as additional evidence of their potential biological function.

Designed to run on Google Colab (recommended, due to GPU requirements) or on a local environment with a compatible GPU.



## Requirements

### For TFM_DorantesBenitez.Rmd

R (≥ 4.2 recommended) and RStudio.

CRAN packages: dplyr, ggplot2, ggplotify, ggrepel, ggvenn, VennDiagram, FactoMineR, factoextra, patchwork, reshape2, stringr, httr, jsonlite, igraph, RMySQL.

Bioconductor packages: NOISeq, biomaRt, clusterProfiler, enrichplot, GO.db, org.Hs.eg.db, AnnotationDbi, mixOmics, MORE, AnnotationDbi.

Quick installation:

```
install.packages(c("dplyr", "ggplot2", "ggplotify", "ggrepel", "ggvenn", 
                   "VennDiagram", "FactoMineR", "factoextra", "patchwork", 
                   "reshape2", "stringr", "httr", "jsonlite", "igraph", "RMySQL"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("NOISeq", "biomaRt", "clusterProfiler", "enrichplot", 
                       "GO.db", "org.Hs.eg.db", "AnnotationDbi", "mixOmics", "MORE"))
```                   

### For AlphaFoldGenerico_DorantesBenitez.ipynb

Python ≥ 3.8

Jupyter Notebook / Google Colab

AlphaFold (or ColabFold) and its associated dependencies (biopython, numpy, etc. depending on the implementation used in the notebook)



## Input data

The R script assumes the presence of the following files in a local folder (not included in the repository due to size/licence):

- `BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt` TCGA-BRCA isoform expression matrix (level 3, RSEM normalised).
- `sample_type.csv` metadata with sample type (primary tumour / normal / metastasis) by TCGA barcode.
- `Galaxy_multiBigwigSummary_bin_counts.tabular` tabular file generated from multiBigwigSummary in GALAXY (Version 3.5.4+galaxy0) using the files generated by the pipeline itself `smORFs_hg19.bed` and the BigWig tracks for hg19 (`PhyloCSF+0.bw, +1.bw, +2.bw`) from the Broad Institute. 
- `filtered_sorfs.csv` list of ORFs from the database https://sorfs.org/database, following the instructions in the manuscript.



## Usage

Download the input data and place them in the working directory or in a data/ folder, following the path instructions in the pipeline itself.
Open `TFM_DorantesBenitez.Rmd` in RStudio and run the chunks sequentially (the longer sections allow saving/loading intermediate objects in .RData to avoid recomputing costly steps).
Export the list of candidate smORFs and their peptide sequences.
Load those sequences into `AlphaFoldGenerico_DorantesBenitez.ipynb` (ideally in Google Colab) to obtain the structural predictions.



## Academic context

Master's Final Project of the Master's Degree in Omic Data Analysis and Systems Biology (MADOBIS).



## Licence

This repository is provided for academic purposes. If you wish to reuse the code, please cite the corresponding Thesis.
