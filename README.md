# TFM
# Introduction
In this repository, we will work in an R environment using RStudio for the analysis of sc-RNA-seq data following Seurat Guided Clustering Tutorial. Data can be accessed on Gene Expression Omnibus (GEO) database under the accession number GSE209792 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209792).
# Pipeline Description
- Step 1: From two cell and gene annotation dataframes and one sparse-matrix, we create a Seurat object. 
- Step 2: Cell filtering.
- Step 3: Normalization, feature selection, scaling, PCA, Clustering and non-linar dimensional reduction (UMAP).
- Step 4: Integration of the data due to presence of bactch effect.
- Step 5: Reclustering and cluster characterization of ExE. GO enrichment analysis of the cluster markers.
- Step 6: Cell proportion of early and later stages across clusters in ExE.
- Step 7: Pseudo-bulk differential expression analysis using NE as a control in ExE clusters. GO enrichment analysis of DEGs.
# Autor
Irene Carrero Castro
# Citation
Amadei, G., Handford, C. E., Qiu, C., De Jonghe, J., Greenfeld, H., Tran, M., et al. (2022b). Embryo model completes gastrulation to neurulation and organogenesis. Nature 610, 143–153. doi: 10.1038/s41586-022-05246-3
# Technical Support
If you have any questions, we recommend commenting on the 'issue' tab for resolution. If the issue is not resolved this way, please do not hesitate to send an email with the incident to irenecarrerocastro@gmail.com.
