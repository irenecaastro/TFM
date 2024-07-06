library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(readxl)
options(future.globals.maxSize = 1e9)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)

#Load of the data: a sparse matrix obtained from tiny-sci-RNA-seq sequencing for each NE and 
#ETiX embryoid analyzed (matrix_ref), along with two dataframes: one containing cell annotations 
#(cell_ann_ref) and the other one containing gene annotations (gene_ann_ref). 

matrix_ref <- readMM("GSE209792_gene_count.mtx.gz")
cell_ann_ref <- read.csv("GSE209792_cell_annotate.csv.gz", text=TRUE, sep=',')
gene_ann_ref <- read.csv("GSE209792_gene_annotate.csv.gz", text=TRUE, sep=',')

#Assigning the gene names to the rows and cell names to the columns.

genes <- gene_ann_ref[,"gene_short_name"]
genes <- make.unique(genes)
row.names(matrix_ref) <- genes

cell <- cell_ann_ref[,"X"]
colnames(matrix_ref) <- cell

#Create a dataframe withe the metadata of interest
metadatos <- cell_ann_ref[, c("sample_resource", "stage", "condition", "sample_ID", "celltype")]
row.names(metadatos) <- cell

#Create a Seurat Object using the count matrix and the metadata dataframe. The parameter min.cells
#indicates the features included must be detected in at least one cell.
mouse_seurat <- CreateSeuratObject(counts = matrix_ref, meta.data = metadatos, min.cells=1)

##QUALITY_CONTROL####

#Visualization of sequencing depth (nCount_RNA) and number of features detected per cell (nFeature_RNA).
VlnPlot(mouse_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, same.y.lims = 3000)
hist(mouse_seurat$nCount_RNA)

#Filtering: unique feature counts above 200 and below 1,500. 
mouse_seurat<- subset(mouse_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 1500)

#Association of sample names to the identity of mouse_seurat Seurat Object.
Idents(mouse_seurat) <- mouse_seurat@meta.data$sample_ID

#Analysis


##NO_INTEGRATION####

#Make a subset with the samples of interest: early and later stages.
d1 <- subset(mouse_seurat, idents = c("NE-E7.5-1", "NE-E7.5-2", "NE-E8.5", "NE-E8.75", "ETiX-Day6-1", "ETiX-Day6-2", "ETiX-Day8-1", "ETiX-Day8-2"))

#The split function divides the ARN column into several groups according to the stage column.
d1[["RNA"]] <- split(d1[["RNA"]], f = d1$stage)

#Normalization. Scale factor of 10,000.
d1 <- NormalizeData(d1, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature selection. 2000 features selected.
d1 <- FindVariableFeatures(d1, selection.method = "vst", nfeatures = 2000)

#Scaling
d1 <- ScaleData(d1)

#Principal Component analysis and Elbowlpot
d1 <- RunPCA(d1)
ElbowPlot(d1, ndims = 100)

#Clustering. Use of first 50 PCs and a resolution parameter of 1.2.
d1 <- FindNeighbors(d1, dims = 1:50, reduction = "pca")
d1 <- FindClusters(d1, resolution = 1.2, cluster.name = "unintegrated_clusters")

#Non-linear dimensional reduction and visualization of UMAP plot
d1 <- RunUMAP(d1, dims=1:50)

DimPlot(d1, group.by = c("stage"))
DimPlot(d1, group.by = "celltype", label = TRUE)


##INTEGRATION####

#Integration using Harmony approach.
d1_integrated<- IntegrateLayers(object = d1, method = CCAIntegration, orig.reduction = "pca", new.reduction = "harmony",
                     verbose = FALSE)
d1_integrated[["RNA"]] <- JoinLayers(d1_integrated[["RNA"]])

#Clustering.
d1_integrated <- FindNeighbors(d1_integrated, reduction = "harmony", dims = 1:50)
d1_integrated <- FindClusters(d1_integrated, resolution = 1.2, cluster.name = "integrated_clusters")

#Non-dimensional reduction and visualization of UMAP plot.
d1_integrated <- RunUMAP(d1_integrated, dims = 1:50, reduction = "harmony",reduction.name = "umap.integrated")


DimPlot(d1_integrated, reduction = "umap.integrated", group.by = c("stage"))
DimPlot(d1_integrated, reduction = "umap.integrated", group.by = c("celltype"), label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", label = TRUE)

#Make a subset within sample_ID agreggating replicates.
D1_NE_75 <- subset(x = d1_integrated,subset = (sample_ID == c("NE-E7.5-1", "NE-E7.5-2")))
D1_cell_NE_75 <- colnames(D1_NE_75)
D1_Etix_D6 <- subset(x = d1_integrated,subset = (sample_ID == c("ETiX-Day6-1", "ETiX-Day6-2")))
D1_cell_Etix_D6 <- colnames(D1_Etix_D6)
D1_NE_85 <- subset(x = d1_integrated,subset = (sample_ID == c("NE-E8.5", "NE-E8.75")))
D1_cell_NE_85 <- colnames(D1_NE_85)
D1_Etix_D8 <- subset(x = d1_integrated,subset = (sample_ID == c("ETiX-Day8-1", "ETiX-Day8-2")))
D1_cell_Etix_D8 <- colnames(D1_Etix_D8)

#Paint in the UMAP plot only the sample of interest.
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_NE_75, group.by = "celltype", label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_NE_75, label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_Etix_D6, group.by = "celltype", label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_Etix_D6, label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_NE_85, group.by = "celltype", label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_NE_85, label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_Etix_D8, group.by = "celltype", label = TRUE)
DimPlot(d1_integrated, reduction = "umap.integrated", cells = D1_cell_Etix_D8, label = TRUE)



##CHARACTERIZATION OF EACH CLUSTER####

#Use function FindMarkers to find dentify markers of the cluster of interest by making a comparison with 
#the other two clusters.
Ex_Ectoderm_markers0 <- FindMarkers(d1_integrated, ident.1 = 0, ident.2 = c(5, 7))
Ex_Ectoderm_markers7 <- FindMarkers(d1_integrated, ident.1 = 7, ident.2 = c(0, 5))
Ex_Ectoderm_markers5 <- FindMarkers(d1_integrated, ident.1 = 5, ident.2 = c(7, 0))

#Selection of significant upregulated markers by using a cutoff threshold of a 0.05 p-adjusted value and a 2 log fold change.
Ex_Ectoderm_markers0 <- Ex_Ectoderm_markers0[Ex_Ectoderm_markers0$p_val_adj<0.05 & Ex_Ectoderm_markers0$avg_log2FC>2, ]
Ex_Ectoderm_markers7 <- Ex_Ectoderm_markers7[Ex_Ectoderm_markers7$p_val_adj<0.05 & Ex_Ectoderm_markers7$avg_log2FC>2, ]
Ex_Ectoderm_markers5 <- Ex_Ectoderm_markers5[Ex_Ectoderm_markers5$p_val_adj<0.05 & Ex_Ectoderm_markers5$avg_log2FC>2, ]

##GO ENRICHMENT ANALYSIS####

# Connect to Ensembl and specify specie: mouse dataset.
mouse_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

#Entrez nomenclature.

activated.genes.deseq2_Cluster0 <- toupper(rownames(Ex_Ectoderm_markers0))

mapping_activation_0 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_Cluster0,
  mart = mouse_mart
)

#Enrichment analysis of Cluster 0

activated.enrich.go_Cluster0 <- enrichGO(gene          = mapping_activation_0$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_Cluster0,showCategory = 17)

#Enrichment analysis of Cluster 7

activated.genes.deseq2_Cluster7 <- toupper(rownames(Ex_Ectoderm_markers7))

mapping_activation_7 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_Cluster7,
  mart = mouse_mart
)

activated.enrich.go_Cluster7 <- enrichGO(gene          = mapping_activation_7$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_Cluster7,showCategory = 17)

#Enrichment analysis of Cluster 5

activated.genes.deseq2_Cluster5 <- toupper(rownames(Ex_Ectoderm_markers5))

mapping_activation_5 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_Cluster5,
  mart = mouse_mart
)

activated.enrich.go_Cluster5 <- enrichGO(gene          = mapping_activation_5$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_Cluster5,showCategory = 17)


##PSEUDO-BULK DIFFERENTIAL EXPRESSION ANALYSIS####

##Subset the data to obtain objects of ExE population within stages.

D1_Stage_1 <- subset(x = d1_integrated,subset = (sample_ID == c("NE-E7.5-1", "NE-E7.5-2","ETiX-Day6-1", "ETiX-Day6-2")))
D1_cell_Stage_1 <- colnames(D1_Stage_1)
DESeq_Stage1_ExE <- subset(D1_Stage_1, idents = c(0, 5, 7))
Cell_ExE_Stage1 <- colnames(DESeq_Stage1_ExE)

Cluster_ID_1 <- Idents(DESeq_Stage1_ExE)

D1_Stage_2 <- subset(x = d1_integrated,subset = (sample_ID == c("NE-E8.5", "NE-E8.75", "ETiX-Day8-1", "ETiX-Day8-2")))
D1_cell_Stage_2 <- colnames(D1_Stage_2)
DESeq_Stage2_ExE <- subset(D1_Stage_2, idents = c(0, 5, 7))
Cell_ExE_Stage2 <- colnames(DESeq_Stage2_ExE)

Cluster_ID_2 <- Idents(DESeq_Stage2_ExE)

##Cluster 0####

#I will keep only early stages from cluster 0.

DESeq0_Stage1 <- subset(d1, idents = c(0))

Ex_Ectoderm_cells_0_Stage_1<- colnames(DESeq0_Stage1)

#I will keep only later stages from cluster 0.

DESeq0_Stage2 <- subset(D1_Stage_2, idents = c(0))

Ex_Ectoderm_cells_0_Stage_2<- colnames(DESeq0_Stage2)


##Cluster 7####

#I will keep only early stages from cluster 7.

DESeq7_Stage_1 <- subset(D1_Stage_1, idents = c(7))

Ex_Ectoderm_cells_7_Stage_1<- colnames(DESeq7_Stage_1)



#I will keep only later stages from cluster 7.

DESeq7_Stage_2 <- subset(D1_Stage_2, idents = c(7))

Ex_Ectoderm_cells_7_Stage_2<- colnames(DESeq7_Stage_2)



##Cluster 5####

#I will keep only early stages from cluster 5.

DESeq5_Stage_1 <- subset(D1_Stage_1, idents = c(5))

Ex_Ectoderm_cells_5_Stage_1<- colnames(DESeq5_Stage_1)


#I will keep only later stages from cluster 5.

DESeq5_Stage_2 <- subset(D1_Stage_2, idents = c(5))

Ex_Ectoderm_cells_5_Stage_2<- colnames(DESeq5_Stage_2)


##DESeq_STAGE1####
##Cluster 0####

metadata_Stage1 <- read.csv("metadata_Stage1.csv")
metadata_Stage1$sample_resource <- factor(metadata_Stage1$sample_resource)

#We use the Seurat object without integration.
DESeq0_Stage1 <- subset(mouse_seurat, cells = Ex_Ectoderm_cells_0_Stage_1)
# Aggregate across sample groups
bulk_0_Stage1 <- AggregateExpression(DESeq0_Stage1, slot = "counts", assays = "RNA", group.by = c("sample_ID","sample_resource"))

##DESeq2 object

countData_0_Stage1 <- as.matrix(bulk_0_Stage1$RNA)

dds_0_Stage1 <- DESeqDataSetFromMatrix(countData_0_Stage1, 
                                       colData = metadata_Stage1, 
                                       design = ~ sample_resource)


# Run DESeq2 differential expression analysis

dds_0_Stage1 <- DESeq(dds_0_Stage1)

# Plot dispersion estimates
plotDispEsts(dds_0_Stage1)

#Comparing in cluster 0 NE to ETiX
levels(metadata_Stage1$sample_resource)[2]
levels(metadata_Stage1$sample_resource)[1]

contrast_0_Stage1 <- c("sample_resource", levels(metadata_Stage1$sample_resource)[1],levels(metadata_Stage1$sample_resource)[2])

# resultsNames(dds)
res_0_Stage1 <- results(dds_0_Stage1,
                        contrast = contrast_0_Stage1,
                        alpha = 0.1)

gene.ids_0_Stage1 <- rownames(DESeq0_Stage1@assays$RNA$counts)


log.fold.change <- res_0_Stage1$log2FoldChange
q.value <- res_0_Stage1$padj
names(log.fold.change) <- gene.ids_0_Stage1
names(q.value) <- gene.ids_0_Stage1

activated.genes.deseq2_0_Stage1 <- gene.ids_0_Stage1[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2_0_Stage1 <- activated.genes.deseq2_0_Stage1[!is.na(activated.genes.deseq2_0_Stage1)]

repressed.genes.deseq2_0_Stage1 <- gene.ids_0_Stage1[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2_0_Stage1 <- repressed.genes.deseq2_0_Stage1[!is.na(repressed.genes.deseq2_0_Stage1)]

length(activated.genes.deseq2_0_Stage1)
length(repressed.genes.deseq2_0_Stage1)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 2 in either direction
res_table_thres_0_Stage1 <- (q.value < 0.05 & abs(log.fold.change) >= 2)

# Assuming res_table_thres, log.fold.change, and q.value are vectors of the same length
volcano_0_Stage1 <- data.frame(
  def = res_table_thres_0_Stage1,
  log2foldChange = log.fold.change,
  padjust = q.value
)

## Volcano plot
ggplot(volcano_0_Stage1) +
  geom_point(aes(x =log2foldChange, y = -log10(padjust), colour = def)) +
  ggtitle("Volcano plot of Extraembryonic ectoderm cells of cluster 0 relative to Natural Embryos E7.5 and ETiX Day 6") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


## Go enrichment analysis. 

# Connect to Ensembl and specify the specie: mouse dataset.
mouse_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

#Entrez nomenclature.

activated.genes.deseq2_0_Stage1 <- toupper(activated.genes.deseq2_0_Stage1)

mapping_activation_0_Stage1 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_0_Stage1,
  mart = mouse_mart
)

repressed.genes.deseq2_0_Stage1 <- toupper(repressed.genes.deseq2_0_Stage1)

mapping_repression_0_Stage1 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = repressed.genes.deseq2_0_Stage1,
  mart = mouse_mart
)

#getBM cannot find the Entrez ID of a gene so we introduce it manually.
mapping_repression_0_Stage1 <- rbind(mapping_repression_0_Stage1, 75745)


activated.enrich.go_0_Stage1 <- enrichGO(gene          = mapping_activation_0_Stage1$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_0_Stage1,showCategory = 17)


repressed.enrich.go_0_Stage1 <- enrichGO(gene          = mapping_repression_0_Stage1$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(repressed.enrich.go_0_Stage1, showCategory = 17)




##Cluster 7####

#We use the Seurat object without integration.
DESeq7_Stage_1 <- subset(mouse_seurat, cells = Ex_Ectoderm_cells_7_Stage_1)
# Aggregate across sample groups
bulk_7_Stage1 <- AggregateExpression(DESeq7_Stage_1, slot = "counts", assays = "RNA", group.by = c("sample_ID","sample_resource"))

##DESeq2 object

countData_7_Stage1 <- as.matrix(bulk_7_Stage1$RNA)

dds_7_Stage1 <- DESeqDataSetFromMatrix(countData_7_Stage1, 
                                       colData = metadata_Stage1, 
                                       design = ~ sample_resource)


# Run DESeq2 differential expression analysis

dds_7_Stage1 <- DESeq(dds_7_Stage1)

# Plot dispersion estimates
plotDispEsts(dds_7_Stage1)

#Comparing in cluster 7 NE to ETiX
levels(metadata_Stage1$sample_resource)[2]
levels(metadata_Stage1$sample_resource)[1]

contrast_7_Stage1 <- c("sample_resource", levels(metadata_Stage1$sample_resource)[1],levels(metadata_Stage1$sample_resource)[2])

# resultsNames(dds)
res_7_Stage1 <- results(dds_7_Stage1,
                        contrast = contrast_7_Stage1,
                        alpha = 0.1)

gene.ids_7_Stage1 <- rownames(DESeq7_Stage_1@assays$RNA$counts)


log.fold.change <- res_7_Stage1$log2FoldChange
q.value <- res_7_Stage1$padj
names(log.fold.change) <- gene.ids_7_Stage1
names(q.value) <- gene.ids_7_Stage1

activated.genes.deseq2_7_Stage1 <- gene.ids_7_Stage1[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2_7_Stage1 <- activated.genes.deseq2_7_Stage1[!is.na(activated.genes.deseq2_7_Stage1)]

repressed.genes.deseq2_7_Stage1 <- gene.ids_7_Stage1[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2_7_Stage1 <- repressed.genes.deseq2_7_Stage1[!is.na(repressed.genes.deseq2_7_Stage1)]

length(activated.genes.deseq2_7_Stage1)
length(repressed.genes.deseq2_7_Stage1)

#There are no DEGs

##Cluster 5####
metadata_5_Stage1 <- read.csv("metadata_5_Stage1.csv")
metadata_5_Stage1$sample_resource <- factor(metadata_5_Stage1$sample_resource)

#We use the Seurat object without integration.
DESeq5_Stage_1 <- subset(mouse_seurat, cells = Ex_Ectoderm_cells_5_Stage_1)

# Aggregate across sample groups
bulk_5_Stage1 <- AggregateExpression(DESeq5_Stage_1, slot = "counts", assays = "RNA", group.by = c("sample_ID","sample_resource"))

##DESeq2 object

countData_5_Stage1 <- as.matrix(bulk_5_Stage1$RNA)


dds_5_Stage1 <- DESeqDataSetFromMatrix(countData_5_Stage1, 
                                       colData = metadata_5_Stage1, 
                                       design = ~ sample_resource)


# Run DESeq2 differential expression analysis

dds_5_Stage1 <- DESeq(dds_5_Stage1)

# Plot dispersion estimates
plotDispEsts(dds_5_Stage1)

#Comparing in cluster 5 NE to ETiX
levels(metadata_5_Stage1$sample_resource)[2]
levels(metadata_5_Stage1$sample_resource)[1]

contrast_5_Stage1 <- c("sample_resource", levels(metadata_5_Stage1$sample_resource)[1],levels(metadata_5_Stage1$sample_resource)[2])

# resultsNames(dds)
res_5_Stage1 <- results(dds_5_Stage1,
                        contrast = contrast_5_Stage1,
                        alpha = 0.1)

gene.ids_5_Stage1 <- rownames(DESeq5_Stage_1@assays$RNA$counts)


log.fold.change <- res_5_Stage1$log2FoldChange
q.value <- res_5_Stage1$padj
names(log.fold.change) <- gene.ids_5_Stage1
names(q.value) <- gene.ids_5_Stage1

activated.genes.deseq2_5_Stage1 <- gene.ids_5_Stage1[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2_5_Stage1 <- activated.genes.deseq2_5_Stage1[!is.na(activated.genes.deseq2_5_Stage1)]

repressed.genes.deseq2_5_Stage1 <- gene.ids_5_Stage1[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2_5_Stage1 <- repressed.genes.deseq2_5_Stage1[!is.na(repressed.genes.deseq2_5_Stage1)]

length(activated.genes.deseq2_5_Stage1)
length(repressed.genes.deseq2_5_Stage1)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 2 in either direction
res_table_thres_5_Stage1 <- (q.value < 0.05 & abs(log.fold.change) >= 2)

# Assuming res_table_thres, log.fold.change, and q.value are vectors of the same length
volcano_5_Stage1 <- data.frame(
  def = res_table_thres_5_Stage1,
  log2foldChange = log.fold.change,
  padjust = q.value
)

## Volcano plot
ggplot(volcano_5_Stage1) +
  geom_point(aes(x =log2foldChange, y = -log10(padjust), colour = def)) +
  ggtitle("Volcano plot of Extraembryonic ectoderm cells of cluster 5 relative to Natural Embryos E7.5 and ETiX Day 6") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


## GO enrichment analysis. 

# Connect to Ensembl and specify the specie: mouse dataset.
mouse_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

#Entrez nomenclature.
activated.genes.deseq2_5_Stage1 <- toupper(activated.genes.deseq2_5_Stage1)

mapping_activation_5_Stage1 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_5_Stage1,
  mart = mouse_mart
)

repressed.genes.deseq2_5_Stage1 <- toupper(repressed.genes.deseq2_5_Stage1)

mapping_repression_5_Stage1 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = repressed.genes.deseq2_5_Stage1,
  mart = mouse_mart
)


activated.enrich.go_5_Stage1 <- enrichGO(gene          = mapping_activation_5_Stage1$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_5_Stage1,showCategory = 17)


repressed.enrich.go_5_Stage1 <- enrichGO(gene          = mapping_repression_5_Stage1$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(repressed.enrich.go_5_Stage1, showCategory = 17)


VlnPlot(d1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



##DESeq_STAGE 2####

##Cluster 0####
metadata_Stage2 <- read.csv("metadata_Stage2.csv")
metadata_Stage2$sample_resource <- factor(metadata_Stage2$sample_resource)

#We use the Seurat object without integration.
DESeq0_Stage2 <- subset(mouse_seurat, cells = Ex_Ectoderm_cells_0_Stage_2)

# Aggregate across sample groups
bulk_0_Stage2 <- AggregateExpression(DESeq0_Stage2, slot = "counts", assays = "RNA", group.by = c("sample_ID","sample_resource"))

##DESeq2 object

countData_0_Stage2 <- as.matrix(bulk_0_Stage2$RNA)

dds_0_Stage2 <- DESeqDataSetFromMatrix(countData_0_Stage2, 
                                       colData = metadata_Stage2, 
                                       design = ~ sample_resource)


# Run DESeq2 differential expression analysis

dds_0_Stage2 <- DESeq(dds_0_Stage2)

# Plot dispersion estimates
plotDispEsts(dds_0_Stage2)

#Comparing in cluster 0 NE to ETiX
levels(metadata_Stage2$sample_resource)[2]
levels(metadata_Stage2$sample_resource)[1]

contrast_0_Stage2 <- c("sample_resource", levels(metadata_Stage2$sample_resource)[1],levels(metadata_Stage2$sample_resource)[2])

# resultsNames(dds)
res_0_Stage2 <- results(dds_0_Stage2,
                        contrast = contrast_0_Stage2,
                        alpha = 0.1)

gene.ids_0_Stage2 <- rownames(DESeq0_Stage2@assays$RNA$counts)


log.fold.change <- res_0_Stage2$log2FoldChange
q.value <- res_0_Stage2$padj
names(log.fold.change) <- gene.ids_0_Stage2
names(q.value) <- gene.ids_0_Stage2

activated.genes.deseq2_0_Stage2 <- gene.ids_0_Stage2[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2_0_Stage2 <- activated.genes.deseq2_0_Stage2[!is.na(activated.genes.deseq2_0_Stage2)]

repressed.genes.deseq2_0_Stage2 <- gene.ids_0_Stage2[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2_0_Stage2 <- repressed.genes.deseq2_0_Stage2[!is.na(repressed.genes.deseq2_0_Stage2)]

length(activated.genes.deseq2_0_Stage2)
length(repressed.genes.deseq2_0_Stage2)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres_0_Stage2 <- (q.value < 0.05 & abs(log.fold.change) >= 2)

# Assuming res_table_thres, log.fold.change, and q.value are vectors of the same length
volcano_0_Stage2 <- data.frame(
  def = res_table_thres_0_Stage2,
  log2foldChange = log.fold.change,
  padjust = q.value
)

## Volcano plot
ggplot(volcano_0_Stage2) +
  geom_point(aes(x =log2foldChange, y = -log10(padjust), colour = def)) +
  ggtitle("Volcano plot of Extraembryonic ectoderm cells of cluster 0 relative to Natural Embryos E7.5 and ETiX Day 6") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


## GO enrichment analysis. 

# Connect to Ensembl and specify the specie: mouse dataset.
mouse_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

#Entrez nomenclature.
activated.genes.deseq2_0_Stage2 <- toupper(activated.genes.deseq2_0_Stage2)

mapping_activation_0_Stage2 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_0_Stage2,
  mart = mouse_mart
)

repressed.genes.deseq2_0_Stage2 <- toupper(repressed.genes.deseq2_0_Stage2)

mapping_repression_0_Stage2 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = repressed.genes.deseq2_0_Stage2,
  mart = mouse_mart
)

activated.enrich.go_0_Stage2 <- enrichGO(gene          = mapping_activation_0_Stage2$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_0_Stage2,showCategory = 17)


repressed.enrich.go_0_Stage2 <- enrichGO(gene          = mapping_repression_0_Stage2$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(repressed.enrich.go_0_Stage2, showCategory = 17)




##Cluster 7####

#We use the Seurat object without integration.
DESeq7_Stage_2 <- subset(mouse_seurat, cells = Ex_Ectoderm_cells_7_Stage_2)

# Aggregate across sample groups
bulk_7_Stage2 <- AggregateExpression(DESeq7_Stage_2, slot = "counts", assays = "RNA", group.by = c("sample_ID","sample_resource"))

##DESeq2 object

countData_7_Stage2 <- as.matrix(bulk_7_Stage2$RNA)

dds_7_Stage2 <- DESeqDataSetFromMatrix(countData_7_Stage2, 
                                       colData = metadata_Stage2, 
                                       design = ~ sample_resource)


# Run DESeq2 differential expression analysis

dds_7_Stage2 <- DESeq(dds_7_Stage2)

# Plot dispersion estimates
plotDispEsts(dds_7_Stage2)

#Comparing in cluster 7 NE to ETiX
levels(metadata_Stage2$sample_resource)[2]
levels(metadata_Stage2$sample_resource)[1]

contrast_7_Stage2 <- c("sample_resource", levels(metadata_Stage2$sample_resource)[1],levels(metadata_Stage2$sample_resource)[2])

# resultsNames(dds)
res_7_Stage2 <- results(dds_7_Stage2,
                        contrast = contrast_7_Stage2,
                        alpha = 0.1)

gene.ids_7_Stage2 <- rownames(DESeq7_Stage_2@assays$RNA$counts)


log.fold.change <- res_7_Stage2$log2FoldChange
q.value <- res_7_Stage2$padj
names(log.fold.change) <- gene.ids_7_Stage2
names(q.value) <- gene.ids_7_Stage2

activated.genes.deseq2_7_Stage2 <- gene.ids_7_Stage2[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2_7_Stage2 <- activated.genes.deseq2_7_Stage2[!is.na(activated.genes.deseq2_7_Stage2)]

repressed.genes.deseq2_7_Stage2 <- gene.ids_7_Stage2[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2_7_Stage2 <- repressed.genes.deseq2_7_Stage2[!is.na(repressed.genes.deseq2_7_Stage2)]

length(activated.genes.deseq2_7_Stage2)
length(repressed.genes.deseq2_7_Stage2)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 2 in either direction
res_table_thres_7_Stage2 <- (q.value < 0.05 & abs(log.fold.change) >= 2)

# Assuming res_table_thres, log.fold.change, and q.value are vectors of the same length
volcano_7_Stage2 <- data.frame(
  def = res_table_thres_7_Stage2,
  log2foldChange = log.fold.change,
  padjust = q.value
)

## Volcano plot
ggplot(volcano_7_Stage2) +
  geom_point(aes(x =log2foldChange, y = -log10(padjust), colour = def)) +
  ggtitle("Volcano plot of Extraembryonic ectoderm cells of cluster 5 relative to Natural Embryos E7.5 and ETiX Day 6") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


## GO enrichment analysis. 

# Connect to Ensembl and specify the specie: mouse dataset.
mouse_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

#Entrez nomenclature.
activated.genes.deseq2_7_Stage2 <- toupper(activated.genes.deseq2_7_Stage2)

mapping_activation_7_Stage2 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_7_Stage2,
  mart = mouse_mart
)

repressed.genes.deseq2_7_Stage2 <- toupper(repressed.genes.deseq2_7_Stage2)

mapping_repression_7_Stage2 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = repressed.genes.deseq2_7_Stage2,
  mart = mouse_mart
)

activated.enrich.go_7_Stage2 <- enrichGO(gene          = mapping_activation_7_Stage2$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_7_Stage2,showCategory = 17)


repressed.enrich.go_7_Stage2 <- enrichGO(gene          = mapping_repression_7_Stage2$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(repressed.enrich.go_7_Stage2, showCategory = 17)


##Cluster 5####
metadata_Stage2 <- read.csv("metadata_Stage2.csv")
metadata_Stage2$sample_resource <- factor(metadata_Stage2$sample_resource)

#We use the Seurat object without integration.
DESeq5_Stage_2 <- subset(mouse_seurat, cells = Ex_Ectoderm_cells_5_Stage_2)

# Aggregate across sample groups
bulk_5_Stage2 <- AggregateExpression(DESeq5_Stage_2, slot = "counts", assays = "RNA", group.by = c("sample_ID","sample_resource"))

##DESeq2 object

countData_5_Stage2 <- as.matrix(bulk_5_Stage2$RNA)


dds_5_Stage2 <- DESeqDataSetFromMatrix(countData_5_Stage2, 
                                       colData = metadata_Stage2, 
                                       design = ~ sample_resource)


# Run DESeq2 differential expression analysis

dds_5_Stage2 <- DESeq(dds_5_Stage2)

# Plot dispersion estimates
plotDispEsts(dds_5_Stage2)

#Comparing in cluster 5 NE to ETiX
levels(metadata_Stage2$sample_resource)[2]
levels(metadata_Stage2$sample_resource)[1]

contrast_5_Stage2 <- c("sample_resource", levels(metadata_Stage2$sample_resource)[1],levels(metadata_Stage2$sample_resource)[2])

# resultsNames(dds)
res_5_Stage2 <- results(dds_5_Stage2,
                        contrast = contrast_5_Stage2,
                        alpha = 0.1)

gene.ids_5_Stage2 <- rownames(DESeq5_Stage_2@assays$RNA$counts)


log.fold.change <- res_5_Stage2$log2FoldChange
q.value <- res_5_Stage2$padj
names(log.fold.change) <- gene.ids_5_Stage2
names(q.value) <- gene.ids_5_Stage2

activated.genes.deseq2_5_Stage2 <- gene.ids_5_Stage2[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2_5_Stage2 <- activated.genes.deseq2_5_Stage2[!is.na(activated.genes.deseq2_5_Stage2)]

repressed.genes.deseq2_5_Stage2 <- gene.ids_5_Stage2[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2_5_Stage2 <- repressed.genes.deseq2_5_Stage2[!is.na(repressed.genes.deseq2_5_Stage2)]

length(activated.genes.deseq2_5_Stage2)
length(repressed.genes.deseq2_5_Stage2)

## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 2 in either direction
res_table_thres_5_Stage2 <- (q.value < 0.05 & abs(log.fold.change) >= 2)

# Assuming res_table_thres, log.fold.change, and q.value are vectors of the same length
volcano_5_Stage2 <- data.frame(
  def = res_table_thres_5_Stage2,
  log2foldChange = log.fold.change,
  padjust = q.value
)

## Volcano plot
ggplot(volcano_5_Stage2) +
  geom_point(aes(x =log2foldChange, y = -log10(padjust), colour = def)) +
  ggtitle("Volcano plot of Extraembryonic ectoderm cells of cluster 5 relative to Natural Embryos E7.5 and ETiX Day 6") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,4.5)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


## GO enrichment analysis. 

# Connect to Ensembl and specify the specie: mouse dataset.
mouse_mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

#Entrez nomenclature.
activated.genes.deseq2_5_Stage2 <- toupper(activated.genes.deseq2_5_Stage2)

mapping_activation_5_Stage2 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = activated.genes.deseq2_5_Stage2,
  mart = mouse_mart
)

repressed.genes.deseq2_5_Stage2 <- toupper(repressed.genes.deseq2_5_Stage2)

mapping_repression_5_Stage2 <- getBM(
  attributes = "entrezgene_id",
  filters = "mgi_symbol",
  values = repressed.genes.deseq2_5_Stage2,
  mart = mouse_mart
)

activated.enrich.go_5_Stage2 <- enrichGO(gene          = mapping_activation_5_Stage2$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(activated.enrich.go_5_Stage2,showCategory = 17)


repressed.enrich.go_5_Stage2 <- enrichGO(gene          = mapping_repression_5_Stage2$entrezgene_id,
                                         OrgDb         = org.Mm.eg.db,
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         readable      = FALSE,
                                         keyType = "ENTREZID")

barplot(repressed.enrich.go_5_Stage2, showCategory = 17)


##Barplot_NE7.5_Day6_ExE####

DESeq_057_1 <- DESeq_Stage1_ExE

# Extract raw counts and metadata to create SingleCellExperiment object
counts_1 <- DESeq_057_1@assays$RNA@layers$counts

#Add to metadata the cluster ID
metadata_1 <- DESeq_057_1@meta.data
metadata_1$Cluster_ID <- Cluster_ID_1
metadata_1$sample_ID <- factor(metadata_1$sample_ID)

#Create Seurat object with the new metadata
DESeq_object_1 <- CreateSeuratObject(counts = counts_1, meta.data = metadata_1)

#Cell Count of each sample
conteo_0_NE7.5_1 <- length(which(DESeq_object_1@meta.data$sample_ID=="NE-E7.5-1" & DESeq_object_1@meta.data$Cluster_ID == "0"))
conteo_0_NE7.5_2 <- length(which(DESeq_object_1@meta.data$sample_ID=="NE-E7.5-2" & DESeq_object_1@meta.data$Cluster_ID == "0"))
conteo_0_D6_1 <- length(which(DESeq_object_1@meta.data$sample_ID=="ETiX-Day6-1" & DESeq_object_1@meta.data$Cluster_ID == "0"))
conteo_0_D6_2 <- length(which(DESeq_object_1@meta.data$sample_ID=="ETiX-Day6-2" & DESeq_object_1@meta.data$Cluster_ID == "0"))

conteo_5_NE7.5_1 <-length(which(DESeq_object_1@meta.data$sample_ID=="NE-E7.5-1" & DESeq_object_1@meta.data$Cluster_ID == "5"))
conteo_5_NE7.5_2 <- length(which(DESeq_object_1@meta.data$sample_ID=="NE-E7.5-2" & DESeq_object_1@meta.data$Cluster_ID == "5"))
conteo_5_D6_1 <- length(which(DESeq_object_1@meta.data$sample_ID=="ETiX-Day6-1" & DESeq_object_1@meta.data$Cluster_ID == "5"))
conteo_5_D6_2 <- length(which(DESeq_object_1@meta.data$sample_ID=="ETiX-Day6-2" & DESeq_object_1@meta.data$Cluster_ID == "5"))

conteo_7_NE7.5_1 <- length(which(DESeq_object_1@meta.data$sample_ID=="NE-E7.5-1" & DESeq_object_1@meta.data$Cluster_ID == "7"))
conteo_7_NE7.5_2 <- length(which(DESeq_object_1@meta.data$sample_ID=="NE-E7.5-2" & DESeq_object_1@meta.data$Cluster_ID == "7"))
conteo_7_D6_1 <- length(which(DESeq_object_1@meta.data$sample_ID=="ETiX-Day6-1" & DESeq_object_1@meta.data$Cluster_ID == "7"))
conteo_7_D6_2 <- length(which(DESeq_object_1@meta.data$sample_ID=="ETiX-Day6-2" & DESeq_object_1@meta.data$Cluster_ID == "7"))

#Cell proportion of replicates within clusters.

conteo_NE_7.5_0 <- numeric(2)
conteo_NE_7.5_0[1] <- conteo_0_NE7.5_1/sum(conteo_0_NE7.5_1, conteo_5_NE7.5_1, conteo_7_NE7.5_1)
conteo_NE_7.5_0[2] <- conteo_0_NE7.5_2/sum(conteo_0_NE7.5_2, conteo_5_NE7.5_2, conteo_7_NE7.5_2)

conteo_D6_0 <- numeric(2)
conteo_D6_0[1] <- conteo_0_D6_1/sum(conteo_0_D6_1, conteo_5_D6_1, conteo_7_D6_1)
conteo_D6_0[2] <- conteo_0_D6_2/sum(conteo_0_D6_2, conteo_5_D6_2, conteo_7_D6_2)


conteo_NE_7.5_5 <- numeric(2)
conteo_NE_7.5_5[1] <- conteo_5_NE7.5_1/sum(conteo_0_NE7.5_1, conteo_5_NE7.5_1, conteo_7_NE7.5_1)
conteo_NE_7.5_5[2] <- conteo_5_NE7.5_2/sum(conteo_0_NE7.5_2, conteo_5_NE7.5_2, conteo_7_NE7.5_2)


conteo_D6_5 <- numeric(2)
conteo_D6_5[1] <- conteo_5_D6_1/sum(conteo_0_D6_1, conteo_5_D6_1, conteo_5_D6_1)
conteo_D6_5[2] <- conteo_5_D6_2/sum(conteo_0_D6_2, conteo_5_D6_2, conteo_5_D6_2)


conteo_NE_7.5_7 <- numeric(2)
conteo_NE_7.5_7[1] <- conteo_7_NE7.5_1/sum(conteo_0_NE7.5_1, conteo_5_NE7.5_1, conteo_7_NE7.5_1)
conteo_NE_7.5_7[2] <- conteo_7_NE7.5_2/sum(conteo_0_NE7.5_2, conteo_5_NE7.5_2, conteo_7_NE7.5_2)

conteo_D6_7 <- numeric(2)
conteo_D6_7[1] <- conteo_7_D6_1/sum(conteo_0_D6_1, conteo_5_D6_1, conteo_7_D6_1)
conteo_D6_7[2] <- conteo_7_D6_2/sum(conteo_0_D6_2, conteo_5_D6_2, conteo_7_D6_2)

#Mean of cell proportion across replicates.

media_NE7.5_0 <- mean(conteo_NE_7.5_0)
media_NE7.5_5 <- mean(conteo_NE_7.5_5)
media_NE7.5_7 <- mean(conteo_NE_7.5_7)

media_D6_0 <- mean(conteo_D6_0)
media_D6_5 <- mean(conteo_D6_5)
media_D6_7 <- mean(conteo_D6_7)

cells_NE_7.5 <- data.frame(C0=media_NE7.5_0, C5=media_NE7.5_5, C7=media_NE7.5_7)

cells_D6 <- data.frame(C0=media_D6_0, C5=media_D6_5, C7=media_D6_7)

#Standard error of the mean (SEM)

S_NE7.5_0 <- sd(conteo_NE_7.5_0)/sqrt(2)
S_NE7.5_5 <- sd(conteo_NE_7.5_5)/sqrt(2)
S_NE7.5_7 <- sd(conteo_NE_7.5_7)/sqrt(2)

S_D6_0 <- sd(conteo_D6_0)/sqrt(2)
S_D6_5 <- sd(conteo_D6_5)/sqrt(2)
S_D6_7 <- sd(conteo_D6_7)/sqrt(2)

S_NE_7.5 <- data.frame(C0=S_NE7.5_0, C8=S_NE7.5_5, C9=S_NE7.5_7)
S_D6 <- data.frame(C0=S_D6_0, C5=S_D6_5, C7=S_D6_7)

#Global dataframe

mediafinal <-   data.frame(Type=rep(c("NE_7.5", "ETiX Day 6"), each=3), Cluster=c("Cluster 0", "Cluster 5", "Cluster 7"),
                           Cell_proportion=c(t(cells_NE_7.5), t(cells_D6)), sd=c(t(S_NE_7.5), t(S_D6)))
#Barplot

ggplot(mediafinal, aes(x=Cluster, y=Cell_proportion, fill=Type)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Cell_proportion-sd, ymax=Cell_proportion+sd), width=.2,
                position=position_dodge(.9))+labs(title = "Cell proportion of Extraembryonic ectoderm")

##Barplot_NE8.5_Day8_ExE####

DESeq_057_2 <- DESeq_Stage2_ExE

# Extract raw counts and metadata to create SingleCellExperiment object
counts_2 <- DESeq_057_2@assays$RNA@layers$counts

#Add to metadata the cluster ID
metadata_2 <- DESeq_057_2@meta.data
metadata_2$Cluster_ID <- Cluster_ID_2
metadata_2$sample_ID <- factor(metadata_2$sample_ID)

#Create Seurat object with the new metadata
DESeq_object_2 <- CreateSeuratObject(counts = counts_2, meta.data = metadata_2)

#Cell Count of each sample
conteo_0_NE8.5 <- length(which(DESeq_object_2@meta.data$sample_ID=="NE-E8.5" & DESeq_object_2@meta.data$Cluster_ID == "0"))
conteo_0_NE8.75 <- length(which(DESeq_object_2@meta.data$sample_ID=="NE-E8.75" & DESeq_object_2@meta.data$Cluster_ID == "0"))
conteo_0_D8_1 <- length(which(DESeq_object_2@meta.data$sample_ID=="ETiX-Day8-1" & DESeq_object_2@meta.data$Cluster_ID == "0"))
conteo_0_D8_2 <- length(which(DESeq_object_2@meta.data$sample_ID=="ETiX-Day8-2" & DESeq_object_2@meta.data$Cluster_ID == "0"))

conteo_5_NE8.5 <-length(which(DESeq_object_2@meta.data$sample_ID=="NE-E8.5" & DESeq_object_2@meta.data$Cluster_ID == "5"))
conteo_5_NE8.75 <- length(which(DESeq_object_2@meta.data$sample_ID=="NE-E8.75" & DESeq_object_2@meta.data$Cluster_ID == "5"))
conteo_5_D8_1 <- length(which(DESeq_object_2@meta.data$sample_ID=="ETiX-Day8-1" & DESeq_object_2@meta.data$Cluster_ID == "5"))
conteo_5_D8_2 <- length(which(DESeq_object_2@meta.data$sample_ID=="ETiX-Day8-2" & DESeq_object_2@meta.data$Cluster_ID == "5"))

conteo_7_NE8.5 <- length(which(DESeq_object_2@meta.data$sample_ID=="NE-E8.5" & DESeq_object_2@meta.data$Cluster_ID == "7"))
conteo_7_NE8.75 <- length(which(DESeq_object_2@meta.data$sample_ID=="NE-E8.75" & DESeq_object_2@meta.data$Cluster_ID == "7"))
conteo_7_D8_1 <- length(which(DESeq_object_2@meta.data$sample_ID=="ETiX-Day8-1" & DESeq_object_2@meta.data$Cluster_ID == "7"))
conteo_7_D8_2 <- length(which(DESeq_object_2@meta.data$sample_ID=="ETiX-Day8-2" & DESeq_object_2@meta.data$Cluster_ID == "7"))

#Cell proportion of replicates within clusters.

conteo_NE_8.5_0 <- numeric(2)
conteo_NE_8.5_0[1] <- conteo_0_NE8.5/sum(conteo_0_NE8.5, conteo_5_NE8.5, conteo_7_NE8.5)
conteo_NE_8.5_0[2] <- conteo_0_NE8.75/sum(conteo_0_NE8.75, conteo_5_NE8.75, conteo_7_NE8.75)

conteo_D8_0 <- numeric(2)
conteo_D8_0[1] <- conteo_0_D8_1/sum(conteo_0_D8_1, conteo_5_D8_1, conteo_7_D8_1)
conteo_D8_0[2] <- conteo_0_D8_2/sum(conteo_0_D8_2, conteo_5_D8_2, conteo_7_D8_2)


conteo_NE_8.5_5 <- numeric(2)
conteo_NE_8.5_5[1] <- conteo_5_NE8.5/sum(conteo_0_NE8.5, conteo_5_NE8.5, conteo_7_NE8.5)
conteo_NE_8.5_5[2] <- conteo_5_NE8.75/sum(conteo_0_NE8.75, conteo_5_NE8.75, conteo_7_NE8.75)


conteo_D8_5 <- numeric(2)
conteo_D8_5[1] <- conteo_5_D8_1/sum(conteo_0_D8_1, conteo_5_D8_1, conteo_5_D8_1)
conteo_D8_5[2] <- conteo_5_D8_2/sum(conteo_0_D8_2, conteo_5_D8_2, conteo_5_D8_2)


conteo_NE_8.5_7 <- numeric(2)
conteo_NE_8.5_7[1] <- conteo_7_NE8.5/sum(conteo_0_NE8.5, conteo_5_NE8.5, conteo_7_NE8.5)
conteo_NE_8.5_7[2] <- conteo_7_NE8.75/sum(conteo_0_NE8.75, conteo_5_NE8.75, conteo_7_NE8.75)

conteo_D8_7 <- numeric(2)
conteo_D8_7[1] <- conteo_7_D8_1/sum(conteo_0_D8_1, conteo_5_D8_1, conteo_7_D8_1)
conteo_D8_7[2] <- conteo_7_D8_2/sum(conteo_0_D8_2, conteo_5_D8_2, conteo_7_D8_2)

#Mean of cell proportion across replicates.

media_NE8.5_0 <- mean(conteo_NE_8.5_0)
media_NE8.5_5 <- mean(conteo_NE_8.5_5)
media_NE8.5_7 <- mean(conteo_NE_8.5_7)

media_D8_0 <- mean(conteo_D8_0)
media_D8_5 <- mean(conteo_D8_5)
media_D8_7 <- mean(conteo_D8_7)

cells_NE_8.5 <- data.frame(C0=media_NE8.5_0, C5=media_NE8.5_5, C7=media_NE8.5_7)

cells_D8 <- data.frame(C0=media_D8_0, C5=media_D8_5, C7=media_D8_7)

#Standard error of the mean (SEM)

S_NE8.5_0 <- sd(conteo_NE_8.5_0)/sqrt(2)
S_NE8.5_5 <- sd(conteo_NE_8.5_5)/sqrt(2)
S_NE8.5_7 <- sd(conteo_NE_8.5_7)/sqrt(2)

S_D8_0 <- sd(conteo_D8_0)/sqrt(2)
S_D8_5 <- sd(conteo_D8_5)/sqrt(2)
S_D8_7 <- sd(conteo_D8_7)/sqrt(2)

S_NE_8.5 <- data.frame(C0=S_NE8.5_0, C8=S_NE8.5_5, C9=S_NE8.5_7)
S_D8 <- data.frame(C0=S_D8_0, C5=S_D8_5, C7=S_D8_7)

#Global dataframe

mediafinal <-   data.frame(Type=rep(c("NE_8.5", "ETiX Day 8"), each=3), Cluster=c("Cluster 0", "Cluster 5", "Cluster 7"),
                           Cell_proportion=c(t(cells_NE_8.5), t(cells_D8)), sd=c(t(S_NE_8.5), t(S_D8)))
#Barplot

ggplot(mediafinal, aes(x=Cluster, y=Cell_proportion, fill=Type)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Cell_proportion-sd, ymax=Cell_proportion+sd), width=.2,
                position=position_dodge(.9))+labs(title = "Cell proportion of Extraembryonic ectoderm")