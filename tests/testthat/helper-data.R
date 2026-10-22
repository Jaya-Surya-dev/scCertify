library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(SummarizedExperiment)

set.seed(123)

##############################
## Simulated count matrix
##############################

counts <- matrix(
  rpois(20 * 10, lambda = 5),
  nrow = 20,
  ncol = 10
)

rownames(counts) <- paste0(
  "Gene",
  seq_len(20)
)

colnames(counts) <- paste0(
  "Cell",
  seq_len(10)
)

##############################
## Seurat object
##############################

toy_seurat <- SeuratObject::CreateSeuratObject(
  counts = counts
)

toy_seurat$predicted_label <- rep(
  c("T cell", "B cell"),
  each = 5
)

##############################
## Development cluster labels
##############################

toy_seurat$cluster <- c(
  rep("Cluster1", 5),
  rep("Cluster2", 5)
)

##############################
## Simple PCA embedding
##############################

pca_embedding <- matrix(
  rnorm(10 * 5),
  nrow = 10,
  ncol = 5
)

rownames(pca_embedding) <- colnames(toy_seurat)

colnames(pca_embedding) <- paste0(
  "PC_",
  1:5
)

toy_seurat[["pca"]] <- SeuratObject::CreateDimReducObject(
  embeddings = pca_embedding,
  key = "PC_",
  assay = SeuratObject::DefaultAssay(toy_seurat)
)

##############################
## SingleCellExperiment object
##############################

toy_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    counts = counts
  )
)

SummarizedExperiment::assay(
  toy_sce,
  "logcounts"
) <- log2(
  counts + 1
)

SummarizedExperiment::colData(
  toy_sce
)$predicted_label <- rep(
  c("T cell", "B cell"),
  each = 5
)

##############################
## Development cluster labels
##############################

SummarizedExperiment::colData(
  toy_sce
)$cluster <- c(
  rep("Cluster1", 5),
  rep("Cluster2", 5)
)

##############################
## PCA embedding
##############################

SingleCellExperiment::reducedDim(
  toy_sce,
  "PCA"
) <- pca_embedding

##############################
## Positive marker genes
##############################

toy_markers <- list(
  "T cell" = c(
    "Gene1",
    "Gene2",
    "Gene3"
  ),
  "B cell" = c(
    "Gene4",
    "Gene5",
    "Gene6"
  )
)

##############################
## Negative marker genes
##############################

toy_negative_markers <- list(
  "T cell" = c(
    "Gene4",
    "Gene5",
    "Gene6"
  ),
  "B cell" = c(
    "Gene1",
    "Gene2",
    "Gene3"
  )
)
