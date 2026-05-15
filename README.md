# cellCertR

## Confidence Scoring Framework for Single-Cell Cell Annotations

`cellCertR` is an R package designed to evaluate the reliability of single-cell annotations using:

* Marker gene consistency
* Neighborhood agreement
* Entropy-based uncertainty
* Ontology-aware label matching
* Explainable confidence scoring

The package aims to provide an interpretable confidence framework for scRNA-seq cell annotations.

---

# Features

## Current Features

* Automatic cell annotation using SingleR
* Marker-based confidence scoring
* kNN neighborhood agreement scoring
* Entropy-based uncertainty estimation
* Automatic ontology matching for inconsistent labels
* Explainable confidence interpretation
* Seurat integration
* UMAP confidence visualization

---

# Installation

## Install Dependencies

```r
install.packages(c(
  "Seurat",
  "ggplot2",
  "FNN",
  "entropy",
  "BiocManager"
))

BiocManager::install(c(
  "SingleR",
  "celldex",
  "SingleCellExperiment"
))
```

## Install Development Version

```r
install.packages("devtools")

devtools::install_github(
  "Jaya-Surya-dev/cellCertR"
)
```

---

# Workflow Example

## Load Libraries

```r
library(Seurat)
library(SingleR)
library(celldex)
library(cellCertR)
```

---

## Load Example Dataset

```r
data("pbmc_small")
```

---

## Preprocess Data

```r
pbmc_small <- NormalizeData(pbmc_small)

pbmc_small <- FindVariableFeatures(pbmc_small)

pbmc_small <- ScaleData(pbmc_small)

pbmc_small <- RunPCA(pbmc_small)

pbmc_small <- RunUMAP(
  pbmc_small,
  dims = 1:10
)
```

---

## Convert to SingleCellExperiment

```r
sce <- as.SingleCellExperiment(pbmc_small)
```

---

## Run SingleR Annotation

```r
ref <- HumanPrimaryCellAtlasData()

pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

pbmc_small$predicted_label <- pred$labels
```

---

## Define Marker Database

```r
markers <- list(

  "B_cell" = c(
    "MS4A1",
    "CD79A"
  ),

  "T_cells" = c(
    "CD3D",
    "IL7R"
  ),

  "Monocyte" = c(
    "LYZ",
    "S100A8"
  ),

  "NK_cell" = c(
    "NKG7",
    "GNLY"
  )
)
```

---

## Calculate Entropy

```r
pbmc_small$entropy_score <-
  entropy_score(
    pred$scores
  )

pbmc_small$entropy_norm <- (
  pbmc_small$entropy_score -
    min(pbmc_small$entropy_score)
) / (
  max(pbmc_small$entropy_score) -
    min(pbmc_small$entropy_score)
)
```

---

## Run cellCertR

```r
pbmc_small <- cell_certify(
  pbmc_small,
  markers
)
```

---

# Visualize Confidence

## Confidence UMAP

```r
FeaturePlot(
  pbmc_small,
  features = "confidence_score"
)
```

---

## Entropy Visualization

```r
FeaturePlot(
  pbmc_small,
  features = "entropy_norm"
)
```

---

## Confidence Distribution

```r
ggplot(
  pbmc_small@meta.data,
  aes(confidence_score)
) +
  geom_histogram(
    bins = 20
  )
```

---

# Explain Cell Confidence

```r
cell_id <- colnames(pbmc_small)[1]

explain_cell(
  pbmc_small,
  cell_id
)
```

Example output:

```txt
Cell: ATGCCAGAACGACT

Predicted Label: T_cells
Confidence Score: 0.82
Marker Score: 0.71
Neighbor Agreement: 0.90
Entropy: 0.12

Interpretation:
- High confidence annotation
```

---

# Current Methodology

The current confidence score combines:

* Marker consistency
* Neighborhood agreement
* Entropy uncertainty

Confidence formula:

```txt
Confidence = 0.4(Marker Score)
             + 0.4(Neighbor Agreement)
             + 0.2(1 - Entropy)
```

---

# Future Directions

Planned features include:

* UCell/AUCell marker enrichment
* Doublet-aware confidence
* Multimodal confidence scoring
* Spatial transcriptomics support
* CITE-seq integration
* ATAC-seq integration
* Trajectory-aware uncertainty
* Confidence calibration
* Automatic marker database retrieval
* Ontology integration
* Benchmarking framework

---

# Package Structure

```txt
cellCertR/
├── R/
│   ├── cell_certify.R
│   ├── entropy_score.R
│   ├── explain_cell.R
│   ├── marker_score.R
│   ├── match_labels.R
│   └── neighbor_score.R
│
├── man/
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

---

# Author

Jaya Surya Doddetipalli

---

# License

MIT License

