# scCertify

## Explainable Confidence Scoring for Single-Cell RNA-seq Cell Annotations

`scCertify` is an R package for evaluating the reliability of single-cell RNA-seq cell annotations using explainable confidence scoring.

The framework integrates:

* Marker enrichment scoring (UCell)
* Neighborhood agreement analysis
* Entropy-based uncertainty estimation
* Doublet-aware confidence scoring
* Ontology-aware label matching
* Confidence calibration
* Explainable confidence attribution

The goal of `scCertify` is to provide a biologically interpretable framework for assessing annotation trustworthiness in single-cell datasets.

---

# Why cellCertR?

Cell annotation tools often provide labels without estimating how reliable those annotations are.

In real biological datasets, annotation uncertainty can arise from:

* Transitional cellular states
* Technical doublets
* Sparse transcriptomic profiles
* Batch effects
* Reference atlas mismatch
* Weak marker expression
* Ambiguous neighborhood structure

`scCertify` addresses this by quantifying annotation confidence and explaining why specific cells are uncertain.

---

# Features

## Current Features

* Single-cell annotation confidence scoring
* UCell-based marker enrichment
* kNN neighborhood agreement scoring
* Entropy-based uncertainty estimation
* Doublet-aware confidence modeling
* Ontology-aware label matching
* Confidence calibration
* Confidence classification
* Explainable confidence attribution
* Seurat integration
* Publication-ready visualization support

---

# Workflow Overview

```text
Single-cell RNA-seq Data
            ↓
      SingleR Annotation
            ↓
     Marker Enrichment
         (UCell)
            ↓
   Neighborhood Agreement
            ↓
    Entropy Uncertainty
            ↓
     Doublet Detection
      (scDblFinder)
            ↓
   Confidence Calibration
            ↓
 Explainable Confidence
            ↓
 Final Confidence Score
   + Confidence Class
```

---

# Installation

## Install Dependencies

```r
install.packages(c(
  "Seurat",
  "ggplot2",
  "FNN",
  "entropy",
  "devtools"
))

install.packages("BiocManager")

BiocManager::install(c(
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "UCell",
  "scDblFinder"
))
```

---

## Install cellCertR

```r
devtools::install_github(
  "Jaya-Surya-dev/scCertify"
)
```

---

# Quick Start

## Load Libraries

```r
library(cellCertR)

library(Seurat)

library(SingleR)

library(celldex)

library(UCell)

library(scDblFinder)
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

## Run SingleR Annotation

```r
sce <- as.SingleCellExperiment(
  pbmc_small
)

ref <- HumanPrimaryCellAtlasData()

pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

pbmc_small$predicted_label <- pred$labels
```

---

## Detect Doublets

```r
sce <- scDblFinder(sce)

pbmc_small$doublet_score <-
  colData(sce)$scDblFinder.score

pbmc_small$doublet_class <-
  colData(sce)$scDblFinder.class
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
  ),

  "DC" = c(
    "FCER1A",
    "CST3"
  ),

  "Platelets" = c(
    "PPBP",
    "PF4"
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

# Example Outputs

## Confidence Score UMAP

```r
FeaturePlot(
  pbmc_small,
  features = "confidence_score"
)
```

![Confidence UMAP](man/figures/confidence_umap.png)

---

## Confidence Classes

```r
DimPlot(
  pbmc_small,
  group.by = "confidence_class"
)
```

![Confidence Classes](man/figures/confidence_classes.png)

---

## Entropy Landscape

```r
FeaturePlot(
  pbmc_small,
  features = "entropy_norm"
)
```

![Entropy Plot](man/figures/entropy_plot.png)

---

# Explain Cell Confidence

```r
cell_id <- colnames(
  pbmc_small
)[1]

explain_confidence(
  pbmc_small,
  cell_id
)
```

Example output:

```text
[1] "Weak marker enrichment"
[2] "High annotation uncertainty"
[3] "Possible doublet contamination"
```

---

# Confidence Framework

The current confidence model integrates:

* Marker enrichment
* Neighborhood agreement
* Entropy certainty
* Doublet probability

Current scoring framework:

Confidence =
0.35(Marker Score) +
0.35(Neighborhood Agreement) +
0.20(Entropy Certainty) -
0.10(Doublet Probability)

---

# Package Structure

```text
cellCertR/

├── R/
│   ├── calibrate_confidence.R
│   ├── cell_certify.R
│   ├── classify_confidence.R
│   ├── entropy_score.R
│   ├── explain_cell.R
│   ├── explain_confidence.R
│   ├── marker_score.R
│   ├── match_labels.R
│   └── neighbor_score.R
│
├── man/
├── DESCRIPTION
├── NAMESPACE
├── README.md
└── LICENSE
```

---

# Planned Features

* Trajectory-aware confidence scoring
* Multimodal confidence integration
* Spatial transcriptomics support
* Automatic marker retrieval
* Cell ontology integration
* Batch-aware confidence estimation
* Calibration models
* Benchmarking framework
* Explainable AI visualization
* Atlas-scale optimization

---

# Citation

If you use `cellCertR` in your work, please cite:

```text
Doddetipalli JS.
cellCertR: Explainable confidence scoring for
single-cell RNA-seq annotations.
```

---

# Author

Jaya Surya Doddetipalli

---

# License

MIT License
