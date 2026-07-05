# scCertify

## Explainable Confidence Scoring for Single-Cell RNA-seq Annotations

`scCertify` is an R package for evaluating the reliability of single-cell RNA-seq annotations using explainable confidence scoring.

The framework integrates:

* UCell-based marker enrichment
* Neighborhood agreement scoring
* Entropy-based uncertainty estimation
* Doublet-aware confidence modeling
* Ontology-aware label matching
* Confidence calibration
* Explainable confidence attribution

`scCertify` aims to provide a biologically interpretable framework for identifying reliable and uncertain cell annotations in single-cell datasets.

---

# Why scCertify?

Most annotation tools assign labels without estimating how trustworthy those labels are.

In real single-cell datasets, uncertainty can arise from:

* Transitional cellular states
* Technical doublets
* Sparse transcriptomic profiles
* Weak marker enrichment
* Reference atlas mismatch
* Ambiguous neighborhood structure
* Batch effects

`scCertify` quantifies annotation reliability and explains why cells are considered uncertain.

---

# Features

## Current Features

* Confidence scoring for single-cell annotations
* UCell-based marker enrichment scoring
* kNN neighborhood agreement analysis
* Entropy-based uncertainty estimation
* Doublet-aware confidence scoring
* Ontology-aware label matching
* Confidence calibration
* Confidence classification
* Explainable confidence attribution
* Seurat integration
* Publication-ready visualization support

---

# Workflow Overview

```text id="k29d8s"
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

**scCertify** is currently under review for inclusion in the Bioconductor project. Until it becomes available through Bioconductor, the development version can be installed from GitHub.

## Install dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SingleCellExperiment",
    "SummarizedExperiment",
    "SingleR",
    "celldex",
    "UCell",
    "scDblFinder"
))

install.packages(c(
    "Seurat",
    "FNN",
    "entropy",
    "ggplot2",
    "remotes"
))
```

## Install scCertify

```r
remotes::install_github(
    "Jaya-Surya-dev/scCertify"
)
```

Once **scCertify** is accepted into Bioconductor, it can be installed using:

```r
BiocManager::install("scCertify")
```

---

# Quick Start

## Load Libraries

```r id="p82d7q"
library(scCertify)

library(Seurat)

library(SingleR)

library(celldex)

library(UCell)

library(scDblFinder)
```

---

## Load Example Dataset

```r id="w71x8p"
data("pbmc_small")
```

---

## Preprocess Data

```r id="m28d9w"
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

```r id="v74d8q"
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

```r id="j91d7x"
sce <- scDblFinder(sce)

pbmc_small$doublet_score <-
  colData(sce)$scDblFinder.score

pbmc_small$doublet_class <-
  colData(sce)$scDblFinder.class
```

---

## Define Marker Database

```r id="q17d8v"
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

```r id="z62d7r"
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

## Run scCertify

```r id="f81d9m"
pbmc_small <- cell_certify(
  pbmc_small,
  markers
)
```

---

# Example Outputs

## Confidence Score UMAP

```r id="g53d7x"
FeaturePlot(
  pbmc_small,
  features = "confidence_score"
)
```

![Confidence UMAP](man/figures/confidence_umap.png)

---

## Confidence Classes

```r id="r84d1q"
DimPlot(
  pbmc_small,
  group.by = "confidence_class"
)
```

![Confidence Classes](man/figures/confidence_classes.png)

---

## Entropy Landscape

```r id="d72x8p"
FeaturePlot(
  pbmc_small,
  features = "entropy_norm"
)
```

![Entropy Plot](man/figures/entropy_plot.png)

---

# Explain Confidence Attribution

```r id="x91d7v"
cell_id <- colnames(
  pbmc_small
)[1]

explain_confidence(
  pbmc_small,
  cell_id
)
```

Example output:

```text id="j37d8r"
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

```text id="u51d8p"
scCertify/

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

If you use `scCertify` in your work, please cite:

```text id="n84d7x"
Doddetipalli JS.
scCertify: Explainable confidence scoring for
single-cell RNA-seq annotations.
```

---

# Author

Jaya Surya Doddetipalli

---

# License

MIT License
