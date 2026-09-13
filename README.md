# scCertify

## Explainable Confidence Scoring for Single-Cell RNA-seq Annotations

`scCertify` is an R package for evaluating the reliability of single-cell RNA-seq annotations using explainable confidence scoring.

The framework integrates:

* UCell-based positive marker enrichment
* Negative marker consistency
* kNN neighborhood agreement
* Entropy-based uncertainty estimation
* Doublet-aware confidence modeling
* Ontology-aware label matching
* Confidence calibration
* Confidence classification
* Explainable confidence attribution
* Discovery-aware annotation analysis
* De novo cluster gene signatures
* Seurat and SingleCellExperiment support

`scCertify` is designed to provide an interpretable framework for identifying reliable, uncertain, contradictory, and potentially unusual cell annotations.

---

# Why scCertify?

Most single-cell annotation methods focus primarily on assigning a cell identity.

However, an assigned label does not necessarily indicate that the annotation is reliable.

Uncertainty can arise from:

* Weak marker enrichment
* Sparse transcriptomic profiles
* Transitional cellular states
* Heterogeneous cell populations
* Technical doublets
* Reference atlas mismatch
* Ambiguous neighborhood structure
* Batch effects
* Understudied or poorly characterized cell states

`scCertify` addresses this problem by quantifying annotation evidence and providing interpretable explanations for why an annotation may be reliable or uncertain.

---

# Core Concept

The main scCertify confidence framework combines multiple sources of evidence:

1. **Marker enrichment** — evaluates whether cells express genes supporting their predicted identity.
2. **Neighborhood agreement** — evaluates whether neighboring cells support the same annotation.
3. **Entropy uncertainty** — measures uncertainty in the annotation score distribution.
4. **Doublet probability** — penalizes confidence when cells show evidence of doublet contamination.

The resulting confidence score is accompanied by an interpretable confidence class and evidence-based explanations.

---

# Features

## Current Features

* Confidence scoring for single-cell annotations
* UCell-based marker enrichment scoring
* Negative marker consistency scoring
* kNN neighborhood agreement scoring
* Entropy-based uncertainty estimation
* Doublet-aware confidence scoring
* Ontology-aware label matching
* Confidence calibration
* Confidence classification
* Explainable confidence attribution
* Discovery-aware annotation status
* De novo cluster gene signatures
* Seurat integration
* SingleCellExperiment integration
* Publication-ready visualization support

---

# Standard Workflow

```text
Single-cell RNA-seq Data
            ↓
      Cell Annotation
        (e.g. SingleR)
            ↓
     Positive Marker
       Enrichment
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

# Discovery-Aware Workflow

For heterogeneous, understudied, or potentially unusual populations, `scCertify` provides an optional discovery-aware layer.

```text
Predicted Cell Identity
            ↓
      Positive Markers
            +
      Negative Markers
            ↓
      Evidence Analysis
            ↓
      Discovery Status
            ↓
   ┌──────────┼──────────────┐
   ↓          ↓              ↓
 Known    Possible Novel   Possible
          State            Transitional
                              State
```

The discovery-aware analysis is intended as an exploratory interpretation layer.

It **does not modify the main scCertify confidence score**.

---

# Installation

`scCertify` is currently available as a development package while development and Bioconductor integration continue.

## Install dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

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
    "Matrix",
    "ggplot2",
    "remotes"
))
```

## Install scCertify from GitHub

```r
remotes::install_github(
    "Jaya-Surya-dev/scCertify"
)
```

Once the package is available through Bioconductor, it can be installed using:

```r
BiocManager::install("scCertify")
```

---

# Quick Start

## Load Libraries

```r
library(scCertify)
library(Seurat)
library(SingleR)
library(celldex)
library(UCell)
library(scDblFinder)
```

## Load Example Dataset

```r
data("pbmc_small")
```

## Preprocess Data

```r
pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- RunPCA(pbmc_small)
pbmc_small <- RunUMAP(pbmc_small, dims = 1:10)
```

## Run SingleR Annotation

```r
sce <- as.SingleCellExperiment(pbmc_small)

ref <- HumanPrimaryCellAtlasData()

pred <- SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.main
)

pbmc_small$predicted_label <- pred$labels
```

## Detect Doublets

```r
sce <- scDblFinder(sce)

pbmc_small$doublet_score <-
    colData(sce)$scDblFinder.score

pbmc_small$doublet_class <-
    colData(sce)$scDblFinder.class
```

## Define Marker Database

```r
markers <- list(
    "B_cell" = c("MS4A1", "CD79A"),
    "T_cells" = c("CD3D", "IL7R"),
    "Monocyte" = c("LYZ", "S100A8"),
    "NK_cell" = c("NKG7", "GNLY"),
    "DC" = c("FCER1A", "CST3"),
    "Platelets" = c("PPBP", "PF4")
)
```

`scCertify` performs ontology-aware label matching to accommodate common differences in marker database and annotation naming conventions.

## Calculate Entropy

If annotation scores are available:

```r
pbmc_small$entropy_score <-
    entropy_score(pred$scores)
```

Entropy can subsequently be normalized:

```r
pbmc_small$entropy_norm <- (
    pbmc_small$entropy_score -
        min(pbmc_small$entropy_score)
) / (
    max(pbmc_small$entropy_score) -
        min(pbmc_small$entropy_score)
)
```

## Run scCertify

```r
pbmc_small <- cell_certify(
    pbmc_small,
    markers
)
```

The main certification framework calculates an interpretable confidence score using marker evidence, neighborhood agreement, entropy uncertainty, and doublet information.

---

# Negative Marker Analysis

Negative markers provide contradictory molecular evidence for a predicted cell identity.

For example, genes characteristic of B cells can be supplied as negative markers for cells predicted to be T cells.

```r
negative_markers <- list(
    "T_cell" = c("MS4A1", "CD79A"),
    "B_cell" = c("CD3D", "CD3E")
)
```

Calculate negative marker consistency:

```r
negative_scores <- negative_marker_score(
    object = pbmc_small,
    negative_markers = negative_markers
)
```

A higher score indicates stronger expression of genes that are inconsistent with the predicted identity.

Negative marker scores are **not incorporated into the main confidence score by default**.

---

# Discovery-Aware Annotation Status

The discovery-aware status combines positive marker evidence with contradictory negative marker evidence.

```r
status <- discovery_status(
    marker_score = positive_scores,
    negative_marker_score = negative_scores
)
```

The function provides four exploratory categories:

| Positive Evidence | Negative Evidence | Status |
|---|---|---|
| Strong | Weak | Known |
| Weak | Weak | Possible novel state |
| Strong | Strong | Possible transitional state |
| Weak | Strong | Insufficient evidence |

The default thresholds are:

```r
positive_threshold = 0.50
negative_threshold = 0.50
```

These thresholds are heuristic and intended for exploratory analysis rather than definitive biological classification.

---

# Combined Discovery-Aware Analysis

The complete discovery-aware analysis can be performed using:

```r
discovery_result <- discovery_aware(
    object = pbmc_small,
    markers = markers,
    negative_markers = negative_markers
)
```

The resulting data frame contains:

```text
positive_marker_score
negative_marker_score
discovery_status
```

Example:

```r
head(discovery_result)
```

The discovery-aware analysis provides an additional interpretive layer without changing the original scCertify confidence model.

---

# De Novo Cluster Gene Signatures

For heterogeneous or understudied datasets, `scCertify` can identify genes relatively enriched within existing clusters.

First, provide cluster identities in object metadata:

```r
pbmc_small$cluster <- Idents(pbmc_small)
```

Then calculate de novo cluster signatures:

```r
signatures <- de_novo_signatures(
    object = pbmc_small,
    cluster_column = "cluster",
    top_n = 10
)
```

The result is a named list containing the top enriched genes for each cluster.

For example:

```r
signatures$`0`
signatures$`1`
```

This functionality is intended to support exploratory investigation of unusual or heterogeneous populations.

It does not assign biological identities automatically.

---

# Seurat Support

`scCertify` supports Seurat objects for the main scoring and discovery-aware workflows.

```r
result <- cell_certify(
    pbmc_small,
    markers
)
```

Discovery-aware analysis:

```r
result <- discovery_aware(
    object = pbmc_small,
    markers = markers,
    negative_markers = negative_markers
)
```

De novo signatures:

```r
signatures <- de_novo_signatures(
    object = pbmc_small,
    cluster_column = "cluster",
    top_n = 10
)
```

---

# SingleCellExperiment Support

`scCertify` also supports SingleCellExperiment objects.

```r
result <- discovery_aware(
    object = sce,
    markers = markers,
    negative_markers = negative_markers
)
```

De novo cluster signatures can similarly be generated:

```r
signatures <- de_novo_signatures(
    object = sce,
    cluster_column = "cluster",
    top_n = 10
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

## Confidence Classes

```r
DimPlot(
    pbmc_small,
    group.by = "confidence_class"
)
```

![Confidence Classes](man/figures/confidence_classes.png)

## Entropy Landscape

```r
FeaturePlot(
    pbmc_small,
    features = "entropy_norm"
)
```

![Entropy Plot](man/figures/entropy_plot.png)

---

# Explain Confidence Attribution

For an individual cell:

```r
cell_id <- colnames(pbmc_small)[1]

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

The explanation identifies the major factors contributing to reduced confidence.

---

# Confidence Framework

The current scCertify confidence model integrates:

* Marker enrichment
* Neighborhood agreement
* Entropy certainty
* Doublet probability

The current scoring framework is:

```text
Confidence =
    0.35 × Marker Score
  + 0.35 × Neighborhood Agreement
  + 0.20 × Entropy Certainty
  - 0.10 × Doublet Probability
```

The main confidence model remains unchanged when discovery-aware functionality is used.

Discovery-aware analysis provides additional evidence and interpretation separately.

---

# Package Structure

```text
scCertify/

├── R/
│   ├── calibrate_confidence.R
│   ├── cell_certify.R
│   ├── classify_confidence.R
│   ├── de_novo_signatures.R
│   ├── discovery_aware.R
│   ├── discovery_status.R
│   ├── entropy_score.R
│   ├── explain_cell.R
│   ├── explain_confidence.R
│   ├── marker_score.R
│   ├── match_labels.R
│   ├── neighbor_score.R
│   └── negative_marker_score.R
│
├── tests/
│   └── testthat/
│
├── man/
├── DESCRIPTION
├── NAMESPACE
├── README.md
└── LICENSE
```

---

# Development Status

Current functionality includes:

* Explainable confidence scoring
* Positive marker evidence
* Negative marker evidence
* Neighborhood agreement
* Entropy-based uncertainty
* Doublet-aware scoring
* Confidence calibration
* Ontology-aware label matching
* Discovery-aware annotation status
* De novo cluster signatures
* Seurat support
* SingleCellExperiment support

The discovery-aware functionality is designed as an exploratory extension for datasets containing heterogeneous, understudied, or potentially unusual cell states.

---

# Planned Features

Future development may include:

* Out-of-distribution detection
* Broader lineage-level certification
* Novel and transitional state detection using additional evidence
* Trajectory-aware confidence scoring
* Multimodal confidence integration
* Spatial transcriptomics support
* Automatic marker retrieval
* Cell ontology integration
* Batch-aware confidence estimation
* Benchmarking framework
* Explainable AI visualization
* Atlas-scale optimization

---

# Citation

If you use `scCertify` in your work, please cite:

```text
Doddetipalli JS.
scCertify: Explainable confidence scoring for
single-cell RNA-seq annotations.
```

A formal publication citation will be added when the associated manuscript is published.

---

# Author

Jaya Surya Doddetipalli

---

# License

MIT License
