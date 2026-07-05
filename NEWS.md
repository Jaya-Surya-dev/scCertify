# scCertify 0.99.2

## Changes in response to Bioconductor review

### New features

* Added support for both Seurat and SingleCellExperiment objects across the package.
* Expanded unit test coverage for exported functions.
* Added package-level documentation (`?scCertify`).
* Added executable examples for exported functions.
* Expanded the vignette with an introduction, workflow overview, methodology, and interpretation of confidence scores.
* Improved compatibility with Bioconductor documentation standards.
* Updated README installation instructions.
* Improved documentation throughout the package.

### Bug fixes

* Improved handling of metadata for both supported object classes.
* Improved dimensionality reduction handling for SingleCellExperiment objects.
* Improved marker score computation for SingleCellExperiment inputs.
* Improved neighborhood agreement compatibility across supported object types.

### Documentation

* Updated function documentation for Seurat and SingleCellExperiment compatibility.
* Added additional usage examples.
* Revised vignette formatting using BiocStyle.

