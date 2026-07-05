#' scCertify: Explainable Confidence Scoring for Single-Cell Annotations
#'
#' scCertify provides an explainable framework for evaluating the confidence
#' of single-cell RNA sequencing (scRNA-seq) cell-type annotations.
#'
#' The package integrates multiple complementary sources of evidence,
#' including marker gene enrichment, neighborhood agreement,
#' entropy-based uncertainty, and doublet information to produce
#' calibrated confidence scores and interpretable confidence classes.
#'
#' Supported object classes include both Seurat and
#' SingleCellExperiment, enabling integration into existing
#' single-cell analysis workflows.
#'
#' Main functions include:
#'
#' \itemize{
#'   \item \code{\link{cell_certify}} — compute confidence scores.
#'   \item \code{\link{marker_score}} — evaluate marker enrichment.
#'   \item \code{\link{neighbor_score}} — compute neighborhood agreement.
#'   \item \code{\link{entropy_score}} — estimate annotation uncertainty.
#'   \item \code{\link{explain_cell}} — summarize confidence for a single cell.
#'   \item \code{\link{explain_confidence}} — explain confidence attribution.
#' }
#'
#' @docType package
#' @name scCertify
#' @keywords
