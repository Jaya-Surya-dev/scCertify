#' Discovery-Aware Annotation Analysis
#'
#' Combines positive and contradictory marker evidence to provide
#' a discovery-aware interpretation of predicted cell identities.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param markers Named list of positive marker genes.
#' @param negative_markers Named list of negative marker genes.
#' @param label_column Metadata column containing predicted labels.
#' @param positive_threshold Threshold for strong positive evidence.
#' @param negative_threshold Threshold for strong contradictory evidence.
#'
#' @return A data.frame containing positive marker scores,
#' negative marker scores, and discovery-aware annotation status.
#'
#' @details
#' This function provides an optional exploratory layer for
#' investigating heterogeneous or potentially unusual cell states.
#' It does not modify the main scCertify confidence score.
#'
#' The function combines positive marker evidence from
#' \code{marker_score()} with contradictory evidence from
#' \code{negative_marker_score()} and classifies each cell using
#' \code{discovery_status()}.
#'
#' @examples
#' set.seed(123)
#'
#' counts <- matrix(
#'   rpois(200, lambda = 5),
#'   nrow = 20,
#'   ncol = 10,
#'   dimnames = list(
#'     paste0("Gene", seq_len(20)),
#'     paste0("Cell", seq_len(10))
#'   )
#' )
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     logcounts = log1p(counts)
#'   )
#' )
#'
#' SummarizedExperiment::colData(sce)$predicted_label <- rep(
#'   c("T cell", "B cell"),
#'   each = 5
#' )
#'
#' markers <- list(
#'   "T cell" = c("Gene1", "Gene2"),
#'   "B cell" = c("Gene3", "Gene4")
#' )
#'
#' negative_markers <- list(
#'   "T cell" = c("Gene3", "Gene4"),
#'   "B cell" = c("Gene1", "Gene2")
#' )
#'
#' result <- discovery_aware(
#'   sce,
#'   markers = markers,
#'   negative_markers = negative_markers
#' )
#'
#' head(result)
#'
#' @export

discovery_aware <- function(
    object,
    markers,
    negative_markers,
    label_column = "predicted_label",
    positive_threshold = 0.50,
    negative_threshold = 0.50
) {

  positive_scores <- marker_score(
    object = object,
    markers = markers,
    label_column = label_column
  )

  negative_scores <- negative_marker_score(
    object = object,
    negative_markers = negative_markers,
    label_column = label_column
  )

  status <- discovery_status(
    marker_score = positive_scores,
    negative_marker_score = negative_scores,
    positive_threshold = positive_threshold,
    negative_threshold = negative_threshold
  )

  result <- data.frame(
    positive_marker_score = positive_scores,
    negative_marker_score = negative_scores,
    discovery_status = status,
    row.names = colnames(object),
    stringsAsFactors = FALSE
  )

  result
}
