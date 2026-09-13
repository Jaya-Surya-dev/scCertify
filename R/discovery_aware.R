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
