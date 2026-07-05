#' Explain Confidence Attribution
#'
#' Explains why a cell received
#' its confidence score.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param cell_id Cell barcode.
#'
#' @return Character vector describing the factors contributing to the
#' confidence score.
#'
#' @examples
#' \dontrun{
#' ## Assume 'object' has already been processed with cell_certify().
#' explain_confidence(
#'   object,
#'   cell_id = "Cell1"
#' )
#' }
#'
#' @export

explain_confidence <- function(
    object,
    cell_id
) {

  if (inherits(object, "Seurat")) {

    meta <- object@meta.data

  } else if (inherits(object, "SingleCellExperiment")) {

    meta <- as.data.frame(
      SummarizedExperiment::colData(object)
    )

  } else {

    stop(
      "'object' must be a Seurat or SingleCellExperiment object.",
      call. = FALSE
    )

  }

  if (!cell_id %in% rownames(meta)) {

    stop(
      sprintf(
        "Cell '%s' not found.",
        cell_id
      ),
      call. = FALSE
    )

  }

  meta <- meta[
    cell_id,
    ,
    drop = FALSE
  ]

  reasons <- character()

  ## Marker contribution
  if (meta$marker_score < 0.3) {

    reasons <- c(
      reasons,
      "Weak marker enrichment"
    )

  }

  ## Neighborhood agreement
  if (meta$neighbor_score < 0.5) {

    reasons <- c(
      reasons,
      "Neighborhood disagreement"
    )

  }

  ## Entropy
  if (meta$entropy_norm > 0.7) {

    reasons <- c(
      reasons,
      "High annotation uncertainty"
    )

  }

  ## Doublet score
  if ("doublet_score" %in% colnames(meta)) {

    if (meta$doublet_score > 0.5) {

      reasons <- c(
        reasons,
        "Possible doublet contamination"
      )

    }

  }

  ## High confidence
  if (length(reasons) == 0) {

    reasons <- c(
      reasons,
      "Strong annotation support"
    )

  }

  reasons

}
