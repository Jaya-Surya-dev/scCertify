#' Explain Confidence Attribution
#'
#' Explains why a cell received
#' its confidence score.
#'
#' @param object Seurat object
#' @param cell_id Cell barcode
#'
#' @return Character vector
#'
#' @export

explain_confidence <- function(
    object,
    cell_id
) {

  meta <- object@meta.data[
    cell_id,
    ,
    drop = FALSE
  ]

  reasons <- c()

  # Marker contribution
  if(meta$marker_score < 0.3) {

    reasons <- c(
      reasons,
      "Weak marker enrichment"
    )
  }

  # Neighborhood agreement
  if(meta$neighbor_score < 0.5) {

    reasons <- c(
      reasons,
      "Neighborhood disagreement"
    )
  }

  # Entropy
  if(meta$entropy_norm > 0.7) {

    reasons <- c(
      reasons,
      "High annotation uncertainty"
    )
  }

  # Doublet score
  if("doublet_score" %in%
     colnames(meta)) {

    if(meta$doublet_score > 0.5) {

      reasons <- c(
        reasons,
        "Possible doublet contamination"
      )
    }
  }

  # High confidence
  if(length(reasons) == 0) {

    reasons <- c(
      reasons,
      "Strong annotation support"
    )
  }

  return(reasons)
}
