#' Cell Annotation Confidence
#'
#' Main confidence scoring framework.
#'
#' @param object Seurat object
#' @param markers Named marker gene list
#' @param label_column Metadata column containing labels
#'
#' @return Seurat object with confidence scores
#'
#' @export

cell_certify <- function(
    object,
    markers,
    label_column = "predicted_label"
) {

  object$marker_score <-
    marker_score(
      object,
      markers,
      label_column
    )

  object$neighbor_score <-
    neighbor_score(
      object,
      label_column =
        label_column
    )

  entropy_component <-
    1 - object$entropy_norm

  confidence <-
    0.35 * scale(
      object$marker_score
    ) +
    0.35 * scale(
      object$neighbor_score
    ) +
    0.2 * scale(
      entropy_component
    ) -
    0.1 * scale(
      object$doublet_score
    )

  confidence <- as.numeric(
    scale(confidence)
  )

  confidence <- (
    confidence -
      min(confidence)
  ) / (
    max(confidence) -
      min(confidence)
  )

  object$confidence_score <-
    calibrate_confidence(
      confidence
    )

  object$confidence_class <-
    classify_confidence(
      object$confidence_score
    )

  return(object)
}
