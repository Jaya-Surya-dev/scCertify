#' Cell Annotation Confidence
#'
#' Main confidence scoring pipeline.
#'
#' @export

cell_certify <- function(
    object,
    markers,
    label_column = "predicted_label"
) {

  object$marker_score <- marker_score(
    object,
    markers,
    label_column
  )

  object$neighbor_score <- neighbor_score(
    object,
    label_column = label_column
  )

  confidence <-
    0.5 * scale(object$marker_score) +
    0.5 * scale(object$neighbor_score)

  confidence <- as.numeric(
    scale(confidence)
  )

  confidence <- (
    confidence - min(confidence)
  ) / (
    max(confidence) -
      min(confidence)
  )

  object$confidence_score <-
    confidence

  return(object)
}
