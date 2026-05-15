#' Cell Annotation Confidence
#'
#' Main confidence framework
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
    0.4 * scale(
      object$marker_score
    ) +
    0.4 * scale(
      object$neighbor_score
    ) +
    0.2 * scale(
      entropy_component
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
    confidence

  return(object)
}
