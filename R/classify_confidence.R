#' Classify Confidence Levels
#'
#' Converts confidence scores into
#' interpretable confidence classes.
#'
#' @param scores Numeric confidence scores
#'
#' @return Character vector
#'
#' @export

classify_confidence <- function(
    scores
) {

  classes <- ifelse(
    scores >= 0.8,

    "High",

    ifelse(
      scores >= 0.5,

      "Moderate",

      "Low"
    )
  )

  return(classes)
}
