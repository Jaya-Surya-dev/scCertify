#' Classify Confidence Levels
#'
#' Converts confidence scores into
#' confidence classes.
#'
#' @param scores Numeric vector of confidence scores.
#' @param moderate_threshold Lower threshold for the "Moderate"
#'   confidence class.
#' @param high_threshold Lower threshold for the "High"
#'   confidence class.
#'
#' @return Character vector of confidence classes.
#'
#' @examples
#' scores <- c(
#'   0.2,
#'   0.7,
#'   0.9
#' )
#'
#' classify_confidence(
#'   scores
#' )
#'
#' classify_confidence(
#'   scores,
#'   moderate_threshold = 0.4,
#'   high_threshold = 0.75
#' )
#'
#' @export

classify_confidence <- function(
    scores,
    moderate_threshold = 0.5,
    high_threshold = 0.8
) {

  if (!is.numeric(scores)) {

    stop(
      "'scores' must be numeric.",
      call. = FALSE
    )

  }

  if (moderate_threshold >= high_threshold) {

    stop(
      "'moderate_threshold' must be smaller than 'high_threshold'.",
      call. = FALSE
    )

  }

  classes <- ifelse(
    scores >= high_threshold,
    "High",
    ifelse(
      scores >= moderate_threshold,
      "Moderate",
      "Low"
    )
  )

  classes

}
