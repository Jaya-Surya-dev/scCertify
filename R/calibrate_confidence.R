#' Calibrate Confidence Scores
#'
#' Converts raw confidence scores
#' into percentile-calibrated scores.
#'
#' @param scores Numeric confidence scores
#'
#' @return Numeric calibrated scores
#'
#' @export

calibrate_confidence <- function(
    scores
) {

  ranks <- rank(
    scores,
    ties.method = "average"
  )

  calibrated <- ranks /
    max(ranks)

  return(calibrated)
}
