#' Discovery-Aware Annotation Status
#'
#' Classifies cells according to positive marker evidence
#' and contradictory negative marker evidence.
#'
#' @param marker_score Numeric vector of positive marker scores.
#' @param negative_marker_score Numeric vector of negative marker scores.
#' @param positive_threshold Threshold for strong positive evidence.
#' Must be a single numeric value between 0 and 1.
#' @param negative_threshold Threshold for strong contradictory evidence.
#' Must be a single numeric value between 0 and 1.
#'
#' @return Character vector containing one discovery-aware status per cell.
#'
#' @details
#' The function provides a simple exploratory interpretation of
#' annotation evidence. It does not modify the main confidence score
#' and does not replace biological validation.
#'
#' Cells with strong positive evidence and weak contradictory evidence
#' are classified as "Known". Cells with weak positive and weak
#' contradictory evidence are classified as "Possible novel state".
#' Cells with both strong positive and contradictory evidence are
#' classified as "Possible transitional state". Remaining cells are
#' classified as "Insufficient evidence".
#'
#' @export

discovery_status <- function(
    marker_score,
    negative_marker_score,
    positive_threshold = 0.50,
    negative_threshold = 0.50
) {

  if (length(marker_score) !=
      length(negative_marker_score)) {
    stop(
      "'marker_score' and 'negative_marker_score' ",
      "must have the same length.",
      call. = FALSE
    )
  }

  if (anyNA(marker_score) ||
      anyNA(negative_marker_score)) {
    stop(
      "Scores must not contain NA values.",
      call. = FALSE
    )
  }

  if (!is.numeric(positive_threshold) ||
      length(positive_threshold) != 1 ||
      is.na(positive_threshold) ||
      positive_threshold < 0 ||
      positive_threshold > 1) {
    stop(
      "'positive_threshold' must be a single numeric value ",
      "between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(negative_threshold) ||
      length(negative_threshold) != 1 ||
      is.na(negative_threshold) ||
      negative_threshold < 0 ||
      negative_threshold > 1) {
    stop(
      "'negative_threshold' must be a single numeric value ",
      "between 0 and 1.",
      call. = FALSE
    )
  }

  status <- rep(
    "Insufficient evidence",
    length(marker_score)
  )

  strong_positive <-
    marker_score >= positive_threshold

  strong_negative <-
    negative_marker_score >= negative_threshold

  status[
    strong_positive &
      !strong_negative
  ] <- "Known"

  status[
    !strong_positive &
      !strong_negative
  ] <- "Possible novel state"

  status[
    strong_positive &
      strong_negative
  ] <- "Possible transitional state"

  status
}
