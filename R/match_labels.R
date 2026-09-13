#' Ontology Label Matching
#'
#' Matches predicted labels
#' to marker database labels.
#'
#' @param label Character label
#' @param marker_names Character vector
#'
#' @return Character vector
#'
#' @examples
#' match_labels(
#'   "T_cells",
#'   c(
#'     "T_cell",
#'     "B_cell"
#'   )
#' )
#'
#' @export

match_labels <- function(
    label,
    marker_names
) {
  label_clean <- tolower(label)

  label_clean <- gsub(
    "_",
    "",
    label_clean
  )

  label_clean <- gsub(
    "s$",
    "",
    label_clean
  )

  marker_clean <- tolower(
    marker_names
  )

  marker_clean <- gsub(
    "_",
    "",
    marker_clean
  )

  marker_clean <- gsub(
    "s$",
    "",
    marker_clean
  )

  exact_match <- marker_names[
    marker_clean ==
      label_clean
  ]

  if (length(exact_match) > 0) {
    return(exact_match[1])
  }

  partial_match <- marker_names[
    vapply(
      marker_clean,
      function(x) {
        grepl(
          label_clean,
          x,
          fixed = TRUE
        ) ||
          grepl(
            x,
            label_clean,
            fixed = TRUE
          )
      },
      logical(1)
    )
  ]

  if (length(partial_match) > 0) {
    return(partial_match[1])
  }

  return(NA)
}
