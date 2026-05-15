#' Match Cell Labels to Marker Ontology
#'
#' Automatically matches predicted
#' labels to available marker sets.
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

  if(length(exact_match) > 0) {

    return(exact_match[1])
  }

  partial_match <- marker_names[
    grepl(
      label_clean,
      marker_clean
    ) |

      grepl(
        marker_clean,
        label_clean
      )
  ]

  if(length(partial_match) > 0) {

    return(partial_match[1])
  }

  return(NA)
}
