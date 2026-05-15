#' Marker Consistency Score
#'
#' Calculates marker confidence
#' using automatic ontology matching.
#'
#' @export

marker_score <- function(
    object,
    markers,
    label_column = "predicted_label"
) {

  scores <- numeric(
    ncol(object)
  )

  labels <- object@meta.data[
    ,
    label_column
  ]

  for(i in seq_along(labels)) {

    label <- as.character(
      labels[i]
    )

    matched_label <- match_labels(
      label,
      names(markers)
    )

    # No ontology match
    if(is.na(matched_label)) {

      warning(
        paste(
          "No ontology match for:",
          label
        )
      )

      scores[i] <- 0

      next
    }

    genes <- markers[[matched_label]]

    genes <- genes[
      genes %in% rownames(object)
    ]

    # No marker genes found
    if(length(genes) == 0) {

      warning(
        paste(
          "No marker genes detected for:",
          matched_label
        )
      )

      scores[i] <- 0

      next
    }

    expr <- FetchData(
      object,
      vars = genes,
      cells = colnames(object)[i]
    )

    expr <- as.numeric(expr)

    expr <- expr[
      !is.na(expr)
    ]

    if(length(expr) == 0) {

      scores[i] <- 0

      next
    }

    scores[i] <- mean(expr)
  }

  return(scores)
}
