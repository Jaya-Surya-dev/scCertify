#' Marker Consistency Score
#'
#' Uses UCell enrichment scoring
#' for robust marker evaluation.
#'
#' @export

marker_score <- function(
    object,
    markers,
    label_column = "predicted_label"
) {

  labels <- object@meta.data[
    ,
    label_column
  ]

  scores <- numeric(
    length(labels)
  )

  for(label in unique(labels)) {

    matched_label <- match_labels(
      label,
      names(markers)
    )

    if(is.na(matched_label)) {

      warning(
        paste(
          "No ontology match for:",
          label
        )
      )

      next
    }

    genes <- markers[[matched_label]]

    genes <- genes[
      genes %in% rownames(object)
    ]

    if(length(genes) == 0) {

      warning(
        paste(
          "No marker genes found for:",
          matched_label
        )
      )

      next
    }

    signature_name <- paste0(
      matched_label,
      "_Signature"
    )

    feature_list <- list()

    feature_list[[signature_name]] <- genes

    object <- AddModuleScore_UCell(
      object,
      features = feature_list
    )

    score_column <- paste0(
      signature_name,
      "_UCell"
    )

    idx <- which(
      labels == label
    )

    scores[idx] <- object@meta.data[
      idx,
      score_column
    ]
  }

  return(scores)
}
