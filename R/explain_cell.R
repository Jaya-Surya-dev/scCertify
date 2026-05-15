#' Explain Cell Annotation Confidence
#'
#' Provides interpretable explanation
#' for confidence scores.
#'
#' @param object Seurat object
#' @param cell_id cell barcode
#'
#' @export

explain_cell <- function(
    object,
    cell_id
) {

  meta <- object@meta.data[
    cell_id,
    ,
    drop = FALSE
  ]

  cat(
    "Cell:",
    cell_id,
    "\n\n"
  )

  cat(
    "Predicted Label:",
    meta$predicted_label,
    "\n"
  )

  cat(
    "Confidence Score:",
    round(
      meta$confidence_score,
      3
    ),
    "\n"
  )

  cat(
    "Marker Score:",
    round(
      meta$marker_score,
      3
    ),
    "\n"
  )

  cat(
    "Neighbor Agreement:",
    round(
      meta$neighbor_score,
      3
    ),
    "\n"
  )

  cat(
    "Entropy:",
    round(
      meta$entropy_norm,
      3
    ),
    "\n\n"
  )

  cat(
    "Interpretation:\n"
  )

  if(meta$confidence_score < 0.4) {

    cat(
      "- Low confidence annotation\n"
    )

  } else if(
    meta$confidence_score < 0.7
  ) {

    cat(
      "- Moderate confidence annotation\n"
    )

  } else {

    cat(
      "- High confidence annotation\n"
    )
  }

  if(meta$entropy_norm > 0.7) {

    cat(
      "- High uncertainty detected\n"
    )
  }

  if(meta$neighbor_score < 0.5) {

    cat(
      "- Local neighborhood disagreement\n"
    )
  }

  if(meta$marker_score < 0.2) {

    cat(
      "- Weak marker support\n"
    )
  }
}
