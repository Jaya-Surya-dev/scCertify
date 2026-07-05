#' Marker Consistency Score
#'
#' Uses UCell enrichment scoring
#' for marker evaluation.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param markers Named list of marker genes.
#' @param label_column Metadata column containing predicted labels.
#'
#' @return Numeric vector of marker consistency scores.
#'
#' @examples
#' NULL
#'
#' @export

marker_score <- function(
    object,
    markers,
    label_column = "predicted_label"
) {

  if (inherits(object, "Seurat")) {

    return(
      marker_score_core(
        object = object,
        markers = markers,
        label_column = label_column,
        object_type = "Seurat"
      )
    )

  }

  if (inherits(object, "SingleCellExperiment")) {

    return(
      marker_score_core(
        object = object,
        markers = markers,
        label_column = label_column,
        object_type = "SingleCellExperiment"
      )
    )

  }

  stop(
    "'object' must be a Seurat or SingleCellExperiment object.",
    call. = FALSE
  )

}
