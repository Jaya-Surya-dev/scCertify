#' Neighbor Agreement Score
#'
#' Computes local neighborhood agreement scores.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param reduction Dimensionality reduction to use.
#' @param dims Dimensions to include.
#' @param k Number of nearest neighbors.
#' @param label_column Metadata column containing predicted labels.
#'
#' @return Numeric vector of neighborhood agreement scores.
#'
#' @examples
#' NULL
#'
#' @export

neighbor_score <- function(
    object,
    reduction = "pca",
    dims = NULL,
    k = 10,
    label_column = "predicted_label"
) {

  if (inherits(object, "Seurat")) {

    return(
      neighbor_score_core(
        object = object,
        reduction = reduction,
        dims = dims,
        k = k,
        label_column = label_column,
        object_type = "Seurat"
      )
    )

  }

  if (inherits(object, "SingleCellExperiment")) {

    return(
      neighbor_score_core(
        object = object,
        reduction = reduction,
        dims = dims,
        k = k,
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
