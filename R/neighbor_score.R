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
#' \dontrun{
#' library(Seurat)
#'
#' counts <- matrix(
#'   rpois(200, lambda = 5),
#'   nrow = 20
#' )
#'
#' rownames(counts) <- paste0("Gene", 1:20)
#' colnames(counts) <- paste0("Cell", 1:10)
#'
#' obj <- CreateSeuratObject(counts)
#' obj$predicted_label <- rep(c("T cell", "B cell"), each = 5)
#'
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj)
#'
#' neighbor_score(obj)
#' }
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
