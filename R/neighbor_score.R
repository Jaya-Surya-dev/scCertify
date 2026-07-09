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
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(
#'     rpois(1000, lambda = 5),
#'     nrow = 50,
#'     ncol = 20
#'   )
#'
#'   rownames(counts) <- paste0("Gene", seq_len(50))
#'   colnames(counts) <- paste0("Cell", seq_len(20))
#'
#'   obj <- Seurat::CreateSeuratObject(
#'     counts = counts
#'   )
#'
#'   obj$predicted_label <- rep(
#'     c("T cell", "B cell"),
#'     each = 10
#'   )
#'
#'   obj <- Seurat::NormalizeData(
#'     obj,
#'     verbose = FALSE
#'   )
#'
#'   obj <- Seurat::FindVariableFeatures(
#'     obj,
#'     verbose = FALSE
#'   )
#'
#'   obj <- Seurat::ScaleData(
#'     obj,
#'     verbose = FALSE
#'   )
#'
#'   obj <- Seurat::RunPCA(
#'     obj,
#'     npcs = 10,
#'     verbose = FALSE
#'   )
#'
#'   scores <- neighbor_score(obj)
#'
#'   head(scores)
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
