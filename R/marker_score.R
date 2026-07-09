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
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(
#'     rpois(200, lambda = 5),
#'     nrow = 20
#'   )
#'
#'   rownames(counts) <- paste0("Gene", seq_len(20))
#'   colnames(counts) <- paste0("Cell", seq_len(10))
#'
#'   obj <- Seurat::CreateSeuratObject(
#'     counts = counts
#'   )
#'
#'   obj <- Seurat::NormalizeData(
#'     obj,
#'     verbose = FALSE
#'   )
#'
#'   obj$predicted_label <- rep(
#'     c("T cell", "B cell"),
#'     each = 5
#'   )
#'
#'   markers <- list(
#'     "T cell" = c("Gene1", "Gene2"),
#'     "B cell" = c("Gene3", "Gene4")
#'   )
#'
#'   scores <- marker_score(
#'     obj,
#'     markers
#'   )
#'
#'   head(scores)
#' }
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
