
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
#' set.seed(123)
#'
#' counts <- matrix(
#'   rpois(200, lambda = 5),
#'   nrow = 20,
#'   ncol = 10,
#'   dimnames = list(
#'     paste0("Gene", seq_len(20)),
#'     paste0("Cell", seq_len(10))
#'   )
#' )
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     logcounts = log1p(counts)
#'   )
#' )
#'
#' SummarizedExperiment::colData(sce)$predicted_label <-
#'   rep(c("T cell", "B cell"), each = 5)
#'
#' markers <- list(
#'   "T cell" = c("Gene1", "Gene2"),
#'   "B cell" = c("Gene3", "Gene4")
#' )
#'
#' scores <- marker_score(
#'   sce,
#'   markers,
#'   label_column = "predicted_label"
#' )
#'
#' head(scores)
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
