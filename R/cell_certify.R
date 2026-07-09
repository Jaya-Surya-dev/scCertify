#' Cell Annotation Confidence
#'
#' Computes confidence scores for cell-type annotations using marker
#' enrichment, neighborhood agreement, entropy-based uncertainty,
#' and doublet information.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param markers Named list of marker genes.
#' @param label_column Metadata column containing predicted labels.
#' @param moderate_threshold Lower threshold for the "Moderate"
#'   confidence class.
#' @param high_threshold Lower threshold for the "High"
#'   confidence class.
#'
#' @return The input object with confidence scores and classifications
#' added to the metadata.
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
#'   obj$entropy_norm <- runif(ncol(obj))
#'   obj$doublet_score <- runif(ncol(obj))
#'
#'   markers <- list(
#'     "T cell" = c("Gene1", "Gene2"),
#'     "B cell" = c("Gene3", "Gene4")
#'   )
#'
#'   obj <- cell_certify(
#'     obj,
#'     markers
#'   )
#'
#'   head(obj$confidence_score)
#' }
#'
#' @export

cell_certify <- function(
    object,
    markers,
    label_column = "predicted_label",
    moderate_threshold = 0.5,
    high_threshold = 0.8
) {

  if (inherits(object, "Seurat")) {

    metadata <- object@meta.data

  } else if (inherits(object, "SingleCellExperiment")) {

    metadata <- as.data.frame(
      SummarizedExperiment::colData(object)
    )

  } else {

    stop(
      "'object' must be a Seurat or SingleCellExperiment object.",
      call. = FALSE
    )

  }

  if (!label_column %in% colnames(metadata)) {

    stop(
      sprintf(
        "'%s' not found in object metadata.",
        label_column
      ),
      call. = FALSE
    )

  }

  if (!"entropy_norm" %in% colnames(metadata)) {

    warning(
      "'entropy_norm' not found. Using zeros.",
      call. = FALSE
    )

    metadata$entropy_norm <- 0

  }

  if (!"doublet_score" %in% colnames(metadata)) {

    warning(
      "'doublet_score' not found. Using zeros.",
      call. = FALSE
    )

    metadata$doublet_score <- 0

  }

  if (inherits(object, "Seurat")) {

    object$entropy_norm <- metadata$entropy_norm
    object$doublet_score <- metadata$doublet_score

  } else {

    SummarizedExperiment::colData(object)$entropy_norm <-
      metadata$entropy_norm

    SummarizedExperiment::colData(object)$doublet_score <-
      metadata$doublet_score

  }

  marker_scores <- marker_score(
    object,
    markers,
    label_column
  )

  neighbor_scores <- neighbor_score(
    object,
    label_column = label_column
  )

  entropy_component <- 1 - metadata$entropy_norm

  confidence <-
    0.35 * scale(marker_scores) +
    0.35 * scale(neighbor_scores) +
    0.20 * scale(entropy_component) -
    0.10 * scale(metadata$doublet_score)

  confidence <- as.numeric(
    scale(confidence)
  )

  confidence <- (
    confidence -
      min(confidence, na.rm = TRUE)
  ) / (
    max(confidence, na.rm = TRUE) -
      min(confidence, na.rm = TRUE)
  )

  confidence_score <-
    calibrate_confidence(
      confidence
    )

  confidence_class <-
    classify_confidence(
      confidence_score,
      moderate_threshold = moderate_threshold,
      high_threshold = high_threshold
    )

  if (inherits(object, "Seurat")) {

    object$marker_score <- marker_scores
    object$neighbor_score <- neighbor_scores
    object$confidence_score <- confidence_score
    object$confidence_class <- confidence_class

  } else {

    SummarizedExperiment::colData(object)$marker_score <-
      marker_scores

    SummarizedExperiment::colData(object)$neighbor_score <-
      neighbor_scores

    SummarizedExperiment::colData(object)$confidence_score <-
      confidence_score

    SummarizedExperiment::colData(object)$confidence_class <-
      confidence_class

  }

  object

}
