
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
#' @param run_discovery Logical. Whether to run the optional
#'   discovery-aware annotation analysis. Defaults to FALSE.
#' @param negative_markers Named list of negative marker genes,
#'   required when run_discovery = TRUE.
#'
#' @return The input object with confidence scores and classifications
#' added to the metadata. When run_discovery = TRUE, discovery-aware
#' scores and statuses are also added to the metadata.
#'
#' @examples
#' set.seed(123)
#'
#' counts <- matrix(
#'   rpois(1000, lambda = 5),
#'   nrow = 50,
#'   ncol = 20,
#'   dimnames = list(
#'     paste0("Gene", seq_len(50)),
#'     paste0("Cell", seq_len(20))
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
#'   rep(c("T cell", "B cell"), each = 10)
#'
#' markers <- list(
#'   "T cell" = c("Gene1", "Gene2"),
#'   "B cell" = c("Gene3", "Gene4")
#' )
#'
#' negative_markers <- list(
#'   "T cell" = c("Gene3", "Gene4"),
#'   "B cell" = c("Gene1", "Gene2")
#' )
#'
#' SingleCellExperiment::reducedDim(sce, "PCA") <- prcomp(
#'   t(counts),
#'   rank. = 5
#' )$x
#'
#' # Original workflow: discovery analysis is disabled by default.
#' sce <- cell_certify(
#'   sce,
#'   markers
#' )
#'
#' head(
#'   SummarizedExperiment::colData(sce)$confidence_score
#' )
#'
#' # Optional discovery-aware workflow.
#' sce_discovery <- cell_certify(
#'   sce,
#'   markers,
#'   run_discovery = TRUE,
#'   negative_markers = negative_markers
#' )
#'
#' head(
#'   SummarizedExperiment::colData(sce_discovery)$discovery_status
#' )
#'
#' @export
cell_certify <- function(
    object,
    markers,
    label_column = "predicted_label",
    moderate_threshold = 0.5,
    high_threshold = 0.8,
    run_discovery = FALSE,
    negative_markers = NULL
) {
  # Validate the input object and extract metadata.
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

  # Validate the discovery option.
  if (
    !is.logical(run_discovery) ||
    length(run_discovery) != 1L ||
    is.na(run_discovery)
  ) {
    stop(
      "'run_discovery' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (run_discovery && is.null(negative_markers)) {
    stop(
      "'negative_markers' must be supplied when ",
      "'run_discovery = TRUE'.",
      call. = FALSE
    )
  }

  # Validate the label column.
  if (!label_column %in% colnames(metadata)) {
    stop(
      sprintf(
        "'%s' not found in object metadata.",
        label_column
      ),
      call. = FALSE
    )
  }

  # Supply defaults for missing uncertainty metadata.
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

  # Synchronize metadata with the input object.
  if (inherits(object, "Seurat")) {
    object$entropy_norm <- metadata$entropy_norm
    object$doublet_score <- metadata$doublet_score
  } else {
    SummarizedExperiment::colData(object)$entropy_norm <-
      metadata$entropy_norm

    SummarizedExperiment::colData(object)$doublet_score <-
      metadata$doublet_score
  }

  # Calculate the original confidence components.
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

  confidence_score <- calibrate_confidence(
    confidence
  )

  confidence_class <- classify_confidence(
    confidence_score,
    moderate_threshold = moderate_threshold,
    high_threshold = high_threshold
  )

  # Preserve the original confidence outputs.
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

  # Run discovery-aware analysis only when explicitly requested.
  # These outputs remain separate from the confidence score.
  if (run_discovery) {
    discovery_results <- discovery_aware(
      object = object,
      markers = markers,
      negative_markers = negative_markers,
      label_column = label_column
    )

    if (inherits(object, "Seurat")) {
      object$discovery_positive_marker_score <-
        discovery_results$positive_marker_score

      object$discovery_negative_marker_score <-
        discovery_results$negative_marker_score

      object$discovery_status <-
        discovery_results$discovery_status
    } else {
      SummarizedExperiment::colData(object)$
        discovery_positive_marker_score <-
        discovery_results$positive_marker_score

      SummarizedExperiment::colData(object)$
        discovery_negative_marker_score <-
        discovery_results$negative_marker_score

      SummarizedExperiment::colData(object)$
        discovery_status <-
        discovery_results$discovery_status
    }
  }

  object
}
