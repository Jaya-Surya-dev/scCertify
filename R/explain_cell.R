#' Explain Cell Annotation Confidence
#'
#' Provides an interpretable explanation for the confidence score
#' assigned to an individual cell.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param cell_id Cell barcode or column name.
#'
#' @return A character vector containing the interpretation.
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
#'   obj$predicted_label <- rep(
#'     c("T cell", "B cell"),
#'     each = 5
#'   )
#'
#'   obj$marker_score <- runif(ncol(obj))
#'   obj$neighbor_score <- runif(ncol(obj))
#'   obj$entropy_norm <- runif(ncol(obj))
#'   obj$confidence_score <- runif(ncol(obj))
#'
#'   explain_cell(
#'     obj,
#'     cell_id = "Cell1"
#'   )
#' }
#'
#' @export

explain_cell <- function(
    object,
    cell_id
) {

  if (inherits(object, "Seurat")) {

    meta <- object@meta.data

  } else if (inherits(object, "SingleCellExperiment")) {

    meta <- as.data.frame(
      SummarizedExperiment::colData(object)
    )

  } else {

    stop(
      "'object' must be a Seurat or SingleCellExperiment object.",
      call. = FALSE
    )

  }

  if (!cell_id %in% rownames(meta)) {

    stop(
      sprintf(
        "Cell '%s' not found.",
        cell_id
      ),
      call. = FALSE
    )

  }

  meta <- meta[
    cell_id,
    ,
    drop = FALSE
  ]

  message("Cell: ", cell_id, "\n")

  message(
    "Predicted Label: ",
    meta$predicted_label
  )

  message(
    "Confidence Score: ",
    round(
      meta$confidence_score,
      3
    )
  )

  message(
    "Marker Score: ",
    round(
      meta$marker_score,
      3
    )
  )

  message(
    "Neighbor Agreement: ",
    round(
      meta$neighbor_score,
      3
    )
  )

  message(
    "Entropy: ",
    round(
      meta$entropy_norm,
      3
    ),
    "\n"
  )

  interpretation <- character()

  if (meta$confidence_score < 0.4) {

    interpretation <- c(
      interpretation,
      "Low confidence annotation"
    )

  } else if (meta$confidence_score < 0.7) {

    interpretation <- c(
      interpretation,
      "Moderate confidence annotation"
    )

  } else {

    interpretation <- c(
      interpretation,
      "High confidence annotation"
    )

  }

  if (meta$entropy_norm > 0.7) {

    interpretation <- c(
      interpretation,
      "High uncertainty detected"
    )

  }

  if (meta$neighbor_score < 0.5) {

    interpretation <- c(
      interpretation,
      "Local neighborhood disagreement"
    )

  }

  if (meta$marker_score < 0.2) {

    interpretation <- c(
      interpretation,
      "Weak marker support"
    )

  }

  message("Interpretation:")

  for (msg in interpretation) {

    message("- ", msg)

  }

  interpretation

}
