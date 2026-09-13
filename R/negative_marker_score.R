#' Negative Marker Consistency Score
#'
#' Uses UCell enrichment scoring to evaluate expression of genes
#' that are inconsistent with a predicted cell identity.
#'
#' A higher score indicates stronger expression of negative markers
#' and therefore greater evidence against the predicted annotation.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param negative_markers Named list of negative marker genes.
#' @param label_column Metadata column containing predicted labels.
#'
#' @return Numeric vector of negative marker scores.
#'
#' @details
#' Negative markers provide contradictory molecular evidence for a
#' predicted cell identity. For example, genes characteristic of B
#' cells can be supplied as negative markers for cells predicted to
#' be T cells.
#'
#' The function returns one UCell score per cell. Scores are not
#' incorporated into the main confidence score by default.
#'
#' If a predicted label has no corresponding negative marker set,
#' or none of the supplied negative marker genes are present in the
#' object, a score of zero is returned for those cells. A zero score
#' therefore indicates that no measurable negative-marker signal
#' was obtained and should not be interpreted as biological evidence
#' supporting the predicted identity.
#'
#' @export
negative_marker_score <- function(
    object,
    negative_markers,
    label_column = "predicted_label"
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

  labels <- metadata[[label_column]]
  scores <- numeric(length(labels))
  gene_names <- rownames(object)

  for (label in unique(labels)) {
    matched_label <- match_labels(
      label,
      names(negative_markers)
    )

    if (is.na(matched_label)) {
      next
    }

    genes <- negative_markers[[matched_label]]
    genes <- genes[genes %in% gene_names]

    if (length(genes) == 0) {
      next
    }

    signature_name <- paste0(
      matched_label,
      "_NegativeSignature"
    )

    feature_list <- list()
    feature_list[[signature_name]] <- genes

    if (inherits(object, "Seurat")) {
      object <- UCell::AddModuleScore_UCell(
        object,
        features = feature_list,
        ncores = 1
      )

      metadata <- object@meta.data

      score_column <- paste0(
        signature_name,
        "_UCell"
      )

      idx <- which(labels == label)

      scores[idx] <- metadata[
        idx,
        score_column
      ]
    } else {
      object <- UCell::ScoreSignatures_UCell(
        object,
        features = feature_list,
        assay = "logcounts",
        ncores = 1
      )

      ucell <- SingleCellExperiment::altExp(
        object,
        "UCell"
      )

      u_scores <- SummarizedExperiment::assay(
        ucell,
        "UCell"
      )

      score_row <- paste0(
        signature_name,
        "_UCell"
      )

      idx <- which(labels == label)

      if (score_row %in% rownames(u_scores)) {
        scores[idx] <- as.numeric(
          u_scores[
            score_row,
            idx
          ]
        )
      }
    }
  }

  scores
}
