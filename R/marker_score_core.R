marker_score_core <- function(
    object,
    markers,
    label_column,
    object_type
) {

  if (object_type == "Seurat") {

    metadata <- object@meta.data

  } else {

    metadata <- as.data.frame(
      SummarizedExperiment::colData(object)
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
      names(markers)
    )

    if (is.na(matched_label)) {

      warning(
        sprintf(
          "No ontology match for: %s",
          label
        ),
        call. = FALSE
      )

      next

    }

    genes <- markers[[matched_label]]

    genes <- genes[
      genes %in% gene_names
    ]

    if (length(genes) == 0) {

      warning(
        sprintf(
          "No marker genes found for: %s",
          matched_label
        ),
        call. = FALSE
      )

      next

    }

    signature_name <- paste0(
      matched_label,
      "_Signature"
    )

    feature_list <- list()
    feature_list[[signature_name]] <- genes

    if (object_type == "Seurat") {

      object <- UCell::AddModuleScore_UCell(
        object,
        features = feature_list,
        ncores = 1
      )

      metadata <- object@meta.data

    } else {

      object <- UCell::ScoreSignatures_UCell(
        object,
        features = feature_list,
        assay = "logcounts"
      )

      metadata <- as.data.frame(
        SummarizedExperiment::colData(object)
      )

    }

    score_column <- paste0(
      signature_name,
      "_UCell"
    )

    idx <- which(
      labels == label
    )

    scores[idx] <- metadata[
      idx,
      score_column
    ]

  }

  scores

}
