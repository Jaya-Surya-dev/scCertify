neighbor_score_core <- function(
    object,
    reduction,
    dims,
    k,
    label_column,
    object_type
) {

  if (object_type == "Seurat") {

    embeddings <- Seurat::Embeddings(
      object,
      reduction = reduction
    )

    metadata <- object@meta.data

  } else {

    rd_names <- SingleCellExperiment::reducedDimNames(
      object
    )

    match_idx <- match(
      tolower(reduction),
      tolower(rd_names)
    )

    if (is.na(match_idx)) {

      stop(
        sprintf(
          "Reduction '%s' not found.",
          reduction
        ),
        call. = FALSE
      )

    }

    embeddings <- SingleCellExperiment::reducedDim(
      object,
      rd_names[match_idx]
    )

    metadata <- as.data.frame(
      SummarizedExperiment::colData(
        object
      )
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

  if (is.null(dims)) {

    dims <- seq_len(
      min(
        10,
        ncol(embeddings)
      )
    )

  }

  embeddings <- embeddings[
    ,
    dims,
    drop = FALSE
  ]

  k <- min(
    k,
    nrow(embeddings) - 1
  )

  nn <- FNN::get.knn(
    embeddings,
    k = k
  )

  labels <- metadata[[label_column]]

  scores <- numeric(
    length(labels)
  )

  for (i in seq_along(labels)) {

    neighbor_labels <- labels[
      nn$nn.index[i, ]
    ]

    scores[i] <- mean(
      neighbor_labels == labels[i]
    )

  }

  scores

}
