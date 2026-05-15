#' Neighbor Agreement Score
#'
#' Computes local label agreement
#'
#' @export

neighbor_score <- function(
    object,
    reduction = "pca",
    dims = NULL,
    k = 10,
    label_column = "predicted_label"
) {

  embeddings <- Embeddings(
    object,
    reduction
  )

  # Automatically determine dimensions
  if(is.null(dims)) {

    max_dims <- ncol(embeddings)

    dims <- 1:min(
      10,
      max_dims
    )
  }

  embeddings <- embeddings[
    ,
    dims,
    drop = FALSE
  ]

  nn <- FNN::get.knn(
    embeddings,
    k = min(
      k,
      nrow(embeddings) - 1
    )
  )

  labels <- object@meta.data[
    ,
    label_column
  ]

  scores <- numeric(
    length(labels)
  )

  for(i in seq_along(labels)) {

    neighbor_labels <- labels[
      nn$nn.index[i, ]
    ]

    scores[i] <- mean(
      neighbor_labels ==
        labels[i]
    )
  }

  return(scores)
}
