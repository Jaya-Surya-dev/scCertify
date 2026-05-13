#' Neighbor Agreement Score
#'
#' Computes local label agreement.
#'
#' @export

neighbor_score <- function(
    object,
    reduction = "pca",
    dims = 1:20,
    k = 20,
    label_column = "predicted_label"
) {

  embeddings <- Embeddings(
    object,
    reduction
  )[, dims]

  nn <- FNN::get.knn(
    embeddings,
    k = k
  )

  labels <- object@meta.data[
    ,
    label_column
  ]

  scores <- numeric(length(labels))

  for(i in 1:length(labels)) {

    neighbor_labels <- labels[
      nn$nn.index[i, ]
    ]

    scores[i] <- mean(
      neighbor_labels == labels[i]
    )
  }

  return(scores)
}
