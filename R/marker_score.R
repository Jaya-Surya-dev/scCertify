#' Marker Consistency Score
#'
#' Calculates marker-based confidence score.
#'
#' @param object Seurat object
#' @param markers named marker list
#' @param label_column metadata column
#'
#' @export

marker_score <- function(
    object,
    markers,
    label_column = "predicted_label"
) {

  scores <- numeric(ncol(object))

  for(i in 1:ncol(object)) {

    label <- object@meta.data[i, label_column]

    if(!label %in% names(markers)) {
      scores[i] <- NA
      next
    }

    genes <- markers[[label]]

    genes <- genes[
      genes %in% rownames(object)
    ]

    if(length(genes) == 0) {
      scores[i] <- 0
      next
    }

    expr <- FetchData(
      object,
      vars = genes,
      cells = colnames(object)[i]
    )

    scores[i] <- mean(
      as.numeric(expr)
    )
  }

  return(scores)
}
