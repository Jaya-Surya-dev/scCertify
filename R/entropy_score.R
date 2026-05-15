#' Entropy Uncertainty Score
#'
#' Calculates uncertainty from
#' SingleR score matrix.
#'
#' @param score_matrix SingleR score matrix
#'
#' @export

entropy_score <- function(score_matrix) {

  entropy_values <- apply(
    score_matrix,
    1,
    function(x) {

      probs <- x / sum(x)

      probs <- probs[probs > 0]

      ent <- entropy::entropy(probs)

      return(ent)
    }
  )

  return(entropy_values)
}
