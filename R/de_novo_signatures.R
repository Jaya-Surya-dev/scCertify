#' De Novo Cluster Gene Signatures
#'
#' Identifies genes that are relatively enriched within each
#' user-specified cluster.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param cluster_column Metadata column containing cluster identities.
#' @param top_n Number of top genes to return for each cluster.
#'
#' @return Named list containing the top genes for each cluster.
#'
#' @details
#' This function provides a simple discovery-aware approach by
#' identifying genes enriched within existing clusters. It does not
#' perform clustering and does not assign biological identities.
#' The resulting signatures can be used to investigate unusual,
#' heterogeneous, or potentially novel cell states.
#'
#' @export

de_novo_signatures <- function(
    object,
    cluster_column,
    top_n = 10
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

  if (!cluster_column %in% colnames(metadata)) {
    stop(
      sprintf(
        "'%s' not found in object metadata.",
        cluster_column
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(top_n) ||
      length(top_n) != 1 ||
      is.na(top_n) ||
      top_n < 1) {
    stop(
      "'top_n' must be a single positive number.",
      call. = FALSE
    )
  }

  top_n <- as.integer(top_n)

  clusters <- metadata[[cluster_column]]

  if (anyNA(clusters)) {
    stop(
      "Cluster identities must not contain NA values.",
      call. = FALSE
    )
  }

  expression_data <- if (
    inherits(object, "Seurat")
  ) {

    assay_name <- Seurat::DefaultAssay(object)

    available_layers <- SeuratObject::Layers(
      object[[assay_name]]
    )

    if ("data" %in% available_layers) {

      Seurat::GetAssayData(
        object,
        assay = assay_name,
        layer = "data"
      )

    } else if ("counts" %in% available_layers) {

      Seurat::GetAssayData(
        object,
        assay = assay_name,
        layer = "counts"
      )

    } else {

      stop(
        "No 'data' or 'counts' layer found in the Seurat object.",
        call. = FALSE
      )
    }

  } else {

    SummarizedExperiment::assay(
      object,
      "logcounts"
    )
  }

  gene_names <- rownames(expression_data)

  signatures <- lapply(
    unique(clusters),
    function(cluster) {

      cells <- which(
        clusters == cluster
      )

      cluster_expression <- Matrix::rowMeans(
        expression_data[
          ,
          cells,
          drop = FALSE
        ]
      )

      other_cells <- which(
        clusters != cluster
      )

      if (length(other_cells) > 0) {

        background_expression <- Matrix::rowMeans(
          expression_data[
            ,
            other_cells,
            drop = FALSE
          ]
        )

        enrichment <- (
          cluster_expression -
            background_expression
        )

      } else {

        enrichment <- cluster_expression
      }

      ranked_genes <- order(
        enrichment,
        decreasing = TRUE
      )

      selected <- ranked_genes[
        seq_len(
          min(
            top_n,
            length(ranked_genes)
          )
        )
      ]

      gene_names[selected]
    }
  )

  names(signatures) <- unique(clusters)

  signatures
}
