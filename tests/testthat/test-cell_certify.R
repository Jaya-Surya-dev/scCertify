test_that(
  "cell_certify works for Seurat",
  {

    obj <- toy_seurat

    obj$entropy_norm <- runif(
      ncol(obj)
    )

    obj$doublet_score <- runif(
      ncol(obj)
    )

    obj <- suppressWarnings(

      cell_certify(
        obj,
        toy_markers
      )

    )

    expect_true(
      "confidence_score" %in%
        colnames(obj@meta.data)
    )

    expect_true(
      "confidence_class" %in%
        colnames(obj@meta.data)
    )

    expect_length(
      obj$confidence_score,
      ncol(obj)
    )

  }
)

test_that(
  "cell_certify works for SingleCellExperiment",
  {

    obj <- toy_sce

    SummarizedExperiment::colData(
      obj
    )$entropy_norm <-
      runif(ncol(obj))

    SummarizedExperiment::colData(
      obj
    )$doublet_score <-
      runif(ncol(obj))

    obj <- suppressWarnings(

      cell_certify(
        obj,
        toy_markers
      )

    )

    expect_true(

      "confidence_score" %in%

        colnames(

          SummarizedExperiment::colData(obj)

        )

    )

    expect_true(

      "confidence_class" %in%

        colnames(

          SummarizedExperiment::colData(obj)

        )

    )

  }
)

test_that(
  "cell_certify errors for missing labels",
  {

    expect_error(

      cell_certify(
        toy_seurat,
        toy_markers,
        label_column = "missing"
      )

    )

  }
)
