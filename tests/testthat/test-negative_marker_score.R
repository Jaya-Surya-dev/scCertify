test_that("negative_marker_score returns numeric vector for Seurat", {

  scores <- negative_marker_score(
    toy_seurat,
    toy_negative_markers
  )

  expect_type(scores, "double")

  expect_length(
    scores,
    ncol(toy_seurat)
  )

  expect_true(
    all(scores >= 0)
  )
})


test_that(
  "negative_marker_score returns numeric vector for SingleCellExperiment",
  {

    scores <- negative_marker_score(
      toy_sce,
      toy_negative_markers
    )

    expect_type(scores, "double")

    expect_length(
      scores,
      ncol(toy_sce)
    )
  }
)


test_that(
  "negative_marker_score errors when label column is missing",
  {

    expect_error(
      negative_marker_score(
        toy_seurat,
        toy_negative_markers,
        label_column = "unknown"
      )
    )
  }
)


test_that(
  "negative_marker_score handles unknown labels",
  {

    tmp <- toy_seurat
    tmp$predicted_label <- "UnknownCell"

    scores <- negative_marker_score(
      tmp,
      toy_negative_markers
    )

    expect_type(scores, "double")

    expect_length(
      scores,
      ncol(tmp)
    )

    expect_true(
      all(scores == 0)
    )
  }
)


test_that(
  "negative_marker_score returns zero for unavailable negative markers",
  {
    unavailable_markers <- list(
      "T cell" = c(
        "NotAGene1",
        "NotAGene2"
      ),
      "B cell" = c(
        "NotAGene3",
        "NotAGene4"
      )
    )

    scores <- negative_marker_score(
      toy_seurat,
      unavailable_markers
    )

    expect_type(
      scores,
      "double"
    )

    expect_length(
      scores,
      ncol(toy_seurat)
    )

    expect_true(
      all(scores == 0)
    )
  }
)


test_that(
  "negative_marker_score returns zero for unmatched labels",
  {
    tmp <- toy_seurat

    tmp$predicted_label <- "UnknownCell"

    scores <- negative_marker_score(
      tmp,
      toy_negative_markers
    )

    expect_type(
      scores,
      "double"
    )

    expect_length(
      scores,
      ncol(tmp)
    )

    expect_true(
      all(scores == 0)
    )
  }
)
