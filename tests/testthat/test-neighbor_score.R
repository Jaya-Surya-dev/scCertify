test_that(
  "neighbor_score returns numeric vector for Seurat",
  {

    scores <- neighbor_score(
      toy_seurat
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
      all(scores >= 0)
    )

    expect_true(
      all(scores <= 1)
    )

  }
)

test_that(
  "neighbor_score returns numeric vector for SingleCellExperiment",
  {

    scores <- neighbor_score(
      toy_sce
    )

    expect_type(
      scores,
      "double"
    )

    expect_length(
      scores,
      ncol(toy_sce)
    )

    expect_true(
      all(scores >= 0)
    )

    expect_true(
      all(scores <= 1)
    )

  }
)

test_that(
  "neighbor_score errors for missing label column",
  {

    expect_error(

      neighbor_score(
        toy_seurat,
        label_column = "unknown"
      )

    )

  }
)

test_that(
  "neighbor_score respects k parameter",
  {

    score1 <- neighbor_score(
      toy_seurat,
      k = 3
    )

    score2 <- neighbor_score(
      toy_seurat,
      k = 5
    )

    expect_length(
      score1,
      ncol(toy_seurat)
    )

    expect_length(
      score2,
      ncol(toy_seurat)
    )

  }
)
