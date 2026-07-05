test_that("marker_score returns numeric vector for Seurat", {

  scores <- marker_score(
    toy_seurat,
    toy_markers
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

test_that("marker_score returns numeric vector for SingleCellExperiment", {

  scores <- marker_score(
    toy_sce,
    toy_markers
  )

  expect_type(scores, "double")

  expect_length(
    scores,
    ncol(toy_sce)
  )

})

test_that("marker_score errors when label column is missing", {

  expect_error(

    marker_score(
      toy_seurat,
      toy_markers,
      label_column = "unknown"
    )

  )

})

test_that("marker_score warns for unknown labels", {

  tmp <- toy_seurat

  tmp$predicted_label <- "UnknownCell"

  expect_warning(

    marker_score(
      tmp,
      toy_markers
    )

  )

})
