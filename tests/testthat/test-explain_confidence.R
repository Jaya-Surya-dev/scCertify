test_that("explain_confidence returns reasons", {

  obj <- toy_seurat

  obj$marker_score <- 0.8
  obj$neighbor_score <- 0.9
  obj$entropy_norm <- 0.2
  obj$doublet_score <- 0.1

  out <- explain_confidence(
    obj,
    "Cell1"
  )

  expect_type(out, "character")
  expect_true(length(out) >= 1)

})

test_that("explain_confidence errors for unknown cell", {

  expect_error(

    explain_confidence(
      toy_seurat,
      "UnknownCell"
    )

  )

})
