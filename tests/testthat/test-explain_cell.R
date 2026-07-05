test_that("explain_cell returns interpretation", {

  obj <- toy_seurat

  obj$confidence_score <- 0.8
  obj$marker_score <- 0.7
  obj$neighbor_score <- 0.9
  obj$entropy_norm <- 0.2

  out <- explain_cell(
    obj,
    "Cell1"
  )

  expect_type(out, "character")
  expect_true(length(out) >= 1)

})

test_that("explain_cell errors for unknown cell", {

  expect_error(

    explain_cell(
      toy_seurat,
      "UnknownCell"
    )

  )

})
