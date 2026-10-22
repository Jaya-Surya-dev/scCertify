test_that(
  "discovery_aware returns results for SingleCellExperiment",
  {
    result <- discovery_aware(
      object = toy_sce,
      markers = toy_markers,
      negative_markers = toy_negative_markers
    )

    expect_true(
      is.data.frame(result)
    )

    expect_equal(
      nrow(result),
      ncol(toy_sce)
    )

    expect_true(
      all(
        c(
          "positive_marker_score",
          "negative_marker_score",
          "discovery_status"
        ) %in% colnames(result)
      )
    )

    expect_type(
      result$positive_marker_score,
      "double"
    )

    expect_type(
      result$negative_marker_score,
      "double"
    )

    expect_type(
      result$discovery_status,
      "character"
    )
  }
)


test_that(
  "discovery_aware returns results for Seurat",
  {
    result <- discovery_aware(
      object = toy_seurat,
      markers = toy_markers,
      negative_markers = toy_negative_markers
    )

    expect_true(
      is.data.frame(result)
    )

    expect_equal(
      nrow(result),
      ncol(toy_seurat)
    )

    expect_true(
      all(
        c(
          "positive_marker_score",
          "negative_marker_score",
          "discovery_status"
        ) %in% colnames(result)
      )
    )

    expect_type(
      result$positive_marker_score,
      "double"
    )

    expect_type(
      result$negative_marker_score,
      "double"
    )

    expect_type(
      result$discovery_status,
      "character"
    )
  }
)


test_that(
  "discovery_aware preserves cell names",
  {
    result <- discovery_aware(
      object = toy_sce,
      markers = toy_markers,
      negative_markers = toy_negative_markers
    )

    expect_identical(
      rownames(result),
      colnames(toy_sce)
    )
  }
)


test_that(
  "discovery_aware accepts custom thresholds",
  {
    result <- discovery_aware(
      object = toy_sce,
      markers = toy_markers,
      negative_markers = toy_negative_markers,
      positive_threshold = 0.70,
      negative_threshold = 0.70
    )

    expect_true(
      is.data.frame(result)
    )

    expect_equal(
      nrow(result),
      ncol(toy_sce)
    )
  }
)


test_that(
  "discovery_aware returns valid discovery statuses",
  {
    result <- discovery_aware(
      object = toy_sce,
      markers = toy_markers,
      negative_markers = toy_negative_markers
    )

    expect_true(
      all(
        result$discovery_status %in% c(
          "Known",
          "Possible novel state",
          "Possible transitional state",
          "Insufficient evidence"
        )
      )
    )
  }
)
