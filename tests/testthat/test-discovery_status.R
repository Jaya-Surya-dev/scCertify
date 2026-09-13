test_that(
  "discovery_status identifies known annotations",
  {

    result <- discovery_status(
      marker_score = 0.8,
      negative_marker_score = 0.2
    )

    expect_identical(
      result,
      "Known"
    )
  }
)


test_that(
  "discovery_status identifies possible novel states",
  {

    result <- discovery_status(
      marker_score = 0.2,
      negative_marker_score = 0.2
    )

    expect_identical(
      result,
      "Possible novel state"
    )
  }
)


test_that(
  "discovery_status identifies possible transitional states",
  {

    result <- discovery_status(
      marker_score = 0.8,
      negative_marker_score = 0.8
    )

    expect_identical(
      result,
      "Possible transitional state"
    )
  }
)


test_that(
  "discovery_status identifies insufficient evidence",
  {

    result <- discovery_status(
      marker_score = 0.2,
      negative_marker_score = 0.8
    )

    expect_identical(
      result,
      "Insufficient evidence"
    )
  }
)


test_that(
  "discovery_status rejects vectors with different lengths",
  {

    expect_error(
      discovery_status(
        marker_score = c(0.8, 0.7),
        negative_marker_score = 0.2
      )
    )
  }
)


test_that(
  "discovery_status rejects NA scores",
  {

    expect_error(
      discovery_status(
        marker_score = NA,
        negative_marker_score = 0.2
      )
    )
  }
)
