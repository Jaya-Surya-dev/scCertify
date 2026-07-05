test_that(
  "classify_confidence works with default thresholds",
  {

    scores <- c(
      0.2,
      0.6,
      0.9
    )

    result <- classify_confidence(
      scores
    )

    expect_length(
      result,
      3
    )

    expect_type(
      result,
      "character"
    )

    expect_equal(
      result,
      c(
        "Low",
        "Moderate",
        "High"
      )
    )

  }
)

test_that(
  "classify_confidence supports custom thresholds",
  {

    scores <- c(
      0.3,
      0.6,
      0.9
    )

    result <- classify_confidence(
      scores,
      moderate_threshold = 0.4,
      high_threshold = 0.8
    )

    expect_equal(
      result,
      c(
        "Low",
        "Moderate",
        "High"
      )
    )

  }
)

test_that(
  "classify_confidence validates thresholds",
  {

    expect_error(

      classify_confidence(
        c(0.5),
        moderate_threshold = 0.8,
        high_threshold = 0.5
      )

    )

  }
)

test_that(
  "classify_confidence requires numeric scores",
  {

    expect_error(

      classify_confidence(
        c("a", "b")
      )

    )

  }
)
