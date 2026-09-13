test_that(
  "de_novo_signatures returns signatures for Seurat",
  {

    signatures <- de_novo_signatures(
      toy_seurat,
      cluster_column = "cluster",
      top_n = 5
    )

    expect_type(
      signatures,
      "list"
    )

    expect_length(
      signatures,
      2
    )

    expect_true(
      all(
        vapply(
          signatures,
          length,
          integer(1)
        ) == 5
      )
    )
  }
)


test_that(
  "de_novo_signatures returns signatures for SingleCellExperiment",
  {

    signatures <- de_novo_signatures(
      toy_sce,
      cluster_column = "cluster",
      top_n = 5
    )

    expect_type(
      signatures,
      "list"
    )

    expect_length(
      signatures,
      2
    )

    expect_true(
      all(
        vapply(
          signatures,
          length,
          integer(1)
        ) == 5
      )
    )
  }
)


test_that(
  "de_novo_signatures errors when cluster column is missing",
  {

    expect_error(
      de_novo_signatures(
        toy_seurat,
        cluster_column = "unknown"
      )
    )
  }
)


test_that(
  "de_novo_signatures errors when top_n is invalid",
  {

    expect_error(
      de_novo_signatures(
        toy_seurat,
        cluster_column = "cluster",
        top_n = 0
      )
    )
  }
)


test_that(
  "de_novo_signatures errors when cluster identities contain NA",
  {

    tmp <- toy_seurat

    tmp$cluster[1] <- NA

    expect_error(
      de_novo_signatures(
        tmp,
        cluster_column = "cluster"
      )
    )
  }
)
