test_that("beta_dissimilarity_plot returns a ggplot for the shared partition", {
  toy <- make_beta_toy_community()

  p <- beta_dissimilarity_plot(
    toy$table, toy$metadata,
    condition1_col = "Group",
    partition      = "shared",
    family         = "sorensen"
  )

  expect_s3_class(p, "ggplot")
})

test_that("beta_dissimilarity_plot facets by condition2_col and adds a stat comparison", {
  toy <- make_beta_toy_community()

  p <- beta_dissimilarity_plot(
    toy$table, toy$metadata,
    condition1_col = "Group",
    condition2_col = "Batch",
    partition      = "turnover",
    family         = "jaccard",
    stat           = "kruskal.test"
  )

  expect_s3_class(p, "ggplot")
})

test_that("beta_dissimilarity_plot's comparison_condition1 filter keeps only the requested pairs", {
  toy <- make_beta_toy_community()

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  beta_dissimilarity_plot(
    toy$table, toy$metadata,
    comparison_condition1 = c("A_vs_B"),
    condition1_col         = "Group",
    partition              = "shared",
    save_table             = TRUE,
    table_filename         = tmp
  )

  saved <- utils::read.delim(tmp, check.names = FALSE)
  expect_true(all(saved$condition1_group == "A_vs_B"))
})

test_that("beta_dissimilarity_plot keeps each pair of samples once and no self-comparisons", {
  toy <- make_beta_toy_community()
  n <- sum(vapply(toy$table, is.numeric, logical(1)))

  for (part in c("shared", "turnover")) {
    tmp <- tempfile(fileext = ".txt")
    beta_dissimilarity_plot(
      toy$table, toy$metadata,
      condition1_col = "Group",
      partition      = part,
      save_table     = TRUE,
      table_filename = tmp
    )
    saved <- utils::read.delim(tmp, check.names = FALSE)
    unlink(tmp)

    expect_equal(nrow(saved), choose(n, 2))
    expect_false(any(saved$site1 == saved$site2))
    pair_key <- paste(pmin(saved$site1, saved$site2), pmax(saved$site1, saved$site2))
    expect_false(any(duplicated(pair_key)))
  }
})
test_that("beta_dissimilarity_plot facets keep only pairs within the same condition2 value", {
  toy <- make_beta_toy_community()

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  beta_dissimilarity_plot(
    toy$table, toy$metadata,
    condition1_col = "Group",
    condition2_col = "Batch",
    partition      = "shared",
    save_table     = TRUE,
    table_filename = tmp
  )

  saved <- utils::read.delim(tmp, check.names = FALSE)
  expect_gt(nrow(saved), 0)
  expect_true(all(saved$Batch.x == saved$Batch.y))
})