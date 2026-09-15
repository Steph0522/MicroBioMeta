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
