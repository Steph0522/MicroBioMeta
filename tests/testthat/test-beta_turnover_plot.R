test_that("beta_turnover_plot returns a plot comparing beta diversity within same-batch pairs", {
  toy <- make_beta_toy_community()

  # comparison_condition1/2 match regardless of order, so each pair only
  # needs to be listed once - not also its "B_vs_A" reverse.
  result <- beta_turnover_plot(
    table                 = toy$table,
    metadata              = toy$metadata,
    comparison_condition1 = c("A_vs_A", "A_vs_B", "B_vs_B"),
    comparison_condition2 = c("B1_vs_B1", "B1_vs_B2", "B2_vs_B2"),
    condition1.x          = "Group.x",
    condition1.y          = "Group.y",
    condition2.x          = "Batch.x",
    condition2.y          = "Batch.y",
    color_facets_x        = c("#E69F00", "#56B4E9"),
    color_axis_x          = c("A" = "#E69F00", "B" = "#56B4E9")
  )

  expect_s3_class(result, "ggplot")
})

test_that("beta_turnover_plot saves a table with the Hill-order column when requested", {
  toy <- make_beta_toy_community()

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  beta_turnover_plot(
    table                 = toy$table,
    metadata              = toy$metadata,
    comparison_condition1 = c("A_vs_A", "A_vs_B", "B_vs_A", "B_vs_B"),
    comparison_condition2 = c("B1_vs_B1", "B1_vs_B2", "B2_vs_B1", "B2_vs_B2"),
    condition1.x          = "Group.x",
    condition1.y          = "Group.y",
    condition2.x          = "Batch.x",
    condition2.y          = "Batch.y",
    color_facets_x        = c("#E69F00", "#56B4E9"),
    color_axis_x          = c("A" = "#E69F00", "B" = "#56B4E9"),
    save_table            = TRUE,
    table_filename        = tmp
  )

  expect_true(file.exists(tmp))
  saved <- utils::read.delim(tmp, check.names = FALSE)
  expect_true("orden" %in% names(saved))
  expect_true(all(saved$orden %in% c("q0", "q1", "q2")))
})
