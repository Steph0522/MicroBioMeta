test_that("y_axis_title changes the y-axis title of the beta boxplots", {
  toy <- make_beta_toy_community()

  p1 <- beta_dissimilarity_plot(toy$table, toy$metadata, condition1_col = "Group",
                                partition = "shared", y_axis_title = "Shared ASVs")
  expect_equal(p1$labels$y, "Shared ASVs")

  # default stays the automatic title
  p2 <- beta_dissimilarity_plot(toy$table, toy$metadata, condition1_col = "Group",
                                partition = "shared")
  expect_equal(p2$labels$y, "Beta diversity (shared features)")

  p3 <- beta_turnover_plot(toy$table, toy$metadata,
                           comparison_condition1 = c("A_vs_B"),
                           condition1_col = "Group", y_axis_title = "Turnover")
  expect_equal(p3$labels$y, "Turnover")
})
