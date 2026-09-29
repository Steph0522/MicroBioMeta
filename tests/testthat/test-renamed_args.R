test_that(".mbm_renamed_args maps an old name to the new one with a warning", {
  expect_warning(
    out <- .mbm_renamed_args(list(col_cond = "Group"), c(col_cond = "group_col"), "f"),
    "`col_cond` was renamed to `group_col`"
  )
  expect_equal(out, list(group_col = "Group"))
})

test_that(".mbm_renamed_args rejects arguments that don't exist", {
  expect_error(
    .mbm_renamed_args(list(grop_col = "Group"), c(col_cond = "group_col"), "f"),
    "Unused argument"
  )
})

test_that("old argument names still work in the plotting functions", {
  toy <- make_beta_toy_community()

  expect_warning(
    p <- beta_partition_ord_plot(table = toy$table, metadata = toy$metadata,
                                 index = "jaccard", group_col = "Group"),
    "`index` was renamed to `family`"
  )
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))

  toy2 <- make_toy_community()
  expect_warning(
    beta_test_table(toy2$table, toy2$metadata, formula_str = "Group",
                    method = "bray", permutations = 99),
    "`method` was renamed to `distance`"
  )
})
