test_that("beta_partition_ord_plot returns a combined ordination panel", {
  toy <- make_beta_toy_community()

  # The function auto-detects and drops a taxonomy column on its own (as
  # documented), so the fixture's table is passed as-is.
  p <- beta_partition_ord_plot(
    table     = toy$table,
    metadata  = toy$metadata,
    index     = "jaccard",
    group_col = "Group"
  )

  expect_s3_class(p, "ggplot")
})

test_that("beta_partition_ord_plot saves jaccard/turnover/nestedness tables when requested", {
  toy <- make_beta_toy_community()

  base <- tempfile()
  on.exit(unlink(paste0(base, c("_jacs.txt", "_jtus.txt", "_jnes.txt"))), add = TRUE)

  beta_partition_ord_plot(
    table          = toy$table,
    metadata       = toy$metadata,
    index          = "jaccard",
    group_col      = "Group",
    save_table     = TRUE,
    table_filename = base
  )

  expect_true(file.exists(paste0(base, "_jacs.txt")))
  expect_true(file.exists(paste0(base, "_jtus.txt")))
  expect_true(file.exists(paste0(base, "_jnes.txt")))
})
