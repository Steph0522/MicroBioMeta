test_that("beta_test_table runs a compositional PERMANOVA and returns a table figure", {
  toy <- make_toy_community()

  res <- beta_test_table(
    toy$table, toy$metadata,
    formula_str  = "Group",
    method       = "compositional",
    test         = "permanova",
    permutations = 99,
    seed         = 1
  )

  # The result is the formatted results table rendered as a ggplot figure
  # (via ggpubr::ggtexttable); the numeric table itself is exposed through
  # save_table (checked below).
  expect_s3_class(res, "ggplot")
})

test_that("beta_test_table saves a PERMANOVA table with the tested term, R2 and p-value", {
  toy <- make_toy_community()

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  beta_test_table(
    toy$table, toy$metadata,
    formula_str    = "Group",
    method         = "compositional",
    test           = "permanova",
    permutations   = 99,
    seed           = 1,
    save_table     = TRUE,
    table_filename = tmp
  )

  expect_true(file.exists(tmp))

  saved <- utils::read.delim(tmp, check.names = FALSE)
  expect_true("Term" %in% names(saved))
  expect_true("Group" %in% saved$Term)
  expect_true(any(grepl("R2", names(saved))))
  expect_true(any(grepl("Pr", names(saved))))
})

test_that("beta_test_table also accepts a plain ecological distance (bray)", {
  toy <- make_toy_community()

  res <- beta_test_table(
    toy$table, toy$metadata,
    formula_str  = "Group",
    method       = "bray",
    test         = "permanova",
    permutations = 99
  )

  expect_s3_class(res, "ggplot")
})
