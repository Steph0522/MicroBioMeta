test_that("aldex_volcano_plot returns a ggplot for the effect-size view", {
  toy <- make_toy_community()

  # aldex()'s Monte Carlo Dirichlet sampling uses R's global RNG; seeding
  # here (the function itself doesn't) makes the test deterministic.
  set.seed(1)
  p <- aldex_volcano_plot(
    table    = toy$table,
    metadata = toy$metadata,
    col_cond = "Group",
    type     = "effect",
    cond     = "A"
  )

  expect_s3_class(p, "ggplot")
})

test_that("aldex_volcano_plot returns a ggplot for the volcano view", {
  toy <- make_toy_community()

  set.seed(1)
  p <- aldex_volcano_plot(
    table    = toy$table,
    metadata = toy$metadata,
    col_cond = "Group",
    type     = "volcano",
    cond     = "A"
  )

  expect_s3_class(p, "ggplot")
})

test_that("aldex_volcano_plot saves the ALDEx2 result table when requested", {
  toy <- make_toy_community()

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(1)
  aldex_volcano_plot(
    table          = toy$table,
    metadata       = toy$metadata,
    col_cond       = "Group",
    type           = "effect",
    cond           = "A",
    save_table     = TRUE,
    table_filename = tmp
  )

  expect_true(file.exists(tmp))
  saved <- utils::read.delim(tmp, check.names = FALSE)
  expect_true(all(c("effect", "wi.ep", "wi.eBH") %in% names(saved)))
})

test_that("aldex_volcano_plot rejects an invalid type", {
  toy <- make_toy_community()

  expect_error(
    aldex_volcano_plot(
      table    = toy$table,
      metadata = toy$metadata,
      col_cond = "Group",
      type     = "not-a-type"
    ),
    "type must be either"
  )
})
