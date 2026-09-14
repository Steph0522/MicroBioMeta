test_that("aldex_heatmap_plot returns a HeatmapList for a two-group comparison", {
  toy <- make_toy_community()

  # A low effect_threshold, plus the fixture's planted group signal
  # (features 1-4 boosted in group A, 5-8 in group B), makes it very likely
  # at least one taxon passes the filter regardless of the Monte Carlo
  # draw - set.seed() makes that draw reproducible either way.
  set.seed(1)
  ht <- aldex_heatmap_plot(
    table            = toy$table,
    metadata         = toy$metadata,
    col_cond         = "Group",
    effect_threshold = 0.1
  )

  expect_s4_class(ht, "HeatmapList")
})

test_that("aldex_heatmap_plot errors when col_cond doesn't have exactly two groups", {
  toy <- make_toy_community()
  toy$metadata$ThreeGroups <- rep(c("A", "B", "C"), 2)

  expect_error(
    aldex_heatmap_plot(
      table    = toy$table,
      metadata = toy$metadata,
      col_cond = "ThreeGroups"
    ),
    "Exactly two conditions"
  )
})

test_that("aldex_heatmap_plot saves the filtered results table when requested", {
  toy <- make_toy_community()

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  set.seed(1)
  aldex_heatmap_plot(
    table            = toy$table,
    metadata         = toy$metadata,
    col_cond         = "Group",
    effect_threshold = 0.1,
    save_table       = TRUE,
    table_filename   = tmp
  )

  expect_true(file.exists(tmp))
})
