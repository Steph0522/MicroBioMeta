test_that("beta_ord_plot returns a ggplot for PCA on compositional (Aitchison) distance", {
  toy <- make_toy_community()

  p <- beta_ord_plot(
    toy$table, toy$metadata,
    distance   = "compositional",
    ordination = "PCA",
    group_col  = "Group",
    seed       = 1
  )

  expect_s3_class(p, "ggplot")
})

test_that("beta_ord_plot rejects PCA combined with a non-compositional distance", {
  toy <- make_toy_community()

  # PCA is only defined on the CLR-transformed (compositional) data here;
  # asking for it with a raw ecological distance should error, not silently
  # produce a misleading plot.
  expect_error(
    beta_ord_plot(
      toy$table, toy$metadata,
      distance   = "bray",
      ordination = "PCA",
      group_col  = "Group"
    ),
    "PCA is only available"
  )
})

test_that("beta_ord_plot is reproducible across runs with the same seed", {
  toy <- make_toy_community()

  p1 <- beta_ord_plot(toy$table, toy$metadata, distance = "compositional",
                      ordination = "PCA", group_col = "Group", seed = 42)
  p2 <- beta_ord_plot(toy$table, toy$metadata, distance = "compositional",
                      ordination = "PCA", group_col = "Group", seed = 42)

  # Same seed -> the CLR Monte-Carlo draw (and thus the plotted coordinates)
  # must match exactly.
  expect_equal(p1$data, p2$data)
})
