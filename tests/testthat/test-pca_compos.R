test_that("beta_ord_plot returns a ggplot for PCA on compositional (Aitchison) distance", {
    toy <- make_toy_community()

    p <- beta_ord_plot(
        toy$table, toy$metadata,
        distance = "compositional",
        ordination = "PCA",
        group_col = "Group"
    )

    expect_s3_class(p, "ggplot")
})

test_that("beta_ord_plot rejects PCA combined with a non-compositional distance", {
    toy <- make_toy_community()

    expect_error(
        beta_ord_plot(
            toy$table, toy$metadata,
            distance = "bray",
            ordination = "PCA",
            group_col = "Group"
        ),
        "PCA is only available"
    )
})

test_that("beta_ord_plot is reproducible across runs with the same set.seed()", {
    toy <- make_toy_community()

    set.seed(42)
    p1 <- beta_ord_plot(toy$table, toy$metadata,
        distance = "compositional",
        ordination = "PCA", group_col = "Group"
    )
    set.seed(42)
    p2 <- beta_ord_plot(toy$table, toy$metadata,
        distance = "compositional",
        ordination = "PCA", group_col = "Group"
    )


    expect_equal(p1$data, p2$data)
})
