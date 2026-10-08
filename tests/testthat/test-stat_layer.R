test_that(".mbm_p_label writes 'p<0.001' for small values and 'p=x' otherwise", {
    expect_equal(.mbm_p_label(0.00001), "p<0.001")
    expect_equal(.mbm_p_label(0.0123), "p=0.012")
})

test_that("a two-sample stat compares every pair of boxes within each facet", {
    toy <- make_beta_toy_community()

    p <- beta_dissimilarity_plot(
        toy$table, toy$metadata,
        condition1_col = "Group",
        partition = "turnover",
        stat = "wilcox.test"
    )
    built <- ggplot2::ggplot_build(p)
    pw <- unique(built$data[[2]][, c("PANEL", "group1", "group2")])

    # Non-significant pairs are hidden (hide.ns), so up to choose(n, 2) of
    # them are drawn, each one a different pair of boxes with its own bracket.
    n_boxes <- length(unique(p$data$condition1_group))
    expect_gt(nrow(pw), 1)
    expect_lte(nrow(pw), choose(n_boxes, 2))
    expect_false(any(duplicated(pw[, c("group1", "group2")])))
})

test_that("a global stat gives one p-value per panel", {
    toy <- make_beta_toy_community()

    p <- beta_dissimilarity_plot(
        toy$table, toy$metadata,
        condition1_col = "Group",
        partition = "turnover",
        stat = "kruskal.test"
    )
    built <- ggplot2::ggplot_build(p)
    expect_equal(nrow(built$data[[2]]), 1)
})
