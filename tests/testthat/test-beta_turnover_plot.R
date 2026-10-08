test_that("beta_turnover_plot returns a plot comparing beta diversity within same-batch pairs", {
    toy <- make_beta_toy_community()

    # comparison_condition1/2 match regardless of order, so each pair only
    # needs to be listed once - not also its "B_vs_A" reverse.
    result <- beta_turnover_plot(
        table = toy$table,
        metadata = toy$metadata,
        comparison_condition1 = c("A_vs_A", "A_vs_B", "B_vs_B"),
        comparison_condition2 = c("B1_vs_B1", "B1_vs_B2", "B2_vs_B2"),
        condition1_col = "Group",
        condition2_col = "Batch",
        facet_colors = c("#E69F00", "#56B4E9"),
        group_colors = c("A_vs_A" = "#E69F00", "A_vs_B" = "#56B4E9", "B_vs_B" = "#009E73")
    )

    expect_s3_class(result, "ggplot")
})

test_that("beta_turnover_plot saves a table with the Hill-order column when requested", {
    toy <- make_beta_toy_community()

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)

    beta_turnover_plot(
        table = toy$table,
        metadata = toy$metadata,
        comparison_condition1 = c("A_vs_A", "A_vs_B", "B_vs_A", "B_vs_B"),
        comparison_condition2 = c("B1_vs_B1", "B1_vs_B2", "B2_vs_B1", "B2_vs_B2"),
        condition1_col = "Group",
        condition2_col = "Batch",
        facet_colors = c("#E69F00", "#56B4E9"),
        group_colors = c("A_vs_A" = "#E69F00", "A_vs_B" = "#56B4E9", "B_vs_B" = "#009E73"),
        save_table = TRUE,
        table_filename = tmp
    )

    expect_true(file.exists(tmp))
    saved <- utils::read.delim(tmp, check.names = FALSE)
    expect_true("orden" %in% names(saved))
    expect_true(all(saved$orden %in% c("q0", "q1", "q2")))
})

test_that("beta_turnover_plot facets by condition2_col alone, keeping same-value pairs", {
    toy <- make_beta_toy_community()

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)

    p <- beta_turnover_plot(
        table                 = toy$table,
        metadata              = toy$metadata,
        comparison_condition1 = c("A_vs_A", "A_vs_B", "B_vs_B"),
        condition1_col        = "Group",
        condition2_col        = "Batch",
        save_table            = TRUE,
        table_filename        = tmp
    )

    expect_s3_class(p, "ggplot")
    expect_setequal(levels(p$data$.facet2_col), c("B1", "B2"))
    saved <- utils::read.delim(tmp, check.names = FALSE)
    expect_true(all(saved$Batch.x == saved$Batch.y))
})
test_that("beta_turnover_plot warns when group_colors names match no comparison", {
    toy <- make_beta_toy_community()
    expect_warning(
        beta_turnover_plot(toy$table, toy$metadata,
            comparison_condition1 = c("A_vs_B"), condition1_col = "Group",
            group_colors = c("A" = "#E69F00", "B" = "#56B4E9")
        ),
        "drawn grey"
    )
})
