test_that("aldex_heatmap_plot returns a printable grob for a two-group comparison", {
    toy <- make_toy_community()

    set.seed(1)
    ht <- aldex_heatmap_plot(
        table = toy$table,
        metadata = toy$metadata,
        group_col = "Group",
        effect_threshold = 0.1,
        pval_threshold = NULL
    )

    expect_s3_class(ht, "gTree")
    expect_s3_class(ht, "mbm_heatmap")
    expect_no_error(print(ht))
})

test_that("aldex_heatmap_plot errors when group_col doesn't have exactly two groups", {
    toy <- make_toy_community()
    toy$metadata$ThreeGroups <- rep(c("A", "B", "C"), 2)

    expect_error(
        aldex_heatmap_plot(
            table = toy$table,
            metadata = toy$metadata,
            group_col = "ThreeGroups"
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
        table = toy$table,
        metadata = toy$metadata,
        group_col = "Group",
        effect_threshold = 0.1,
        pval_threshold = NULL,
        save_table = TRUE,
        table_filename = tmp
    )

    expect_true(file.exists(tmp))
})

test_that("aldex_heatmap_plot accepts plain color vectors for effect and p-value strips", {
    toy <- make_toy_community()
    set.seed(1)
    expect_s3_class(
        aldex_heatmap_plot(
            table = toy$table, metadata = toy$metadata, group_col = "Group",
            effect_threshold = 0.1, pval_threshold = NULL, draw = FALSE,
            effect_colors = c("#0072B2", "white", "#E69F00"),
            pvalue_colors = c(
                "<0.001" = "black", "<0.01" = "grey30",
                "<0.05" = "grey60", ">0.05" = "grey90"
            )
        ),
        "mbm_heatmap"
    )
    expect_error(
        aldex_heatmap_plot(
            table = toy$table, metadata = toy$metadata,
            group_col = "Group", effect_threshold = 0.1,
            draw = FALSE, effect_colors = c("red", "blue")
        ),
        "3 colors"
    )
})
