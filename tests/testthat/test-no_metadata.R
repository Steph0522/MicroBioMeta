test_that("abundance_bar_plot draws one bar per sample without metadata", {
    toy <- make_toy_community()
    p <- abundance_bar_plot(toy$table, level = "phylum")
    expect_s3_class(p, "ggplot")
    expect_setequal(as.character(unique(p$data$SAMPLEID)), toy$metadata$SampleID)
    expect_error(abundance_bar_plot(toy$table, facet_by = "Group"), "needs `metadata`")
})

test_that("abundance_bar_plot without x_col uses one bar per sample", {
    toy <- make_toy_community()
    p <- abundance_bar_plot(toy$table, toy$metadata, level = "phylum")
    expect_setequal(as.character(unique(p$data$SAMPLEID)), toy$metadata$SampleID)
    expect_error(abundance_bar_plot(toy$table, toy$metadata, x_col = "Nope"), "not a column")
})

test_that("abundance_heatmap_plot works without metadata", {
    toy <- make_toy_community()
    expect_s3_class(
        suppressWarnings(abundance_heatmap_plot(toy$table, top_n = 5, draw = FALSE)),
        "mbm_heatmap"
    )
    expect_error(abundance_heatmap_plot(toy$table, condition1 = "Group"), "need `metadata`")
})

test_that("a function passed as `table` gives a clear error", {
    expect_error(abundance_bar_plot(table = table), "must be a data frame.*a function")
    expect_error(abundance_heatmap_plot(table = table), "must be a data frame")
})
