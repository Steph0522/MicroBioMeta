load_bacteria <- function() {
    list(
        table = read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
            row.names = 1, check.names = FALSE
        ),
        metadata = read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
            check.names = FALSE
        )
    )
}

test_that("beta_ord_plot works without group_col", {
    d <- load_bacteria()
    p <- suppressMessages(beta_ord_plot(d$table, d$metadata, distance = "bray", ordination = "PCoA"))
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
    p2 <- suppressMessages(beta_ord_plot(d$table, d$metadata,
        distance = "bray",
        ordination = "PCoA", shape_col = "Treatment"
    ))
    expect_s3_class(p2, "ggplot")
})

test_that("functions that need groups say so clearly", {
    d <- load_bacteria()
    expect_error(beta_partition_ord_plot(d$table, d$metadata), "`group_col` is required")
    expect_error(venn_plot(d$table, d$metadata), "merge_by")
})

test_that("abundance_heatmap_plot has a default top_n", {
    d <- load_bacteria()
    expect_equal(formals(abundance_heatmap_plot)$top_n, 15)
})
