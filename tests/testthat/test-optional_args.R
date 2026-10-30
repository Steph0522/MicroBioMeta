load_bacteria <- function() {
    list(
        table = read.delim(system.file("extdata", "subset_table_bacteria.txt", package = "MicroBioMeta"),
            row.names = 1, check.names = FALSE
        ),
        metadata = read.delim(system.file("extdata", "subset_metadata_bacteria.txt", package = "MicroBioMeta"),
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

test_that("abundance_bar_plot taxa_order flips the stack, not the colors", {
    table <- read.delim(system.file("extdata", "subset_table_bacteria.txt", package = "MicroBioMeta"),
        row.names = 1, check.names = FALSE
    )
    metadata <- read.delim(system.file("extdata", "subset_metadata_bacteria.txt", package = "MicroBioMeta"),
        check.names = FALSE
    )
    for (pal in c("default", "safe")) {
        p <- lapply(c("ascending", "descending"), function(o) {
            abundance_bar_plot(table, metadata,
                level = "phylum", x_col = "Location", top_n = 5,
                add_remained = TRUE, palette = pal, taxa_order = o
            )
        })
        lev <- lapply(p, function(x) levels(x$data$taxonomy))
        expect_equal(lev[[2]], rev(lev[[1]]))
        cols <- lapply(p, function(x) x$scales$get_scales("fill")$palette(length(lev[[1]])))
        expect_equal(cols[[1]][lev[[1]]], cols[[2]][lev[[1]]])
    }
    expect_error(abundance_bar_plot(table, metadata, taxa_order = "x"))
})
