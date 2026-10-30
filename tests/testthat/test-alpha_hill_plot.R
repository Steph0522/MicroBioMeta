hill_values <- function(...) {
    f <- tempfile(fileext = ".txt")
    suppressMessages(alpha_hill_plot(..., save_table = TRUE, table_filename = f))
    utils::read.delim(f)[, c("q0", "q1", "q2")]
}

comm <- data.frame(
    s1 = c(5, 3, 2, 0), s2 = c(1, 1, 4, 4), s3 = c(0, 6, 1, 3),
    row.names = paste0("t", 1:4)
)
meta <- data.frame(SAMPLEID = c("s1", "s2", "s3"), group = c("a", "a", "b"))
args <- list(table = comm, metadata = meta, x_col = "group", fill_col = "group")

test_that("equidistant taxa make phylogenetic and functional Hill numbers taxonomic", {
    taxonomic <- do.call(hill_values, args)
    star <- ape::read.tree(text = "(t1:1,t2:1,t3:1,t4:1,t5:1);")
    traits <- as.data.frame(diag(5), row.names = paste0("t", 1:5))

    phylo <- do.call(hill_values, c(args, list(diversity = "phylogenetic", tree = star)))
    func <- do.call(hill_values, c(args, list(diversity = "functional", traits = traits)))

    expect_equal(phylo, taxonomic, tolerance = 1e-8)
    expect_equal(func, taxonomic, tolerance = 1e-8)
})

test_that("diversity needs its input and a valid name", {
    expect_error(do.call(alpha_hill_plot, c(args, list(diversity = "phylogenetic"))), "needs `tree`")
    expect_error(do.call(alpha_hill_plot, c(args, list(diversity = "functional"))), "needs `traits`")
    expect_error(do.call(alpha_hill_plot, c(args, list(diversity = "beta"))), "Invalid diversity")
})

test_that("taxa missing from the tree or traits are dropped with a warning", {
    short <- ape::read.tree(text = "(t1:1,t2:1,t3:1);")
    expect_warning(
        do.call(hill_values, c(args, list(diversity = "phylogenetic", tree = short))),
        "1 of 4 taxa dropped"
    )
    other <- data.frame(a = 1:2, row.names = c("x", "y"))
    expect_error(
        do.call(hill_values, c(args, list(diversity = "functional", traits = other))),
        "No matching taxa"
    )
})
