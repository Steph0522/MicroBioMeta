toy <- make_beta_toy_community()
counts <- as.matrix(toy$table[, colnames(toy$table) != "taxonomy"])
tree <- ape::read.tree(text = paste0(
    "((F1:1,F2:1):1,(F3:1,F4:1):2,((F5:1,F6:1):1,(F7:1,F8:1):1):1,",
    "(F9:1,(F10:1,(F11:1,F12:1):1):1):1);"
))
traits <- data.frame(a = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), row.names = rownames(counts))

test_that("U and C are Jaccard and Sorensen at q = 0", {
    expect_equal(
        as.vector(.mbm_hill_dist(counts, "U", 0, "taxonomic", NULL, NULL)),
        as.vector(vegan::vegdist(t(counts), "jaccard", binary = TRUE))
    )
    expect_equal(
        as.vector(.mbm_hill_dist(counts, "C", 0, "taxonomic", NULL, NULL)),
        as.vector(vegan::vegdist(t(counts), "bray", binary = TRUE))
    )
})

test_that("the Hill metrics reach the ordination and the tests", {
    p <- suppressMessages(beta_ord_plot(toy$table, toy$metadata,
        distance = "U", ordination = "PCoA", group_col = "Group",
        q = 0, diversity = "phylogenetic", tree = tree
    ))
    expect_s3_class(p, "ggplot")

    test_values <- function(...) {
        file <- tempfile(fileext = ".txt")
        set.seed(1)
        suppressWarnings(suppressMessages(beta_test_table(..., save_table = TRUE, table_filename = file)))
        utils::read.delim(file)
    }
    direct <- vegan::vegdist(t(counts), "jaccard", binary = TRUE)
    expect_equal(
        test_values(toy$table, toy$metadata, "Group", distance = "U", q = 0),
        test_values(direct, toy$metadata, "Group")
    )

    disp <- test_values(toy$table, toy$metadata, "Group",
        distance = "U", test = "betadisper",
        diversity = "functional", traits = traits
    )
    expect_gt(nrow(disp), 0)
    expect_error(
        beta_test_table(toy$table, toy$metadata, "Group", distance = "U", q = c(0, 1)),
        "single order"
    )
})

test_that("Faith PD counts each branch once up to the root", {
    small <- ape::read.tree(text = "((a:1,b:2):3,c:4);")
    tab <- data.frame(s1 = c(1, 1, 0), s2 = c(1, 0, 1), row.names = c("a", "b", "c"))
    expect_equal(.mbm_classic_weighted(tab, "phylogenetic", small, NULL)[[1]], c(6, 8))
})

test_that("Rao Q fits the functional Hill number of order 2", {
    rao <- .mbm_classic_weighted(counts, "functional", NULL, traits)[[1]]
    tau <- max(hilldiv3::traits2dist(traits))
    q2 <- .mbm_hill_weighted(counts, "functional", NULL, traits, tau = tau)$q2
    expect_equal(1 / (1 - rao / tau), q2)
})

test_that("tau defaults to the mean distance of the pooled samples", {
    p <- rowSums(counts) / sum(counts)
    d <- as.matrix(hilldiv3::traits2dist(traits))
    expect_equal(
        .mbm_hill_weighted(counts, "functional", NULL, traits),
        .mbm_hill_weighted(counts, "functional", NULL, traits, tau = sum(p %*% d %*% p))
    )
})

test_that("alpha_diversity_plot draws Faith PD and Rao Q and keeps hashed taxa", {
    hashed <- counts
    rownames(hashed) <- paste0(seq_len(nrow(hashed)), "abc")
    tips <- tree
    tips$tip.label <- rownames(hashed)[match(tree$tip.label, rownames(counts))]
    file <- tempfile(fileext = ".txt")
    p <- suppressMessages(alpha_diversity_plot(as.data.frame(hashed), toy$metadata,
        x_col = "Group", fill_col = "Group",
        diversity = "phylogenetic", tree = tips,
        save_table = TRUE, table_filename = file
    ))
    expect_s3_class(p, "ggplot")
    expect_true("Faith PD" %in% names(utils::read.delim(file, check.names = FALSE)))

    q <- alpha_diversity_plot(toy$table[, colnames(counts)], toy$metadata,
        x_col = "Group", fill_col = "Group",
        diversity = "functional", traits = traits
    )
    expect_s3_class(q, "ggplot")
    expect_error(
        alpha_diversity_plot(toy$table, toy$metadata, x_col = "Group", fill_col = "Group", diversity = "x"),
        "Invalid diversity"
    )
})
