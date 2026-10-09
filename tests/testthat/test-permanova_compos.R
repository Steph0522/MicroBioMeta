test_that("beta_test_table runs a compositional PERMANOVA and returns a table figure", {
    toy <- make_toy_community()

    set.seed(1)
    res <- beta_test_table(
        toy$table, toy$metadata,
        formula_str = "Group",
        distance = "compositional",
        test = "permanova",
        permutations = 99
    )

    expect_s3_class(res, "ggplot")
    expect_no_error(print(res))
})

test_that("beta_test_table saves a PERMANOVA table with the tested term, R2 and p-value", {
    toy <- make_toy_community()

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)

    set.seed(1)
    beta_test_table(
        toy$table, toy$metadata,
        formula_str = "Group",
        distance = "compositional",
        test = "permanova",
        permutations = 99,
        save_table = TRUE,
        table_filename = tmp
    )

    expect_true(file.exists(tmp))

    saved <- utils::read.delim(tmp, check.names = FALSE)
    expect_true("Term" %in% names(saved))
    expect_true("Group" %in% saved$Term)
    expect_true(any(grepl("R2", names(saved))))
    expect_true(any(grepl("Pr", names(saved))))
})

test_that("beta_test_table also accepts a plain ecological distance (bray)", {
    toy <- make_toy_community()

    res <- beta_test_table(
        toy$table, toy$metadata,
        formula_str = "Group",
        distance = "bray",
        test = "permanova",
        permutations = 99
    )

    expect_s3_class(res, "ggplot")
})

test_that("beta_test_table uses a precomputed distance as-is", {
    toy <- make_toy_community()
    d <- vegan::vegdist(t(toy$table[, toy$metadata$SampleID]), method = "bray")

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)
    set.seed(1)
    beta_test_table(d, toy$metadata,
        formula_str = "Group", permutations = 99,
        save_table = TRUE, table_filename = tmp
    )
    saved <- utils::read.delim(tmp, check.names = FALSE)

    set.seed(1)
    ref <- vegan::adonis2(d ~ Group, data = toy$metadata, permutations = 99, by = "terms")
    expect_equal(as.numeric(saved$R2[saved$Term == "Group"]), round(ref$R2[1], 3))
})

test_that("beta_test_table matches metadata to samples by ID, not by row order", {
    toy <- make_toy_community()
    shuffled <- toy$metadata[c(6, 1, 4, 2, 5, 3), ]

    run <- function(meta) {
        tmp <- tempfile(fileext = ".txt")
        on.exit(unlink(tmp))
        set.seed(1)
        beta_test_table(toy$table, meta,
            formula_str = "Group", distance = "bray",
            permutations = 99, save_table = TRUE, table_filename = tmp
        )
        utils::read.delim(tmp, check.names = FALSE)
    }

    expect_equal(run(shuffled), run(toy$metadata))
})
test_that("mc_samples averages the clr values over the Monte Carlo instances", {
    toy <- make_toy_community()
    counts <- as.matrix(toy$table[, toy$metadata$SampleID])

    set.seed(1)
    one <- .mbm_aldex_clr(counts, mc_samples = 1)
    expect_equal(dim(one), rev(dim(counts)))
    expect_no_warning({
        set.seed(1)
        .mbm_aldex_clr(counts, mc_samples = 1)
    })

    set.seed(1)
    avg <- .mbm_aldex_clr(counts, mc_samples = 16)
    set.seed(1)
    obj <- suppressWarnings(ALDEx2::aldex.clr(counts,
        mc.samples = 16, denom = "all",
        verbose = FALSE, useMC = FALSE
    ))
    manual <- t(sapply(ALDEx2::getMonteCarloInstances(obj), rowMeans))
    expect_equal(avg, manual)
    expect_equal(colnames(avg), rownames(counts))

    expect_error(.mbm_aldex_clr(counts, mc_samples = 0), "mc_samples")
})

test_that("beta_test_table and beta_ord_plot take mc_samples", {
    toy <- make_toy_community()
    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)

    set.seed(1)
    beta_test_table(toy$table, toy$metadata,
        formula_str = "Group",
        distance = "compositional", permutations = 99, mc_samples = 16,
        save_table = TRUE, table_filename = tmp
    )
    expect_true("Group" %in% utils::read.delim(tmp, check.names = FALSE)$Term)

    set.seed(1)
    expect_s3_class(
        beta_ord_plot(toy$table, toy$metadata, group_col = "Group", mc_samples = 16, top_n = 2),
        "ggplot"
    )
})
