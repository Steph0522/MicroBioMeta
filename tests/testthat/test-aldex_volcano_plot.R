test_that("aldex_volcano_plot returns a ggplot for the effect-size view", {
    toy <- make_toy_community()

    set.seed(1)
    p <- aldex_volcano_plot(
        table = toy$table,
        metadata = toy$metadata,
        group_col = "Group",
        type = "effect",
        cond = "A"
    )

    expect_s3_class(p, "ggplot")
})

test_that("aldex_volcano_plot returns a ggplot for the volcano view", {
    toy <- make_toy_community()

    set.seed(1)
    p <- aldex_volcano_plot(
        table = toy$table,
        metadata = toy$metadata,
        group_col = "Group",
        type = "volcano",
        cond = "A"
    )

    expect_s3_class(p, "ggplot")
})

test_that("aldex_volcano_plot saves the ALDEx2 result table when requested", {
    toy <- make_toy_community()

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)

    set.seed(1)
    aldex_volcano_plot(
        table = toy$table,
        metadata = toy$metadata,
        group_col = "Group",
        type = "effect",
        cond = "A",
        save_table = TRUE,
        table_filename = tmp
    )

    expect_true(file.exists(tmp))
    saved <- utils::read.delim(tmp, check.names = FALSE)
    expect_true(all(c("effect", "wi.ep", "wi.eBH") %in% names(saved)))
})

test_that("aldex_volcano_plot rejects an invalid type", {
    toy <- make_toy_community()

    expect_error(
        aldex_volcano_plot(
            table = toy$table,
            metadata = toy$metadata,
            group_col = "Group",
            type = "not-a-type"
        ),
        "`type` must be either"
    )
})

test_that("aldex_volcano_plot: in the plot, positive effect means higher in cond", {
    toy <- make_toy_community()

    run <- function(cond) {
        tmp <- tempfile(fileext = ".txt")
        on.exit(unlink(tmp))
        set.seed(1)
        p <- aldex_volcano_plot(toy$table, toy$metadata,
            group_col = "Group",
            type = "effect", cond = cond,
            save_table = TRUE, table_filename = tmp
        )
        saved <- utils::read.delim(tmp, check.names = FALSE)
        list(
            plot = p$data$effect[p$data$Feature.ID == "OTU1"],
            saved = saved$effect[saved$Feature.ID == "OTU1"]
        )
    }

    a <- run("A")
    b <- run("B")
    expect_gt(a$plot, 0)
    expect_lt(b$plot, 0)
    expect_lt(a$saved, 0)
    expect_equal(a$saved, b$saved)
})

test_that("aldex_volcano_plot: p_adjust_method picks the adjusted or raw p-value", {
    toy <- make_toy_community()

    set.seed(1)
    p_bh <- aldex_volcano_plot(toy$table, toy$metadata,
        group_col = "Group",
        type = "effect", cond = "A", p_adjust_method = "BH"
    )
    set.seed(1)
    p_raw <- aldex_volcano_plot(toy$table, toy$metadata,
        group_col = "Group",
        type = "effect", cond = "A", p_adjust_method = "none"
    )

    expect_equal(p_bh$data$.p, p_bh$data$wi.eBH)
    expect_equal(p_raw$data$.p, p_raw$data$wi.ep)
    expect_error(aldex_volcano_plot(toy$table, toy$metadata,
        group_col = "Group",
        type = "effect", p_adjust_method = "holm"
    ))
})
