test_that("corr_env_abund_plot corrects correlation p-values before filtering", {
    table <- read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
        row.names = 1, check.names = FALSE
    )
    metadata <- read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
        check.names = FALSE
    )
    colnames(metadata)[1] <- "SampleID"
    env_data <- metadata
    rownames(env_data) <- env_data$SampleID

    run <- function(adj) {
        corr_env_abund_plot(
            table = table, env_data = env_data, metadata = metadata,
            env_vars = c("pH", "TOC", "FW", "DW"), method = "spearman",
            geom = "tile", hc.order = FALSE, show_labels = FALSE,
            level = "genus", taxonomy_db = "silva",
            pval_threshold = 0.05, p_adjust_method = adj
        )
    }

    p_raw <- run("none")
    expect_s3_class(p_raw, "ggplot")
    p_bh <- tryCatch(run("BH"), error = function(e) e)
    if (inherits(p_bh, "error")) {
        expect_match(conditionMessage(p_bh), "No taxa have a significant correlation")
    } else {
        expect_lte(nrow(p_bh$data), nrow(p_raw$data))
    }
})

test_that("corr_env_abund_plot saves the plotted correlations with taxon names", {
    table <- read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
        row.names = 1, check.names = FALSE
    )
    metadata <- read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
        check.names = FALSE
    )
    colnames(metadata)[1] <- "SampleID"
    env_data <- metadata
    rownames(env_data) <- env_data$SampleID

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)
    p <- suppressMessages(corr_env_abund_plot(
        table = table, env_data = env_data, metadata = metadata,
        env_vars = c("pH", "TOC"), method = "spearman", level = "phylum",
        taxonomy_db = "silva", save_table = TRUE, table_filename = tmp
    ))

    saved <- utils::read.delim(tmp, check.names = FALSE)
    expect_named(saved, c("Taxon", "Variable", "Correlation", "p_value", "p_adj"))
    expect_equal(nrow(saved), nrow(p$data))
    expect_setequal(saved$Taxon, unique(as.character(p$data$Taxon)))
    expect_true(all(saved$Variable %in% c("pH", "TOC")))
    expect_true(all(saved$p_value >= 0 & saved$p_value <= 1))
})
