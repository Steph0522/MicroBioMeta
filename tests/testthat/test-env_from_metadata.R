load_bacteria <- function() {
    table <- read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
        row.names = 1, check.names = FALSE
    )
    metadata <- read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
        check.names = FALSE
    )
    colnames(metadata)[1] <- "SampleID"
    list(table = table, metadata = metadata)
}

test_that("environmental variables taken from metadata give the same result as a separate env_data", {
    d <- load_bacteria()
    env_data <- d$metadata
    rownames(env_data) <- env_data$SampleID
    vars <- c("pH", "TOC", "FW", "DW")

    p_meta <- suppressMessages(corr_env_abund_plot(
        table = d$table, metadata = d$metadata, env_vars = vars,
        level = "phylum", hc.order = FALSE
    ))
    p_sep <- suppressMessages(corr_env_abund_plot(
        table = d$table, env_data = env_data, metadata = d$metadata, env_vars = vars,
        level = "phylum", hc.order = FALSE
    ))
    expect_equal(p_meta$data, p_sep$data)

    set.seed(1)
    b_meta <- suppressMessages(cca_rda_biplot(
        table = d$table, metadata = d$metadata, env_vars = vars,
        analysis = "RDA", show_all_env_vectors = TRUE
    ))
    set.seed(1)
    b_sep <- suppressMessages(cca_rda_biplot(
        table = d$table, env_data = env_data, metadata = d$metadata, env_vars = vars,
        analysis = "RDA", show_all_env_vectors = TRUE
    ))
    expect_s3_class(b_meta, "ggplot")
    expect_equal(b_meta$data, b_sep$data)
})

test_that("taking the variables from metadata needs env_vars that exist", {
    d <- load_bacteria()
    expect_error(
        corr_env_abund_plot(table = d$table, metadata = d$metadata),
        "`env_vars` is required"
    )
    expect_error(
        corr_env_abund_plot(table = d$table, env_vars = "pH"),
        "Give `metadata`"
    )
    expect_error(
        cca_rda_biplot(
            table = d$table, metadata = d$metadata,
            env_vars = c("pH", "not_a_var")
        ),
        "not found in metadata: not_a_var"
    )
})
