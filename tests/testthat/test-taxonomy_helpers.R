test_that("taxonomy_db and level accept the same values in every function", {
    expect_equal(.mbm_taxonomy_db("SILVA"), "silva")
    expect_equal(.mbm_taxonomy_db("gg"), "gg2")
    expect_equal(.mbm_taxonomy_db("kraken"), "Kraken2")
    expect_error(.mbm_taxonomy_db("ncbi"), "Invalid `taxonomy_db`")
    expect_equal(.mbm_check_level("Genus"), "genus")
    expect_error(.mbm_check_level("specie"), "Invalid `level`")
})

test_that("ratios_bubble_plot names taxa at every level, like abundance_bar_plot", {
    table <- read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
        row.names = 1, check.names = FALSE
    )
    metadata <- read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
        check.names = FALSE
    )
    for (lv in c("phylum", "class", "order", "family", "genus", "species")) {
        p <- ratios_bubble_plot(table, metadata,
            group_col = "Location",
            condition_A = "Rhizosphere", condition_B = "Roots",
            level = lv, top_n = 10
        )
        expect_s3_class(p, "ggplot")
        # short names, not full taxonomy strings
        expect_false(any(grepl(";", p$data$taxonomy)), info = lv)
    }
    expect_error(
        ratios_bubble_plot(table, metadata,
            group_col = "Location",
            condition_A = "Rhizosphere", condition_B = "Roots",
            level = "specie"
        ),
        "Invalid `level`"
    )
})
