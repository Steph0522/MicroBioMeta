test_that("taxonomy_db and level accept the same values in every function", {
    expect_equal(.mbm_taxonomy_db("SILVA"), "silva")
    expect_equal(.mbm_taxonomy_db("gg"), "gg2")
    expect_equal(.mbm_taxonomy_db("kraken"), "Kraken2")
    expect_error(.mbm_taxonomy_db("ncbi"), "Invalid `taxonomy_db`")
    expect_equal(.mbm_check_level("Genus"), "genus")
    expect_error(.mbm_check_level("specie"), "Invalid `level`")
})

test_that("ratios_bubble_plot names taxa at every level, like abundance_bar_plot", {
    table <- read.delim(system.file("extdata", "subset_table_bacteria.txt", package = "MicroBioMeta"),
        row.names = 1, check.names = FALSE
    )
    metadata <- read.delim(system.file("extdata", "subset_metadata_bacteria.txt", package = "MicroBioMeta"),
        check.names = FALSE
    )
    for (lv in c("phylum", "class", "order", "family", "genus", "species")) {
        p <- ratios_bubble_plot(table, metadata,
            group_col = "Location",
            condition_A = "Rhizosphere", condition_B = "Roots",
            level = lv, top_n = 10
        )
        expect_s3_class(p, "ggplot")
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

test_that("phylum names are updated only at the phylum rank, with one message", {
    .mbm_env$phyla_shown <- NULL
    tax <- c(
        "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria",
        "d__Bacteria;p__Firmicutes;c__Bacilli",
        "D_0__Bacteria;D_1__Cyanobacteria;D_2__Cyanobacteriia",
        "k__Bacteria; p__Deinococcus-Thermus",
        "d__Bacteria; p__Acidobacteriota; c__Actinobacteria",
        "k__Fungi; p__Ascomycota",
        NA
    )
    expect_message(out <- .mbm_update_phyla(tax), "Proteobacteria -> Pseudomonadota")
    expect_equal(out, c(
        "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria",
        "d__Bacteria;p__Bacillota;c__Bacilli",
        "D_0__Bacteria;D_1__Cyanobacteriota;D_2__Cyanobacteriia",
        "k__Bacteria; p__Deinococcota",
        "d__Bacteria; p__Acidobacteriota; c__Actinobacteria",
        "k__Fungi; p__Ascomycota",
        NA
    ))
    expect_silent(.mbm_update_phyla(tax))
    expect_identical(.mbm_update_phyla(tax[6]), tax[6])
})
