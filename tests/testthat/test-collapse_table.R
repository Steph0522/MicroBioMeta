test_that("collapse_table sums counts of features sharing taxonomy at the requested level", {
    table <- data.frame(
        taxonomy = c(
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__casei",
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__plantarum",
            "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__fragilis"
        ),
        Sample1 = c(10, 5, 20),
        Sample2 = c(2, 3, 8),
        row.names = c("OTU1", "OTU2", "OTU3")
    )

    result <- collapse_table(table, level = "genus")

    expect_type(result, "list")
    expect_named(result, c("collapsed_table", "long_format"))

    collapsed <- result$collapsed_table
    # The two Lactobacillus features (OTU1, OTU2) should be summed into one row
    expect_equal(nrow(collapsed), 2)
    lacto_row <- collapsed[grepl("Lactobacillus$", collapsed$taxonomy), ]
    expect_equal(lacto_row$Sample1, 15)
    expect_equal(lacto_row$Sample2, 5)
})

test_that("collapse_table converts to relative abundance when rel_abun = TRUE", {
    table <- data.frame(
        taxonomy = c(
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__casei",
            "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__fragilis"
        ),
        Sample1 = c(30, 70),
        row.names = c("OTU1", "OTU2")
    )

    result <- collapse_table(table, level = "genus", rel_abun = TRUE)

    expect_equal(sum(result$collapsed_table$Sample1), 100)
})

test_that("collapse_table writes a table to disk when save_table = TRUE", {
    table <- data.frame(
        taxonomy = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__casei",
        Sample1 = 10,
        row.names = "OTU1"
    )

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)
    collapse_table(table, level = "genus", save_table = TRUE, table_filename = tmp)

    expect_true(file.exists(tmp))
})

test_that("collapse_table keeps every sample column, in table order", {
    table <- data.frame(
        taxonomy = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__casei",
        S2 = 1, S1 = 2, S3 = 3,
        row.names = "OTU1"
    )
    result <- collapse_table(table, level = "genus")
    expect_named(result$collapsed_table, c("taxonomy", "S2", "S1", "S3"))
    expect_error(
        collapse_table(table, data.frame(SampleID = "S1")),
        "no longer takes `metadata`"
    )
})
