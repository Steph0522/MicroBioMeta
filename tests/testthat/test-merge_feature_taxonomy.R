test_that("merge_feature_taxonomy joins matching features and renames single taxonomy column", {
    table <- data.frame(
        Sample1 = c(10, 20, 30),
        Sample2 = c(5, 15, 25),
        row.names = c("OTU1", "OTU2", "OTU3")
    )
    taxonomy <- data.frame(
        Taxon = c("d__Bacteria;p__Firmicutes", "d__Bacteria;p__Bacteroidota", "d__Bacteria;p__Proteobacteria"),
        row.names = c("OTU1", "OTU2", "OTU3")
    )

    result <- merge_feature_taxonomy(table, taxonomy)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 3)
    expect_true("taxonomy" %in% colnames(result))
    expect_equal(rownames(result), c("OTU1", "OTU2", "OTU3"))
    expect_equal(result["OTU1", "Sample1"], 10)
    expect_equal(result["OTU2", "taxonomy"], "d__Bacteria;p__Bacteroidota")
})

test_that("merge_feature_taxonomy keeps only features common to both inputs", {
    table <- data.frame(
        Sample1 = c(10, 20),
        row.names = c("OTU1", "OTU2")
    )
    taxonomy <- data.frame(
        Taxon = c("d__Bacteria;p__Firmicutes", "d__Bacteria;p__Bacteroidota", "d__Bacteria;p__Proteobacteria"),
        row.names = c("OTU1", "OTU2", "OTU3")
    )

    expect_warning(
        result <- merge_feature_taxonomy(table, taxonomy),
        "not present in the table"
    )
    expect_equal(nrow(result), 2)
    expect_false("OTU3" %in% rownames(result))
})

test_that("merge_feature_taxonomy errors when there are no matching IDs at all", {
    table <- data.frame(Sample1 = 10, row.names = "OTU1")
    taxonomy <- data.frame(Taxon = "d__Bacteria", row.names = "OTU_other")

    expect_error(
        merge_feature_taxonomy(table, taxonomy),
        "No matching IDs found"
    )
})

test_that("merge_feature_taxonomy matches by ID when IDs are in the first column", {
    table <- data.frame(
        `#OTU ID` = c("OTU1", "OTU2", "OTU3"),
        Sample1 = c(10, 20, 30),
        check.names = FALSE
    )
    taxonomy <- data.frame(
        `Feature ID` = c("OTU3", "OTU1", "OTU2"),
        Taxon = c("p__Proteobacteria", "p__Firmicutes", "p__Bacteroidota"),
        check.names = FALSE
    )

    expect_message(
        expect_message(
            result <- merge_feature_taxonomy(table, taxonomy),
            "first column"
        ),
        "first column"
    )
    expect_equal(colnames(result), c("Sample1", "taxonomy"))
    expect_equal(result["OTU1", "taxonomy"], "p__Firmicutes")
    expect_equal(result["OTU3", "taxonomy"], "p__Proteobacteria")
    expect_equal(result["OTU3", "Sample1"], 30)
})

test_that("merge_feature_taxonomy errors when IDs cannot be found", {
    table <- data.frame(Sample1 = c(10, 20), Sample2 = c(1, 2))
    taxonomy <- data.frame(Taxon = c("p__A", "p__B"), row.names = c("OTU1", "OTU2"))

    expect_error(
        merge_feature_taxonomy(table, taxonomy),
        "no feature IDs as row names"
    )
})

test_that("merge_feature_taxonomy writes a table to disk when save_table = TRUE", {
    table <- data.frame(Sample1 = 10, row.names = "OTU1")
    taxonomy <- data.frame(Taxon = "d__Bacteria", row.names = "OTU1")

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)
    merge_feature_taxonomy(table, taxonomy, save_table = TRUE, table_filename = tmp)

    expect_true(file.exists(tmp))
})
