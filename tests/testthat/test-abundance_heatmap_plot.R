test_that("abundance_heatmap_plot returns a grob, with or without drawing", {
  table <- read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
                      row.names = 1, check.names = FALSE)
  metadata <- read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
                         check.names = FALSE)
  colnames(metadata)[1] <- "SampleID"

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  for (d in c(TRUE, FALSE)) {
    g <- suppressWarnings(abundance_heatmap_plot(
      table = table, metadata = metadata, condition1 = "Location",
      top_n = 10, show_column_names = FALSE, draw = d))
    expect_s3_class(g, "gTree")
    expect_s3_class(g, "mbm_heatmap")
  }
  # printing the returned object draws the heatmap
  expect_no_error(print(g))
})

test_that("composite SILVA names become '<last genus> group'", {
  labs <- MicroBioMeta:::.mbm_shorten_labels(c(
    "ASV8_Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
    "ASV3_Escherichia-Shigella",
    "ASV20_Gitt-GS-136"))
  expect_equal(labs, c("ASV8_Rhizobium group", "ASV3_Escherichia-Shigella",
                       "ASV20_Gitt-GS-136"))
})

test_that(".mbm_composite_genus keeps prefixes and suffixes and leaves NA", {
  expect_equal(
    MicroBioMeta:::.mbm_composite_genus(c(
      "other Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
      "Burkholderia-Caballeronia-Paraburkholderia",
      "Escherichia-Shigella", "Bacillota-D", NA)),
    c("other Rhizobium group", "Paraburkholderia group",
      "Escherichia-Shigella", "Bacillota-D", NA))
})
