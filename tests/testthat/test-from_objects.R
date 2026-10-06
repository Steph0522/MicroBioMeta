make_parts <- function() {
  counts <- matrix(c(10, 0, 5, 3, 8, 1), nrow = 3,
                   dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2")))
  tax <- data.frame(Kingdom = "Bacteria",
                    Phylum  = c("p__Firmicutes", "Proteobacteria", NA),
                    Genus   = c("Bacillus", NA, NA),
                    row.names = rownames(counts))
  samples <- data.frame(Soil = c("Rhizosphere", "Bulk soil"),
                        row.names = c("S1", "S2"))
  list(counts = counts, tax = tax, samples = samples)
}

expected_taxonomy <- c("k__Bacteria; p__Firmicutes; g__Bacillus",
                       "k__Bacteria; p__Proteobacteria",
                       "k__Bacteria")

test_that("from_tse returns MicroBioMeta's table and metadata", {
  skip_if_not_installed("SummarizedExperiment")
  x <- make_parts()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = x$counts), rowData = x$tax, colData = x$samples)

  out <- from_tse(se)
  expect_named(out, c("table", "metadata"))
  expect_equal(names(out$table), c("S1", "S2", "taxonomy"))
  expect_equal(out$table$taxonomy, expected_taxonomy)
  expect_equal(as.matrix(out$table[, 1:2]), x$counts)
  expect_equal(out$metadata$SAMPLEID, c("S1", "S2"))
  expect_equal(out$metadata$Soil, c("Rhizosphere", "Bulk soil"))
  expect_error(from_tse(se, assay_name = "relabundance"), "not found")
})

test_that("from_tse works with a TreeSummarizedExperiment", {
  skip_if_not_installed("TreeSummarizedExperiment")
  x <- make_parts()
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = x$counts), rowData = x$tax, colData = x$samples)
  expect_equal(from_tse(tse)$table$taxonomy, expected_taxonomy)
})

test_that("from_phyloseq returns the same as from_tse", {
  skip_if_not_installed("phyloseq")
  x <- make_parts()
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(t(x$counts), taxa_are_rows = FALSE),
    phyloseq::tax_table(as.matrix(x$tax)),
    phyloseq::sample_data(x$samples))

  out <- from_phyloseq(ps)
  expect_equal(out$table$taxonomy, expected_taxonomy)
  expect_equal(as.matrix(out$table[, 1:2]), x$counts)
  expect_equal(out$metadata$SAMPLEID, c("S1", "S2"))
  expect_error(from_phyloseq(x$counts), "phyloseq object")
})

test_that("converted tables work in MicroBioMeta functions", {
  skip_if_not_installed("SummarizedExperiment")
  toy <- make_toy_community()
  tax <- data.frame(
    Kingdom = "Bacteria",
    Phylum  = rep(c("Firmicutes", "Proteobacteria"), length.out = nrow(toy$table)),
    row.names = rownames(toy$table))
  counts <- as.matrix(toy$table[, setdiff(names(toy$table), "taxonomy")])
  samples <- toy$metadata
  rownames(samples) <- samples[[1]]
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), rowData = tax, colData = samples[, -1, drop = FALSE])

  mbm <- from_tse(se)
  expect_s3_class(
    abundance_bar_plot(mbm$table, mbm$metadata, level = "phylum", x_col = "Group"),
    "ggplot")
})
