test_that("abundance_sankey_plot doesn't turn SILVA placeholders into fake taxa", {
  # saveNetwork() writes a self-contained HTML, which needs pandoc
  skip_if_not(rmarkdown::pandoc_available())

  tax <- c(
    "d__Eukaryota; p__Mucoromycota; c__Incertae_Sedis; o__Mortierellales; f__Mortierellaceae; g__Mortierella",
    "d__Eukaryota; p__Mucoromycota; c__Mucoromycetes; o__Mucorales; f__Mucoraceae; g__Mucor",
    "d__Eukaryota; p__Ascomycota; c__Sordariomycetes; o__uncultured; f__uncultured; g__uncultured",
    "d__Eukaryota; p__Nucleariidae_and_Fonticula_group; c__Nucleariidae; o__Nucleariida; f__Nucleariidae; g__Nuclearia"
  )
  counts <- matrix(c(50, 30, 20, 10, 40, 25, 35, 15, 45, 20, 30, 12), nrow = 4,
                   dimnames = list(paste0("OTU", 1:4), paste0("S", 1:3)))
  tab <- as.data.frame(counts)
  tab$taxonomy <- tax

  out <- tempfile(fileext = ".txt")
  html <- tempfile(fileext = ".html")
  on.exit(unlink(c(out, html)), add = TRUE)
  suppressMessages(abundance_sankey_plot(
    table = tab, taxonomy_db = "silva", maxn = 10,
    taxRanks = c("P", "C", "O", "F", "G"), output_file = html,
    save_table = TRUE, table_filename = out))

  saved <- utils::read.delim(out, check.names = FALSE)
  nodes <- saved$name[saved$table_type == "node"]
  expect_false(any(nodes %in% c("Sedis", "Incertae Sedis", "uncultured", "group")))
  expect_true("other Mucoromycota" %in% nodes)
  expect_true("Nucleariidae and Fonticula group" %in% nodes)
})
