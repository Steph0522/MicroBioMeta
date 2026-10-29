test_that("abundance_sankey_plot doesn't turn SILVA placeholders into fake taxa", {
    skip_if_not(rmarkdown::pandoc_available())

    tax <- c(
        "d__Eukaryota; p__Mucoromycota; c__Incertae_Sedis; o__Mortierellales; f__Mortierellaceae; g__Mortierella",
        "d__Eukaryota; p__Mucoromycota; c__Mucoromycetes; o__Mucorales; f__Mucoraceae; g__Mucor",
        "d__Eukaryota; p__Ascomycota; c__Sordariomycetes; o__uncultured; f__uncultured; g__uncultured",
        "d__Eukaryota; p__Nucleariidae_and_Fonticula_group; c__Nucleariidae; o__Nucleariida; f__Nucleariidae; g__Nuclearia"
    )
    counts <- matrix(c(50, 30, 20, 10, 40, 25, 35, 15, 45, 20, 30, 12),
        nrow = 4,
        dimnames = list(paste0("OTU", 1:4), paste0("S", 1:3))
    )
    tab <- as.data.frame(counts)
    tab$taxonomy <- tax

    out <- tempfile(fileext = ".txt")
    html <- tempfile(fileext = ".html")
    on.exit(unlink(c(out, html)), add = TRUE)
    suppressMessages(abundance_sankey_plot(
        table = tab, taxonomy_db = "silva", maxn = 10,
        taxRanks = c("P", "C", "O", "F", "G"), output_file = html,
        save_table = TRUE, table_filename = out
    ))

    saved <- utils::read.delim(out, check.names = FALSE)
    nodes <- saved$name[saved$table_type == "node"]
    expect_false(any(nodes %in% c("Sedis", "Incertae Sedis", "uncultured", "group")))
    expect_true("other Mucoromycota" %in% nodes)
    expect_true("Nucleariidae and Fonticula group" %in% nodes)
})

test_that("abundance_sankey_plot writes no file unless output_file is given", {
    tab <- data.frame(
        S1 = c(10, 5), S2 = c(4, 8),
        taxonomy = c(
            "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas",
            "d__Bacteria; p__Bacillota; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus"
        ),
        row.names = c("OTU1", "OTU2")
    )
    wd <- tempfile("sankey_wd")
    dir.create(wd)
    old <- setwd(wd)
    on.exit(
        {
            setwd(old)
            unlink(wd, recursive = TRUE)
        },
        add = TRUE
    )
    s <- abundance_sankey_plot(
        table = tab, taxonomy_db = "silva",
        taxRanks = c("P", "C", "G")
    )
    expect_s3_class(s, "htmlwidget")
    expect_length(list.files(wd), 0)
})

test_that("abundance_sankey_plot leaves no _files folder next to the HTML", {
    skip_if_not(rmarkdown::pandoc_available())
    tab <- data.frame(
        S1 = c(10, 5), S2 = c(4, 8),
        taxonomy = c(
            "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas",
            "d__Bacteria; p__Bacillota; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus"
        ),
        row.names = c("OTU1", "OTU2")
    )
    out_dir <- tempfile("sankey_out")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    html <- file.path(out_dir, "s.html")
    suppressMessages(abundance_sankey_plot(
        table = tab, taxonomy_db = "silva",
        taxRanks = c("P", "C", "G"), output_file = html
    ))
    expect_true(file.exists(html))
    expect_false(dir.exists(file.path(out_dir, "s_files")))
})

test_that("abundance_sankey_plot links follow the hierarchy and leave no orphan nodes", {
    tab <- data.frame(
        S1 = c(50, 30, 5, 2), S2 = c(40, 35, 4, 3),
        taxonomy = c(
            "d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas",
            "d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Vicinamibacterales; f__uncultured; g__uncultured",
            "d__Bacteria; p__Bacillota; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus",
            "d__Bacteria; p__Chloroflexota; c__KD4-96; o__KD4-96; f__KD4-96; g__KD4-96"
        ),
        row.names = paste0("OTU", 1:4)
    )
    out <- tempfile(fileext = ".txt")
    on.exit(unlink(out), add = TRUE)
    suppressMessages(abundance_sankey_plot(
        table = tab, taxonomy_db = "silva", maxn = 2,
        taxRanks = c("P", "C", "O", "F", "G"),
        save_table = TRUE, table_filename = out
    ))
    saved <- utils::read.delim(out, check.names = FALSE)
    nodes <- saved[saved$table_type == "node", ]
    links <- saved[saved$table_type == "link", ]
    depth <- nodes$depth
    # every link goes to the next column
    expect_true(all(depth[links$target + 1] - depth[links$source + 1] == 1))
    # every node has at least one link
    expect_true(all((seq_len(nrow(nodes)) - 1) %in% c(links$source, links$target)))
    # "other Vicinamibacterales" (family and genus) are two different nodes
    expect_equal(sum(nodes$name == "other Vicinamibacterales"), 2)
})
