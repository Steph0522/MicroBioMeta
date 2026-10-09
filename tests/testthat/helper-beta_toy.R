make_beta_toy_community <- function() {
    incidence <- rbind(
        F1  = c(1, 1, 1, 1, 1, 1),
        F2  = c(1, 1, 1, 1, 1, 1),
        F3  = c(1, 1, 0, 1, 1, 1),
        F4  = c(1, 1, 1, 1, 0, 1),
        F5  = c(1, 1, 1, 0, 0, 0),
        F6  = c(1, 1, 0, 0, 0, 0),
        F7  = c(1, 0, 1, 0, 0, 0),
        F8  = c(0, 0, 0, 1, 1, 1),
        F9  = c(0, 0, 0, 1, 1, 0),
        F10 = c(0, 0, 0, 0, 1, 1),
        F11 = c(1, 0, 0, 1, 0, 0),
        F12 = c(0, 1, 0, 0, 1, 0)
    )
    colnames(incidence) <- paste0("S", seq_len(6))

    mult <- matrix(rep(10 + 3 * seq_len(6), each = nrow(incidence)),
        nrow = nrow(incidence)
    )
    counts <- incidence * mult

    tab <- as.data.frame(counts)
    tab$taxonomy <- paste0("d__Bacteria;p__Phylum", seq_len(nrow(incidence)))

    meta <- data.frame(
        SampleID = colnames(incidence),
        Group    = rep(c("A", "B"), each = 3),
        Batch    = rep(c("B1", "B2", "B1"), times = 2)
    )

    list(table = tab, metadata = meta)
}
