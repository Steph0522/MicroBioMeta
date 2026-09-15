# Shared toy community for beta-diversity partition tests
# (beta_dissimilarity_plot, beta_turnover_plot, beta_partition_ord_plot).
#
# Unlike make_toy_community() (whose features are present in every sample,
# varying only in abundance - fine for compositional/CLR-based distances,
# but degenerate for presence/absence partitioning since every pairwise
# Jaccard/Sorensen distance would come out identical), this fixture has a
# hand-built incidence pattern so beta diversity has real structure to
# partition:
#   - F1-F4: a shared "core", present in every sample.
#   - F5-F7: present only in group A (S1-S3).
#   - F8-F10: present only in group B (S4-S6).
#   - F11-F12: present in exactly one sample per group, for a little
#     within-group texture so distances aren't all tied at zero.
# Two categorical metadata columns (Group, Batch) are included, since
# beta_turnover_plot needs a second condition to facet/filter same-batch
# pairs on.
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

  # Turn incidence into counts: present cells get a modest, sample-varying
  # count so abundance-based metrics (e.g. Hill numbers at q = 1, 2) have
  # something to work with too, not just presence/absence.
  mult <- matrix(rep(10 + 3 * seq_len(6), each = nrow(incidence)),
                 nrow = nrow(incidence))
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
