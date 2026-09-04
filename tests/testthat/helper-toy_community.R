# Shared toy community used by the compositional ordination / PERMANOVA tests.
# Two groups (A: S1-S3, B: S4-S6) of 8 features, with features 1-4 enriched in
# group A and 5-8 enriched in group B, so the compositional distance has clear,
# non-degenerate structure for the ordination and PERMANOVA to pick up.
make_toy_community <- function() {
  set.seed(1)
  n_feat <- 8
  n_samp <- 6
  counts <- matrix(
    stats::rpois(n_feat * n_samp, lambda = 50),
    nrow = n_feat,
    dimnames = list(paste0("OTU", seq_len(n_feat)), paste0("S", seq_len(n_samp)))
  )
  counts[seq_len(4), seq_len(3)] <- counts[seq_len(4), seq_len(3)] + 100
  counts[5:8, 4:6] <- counts[5:8, 4:6] + 100

  tab <- as.data.frame(counts)
  tab$taxonomy <- paste0("d__Bacteria;p__Phylum", seq_len(n_feat))

  meta <- data.frame(
    SampleID = paste0("S", seq_len(n_samp)),
    Group = rep(c("A", "B"), each = 3)
  )

  list(table = tab, metadata = meta)
}
