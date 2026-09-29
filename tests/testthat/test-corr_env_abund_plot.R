test_that("corr_env_abund_plot corrects correlation p-values before filtering", {
  table <- read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
                      row.names = 1, check.names = FALSE)
  metadata <- read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
                         check.names = FALSE)
  colnames(metadata)[1] <- "SampleID"
  env_data <- metadata
  rownames(env_data) <- env_data$SampleID

  run <- function(adj) {
    corr_env_abund_plot(
      table = table, env_data = env_data, metadata = metadata,
      env_vars = c("pH", "TOC", "FW", "DW"), method = "spearman",
      geom = "tile", hc.order = FALSE, show_labels = FALSE,
      level = "genus", taxonomy_db = "silva",
      pval_threshold = 0.05, p_adjust_method = adj
    )
  }

  p_raw <- run("none")
  expect_s3_class(p_raw, "ggplot")

  # FDR correction can only keep the same taxa or fewer; if none survive,
  # the function says so clearly instead of failing inside cor().
  p_bh <- tryCatch(run("BH"), error = function(e) e)
  if (inherits(p_bh, "error")) {
    expect_match(conditionMessage(p_bh), "No taxa have a significant correlation")
  } else {
    expect_lte(nrow(p_bh$data), nrow(p_raw$data))
  }
})
