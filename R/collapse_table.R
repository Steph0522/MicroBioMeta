#' Collapse table
#' Collapse an OTU/ASV table by taxonomic level
#'
#' Collapses an abundance table by a specified taxonomic rank (e.g. genus,
#' family, phylum), summing counts across features that share the same taxonomy.
#' Features with lower taxonomic resolution than the selected level are retained
#' unchanged. Optionally converts counts to relative abundance and exports the
#' collapsed table to a tab-delimited file.
#'
#' @param table A data frame containing an OTU/ASV abundance table. Must include
#'   a column with OTUID, one with taxonomy and sample columns with numeric counts.
#' @param metadata A data frame containing sample metadata. The first column is
#'   used to define the order of samples in the output table.
#' @param level Character string specifying the taxonomic level to collapse to.
#'   ("kingdom","phylum","class","order","family","genus" or"species".
#' @param rel_abun Logical. If TRUE, converts counts to relative abundance
#'   (%) per sample.
#' @param save_table Logical. If \code{TRUE}, saves the collapsed table to
#'   disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"collapsed_table.txt"}.
#' @return A list with two elements: \code{collapsed_table} (wide format,
#'   one row per collapsed taxon) and \code{long_format} (the same data
#'   pivoted to one row per taxon/sample combination).
#' @export
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' collapse_table(
#'   table      = table,
#'   metadata   = metadata,
#'   level      = "genus",
#'   rel_abun   = FALSE,
#'   save_table = FALSE
#' )

collapse_table <- function(table,
                          metadata,
                          level = "genus",
                          rel_abun = FALSE,
                          save_table = FALSE,
                          table_filename = "collapsed_table.txt") {
  # --- Order columns ---
  table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
  ordered_samples <- intersect(metadata[, 1], colnames(table)[-1])
  table <- table[, c("taxonomy", ordered_samples)]

  # Save rownames as a temporary column, if present
  if (!is.null(rownames(table))) {
    table <- tibble::rownames_to_column(table, var = "OTU_ID")
  } else {
    table$OTU_ID <- paste0("OTU_", seq_len(nrow(table)))
  }

  # --- Strip trailing empty taxonomy levels ---
  table$taxonomy <- gsub("(;__)+$", "", table$taxonomy)

  # --- Map the requested level to its position in the taxonomy string ---
  level_idx <- switch(
    level,
    kingdom = 1,
    phylum = 2,
    class = 3,
    order = 4,
    family = 5,
    genus = 6,
    species = 7
  )

  # --- Compute the resolution depth of each row ---
  get_depth <- function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    sum(grepl("__", levels))
  }
  table$depth <- vapply(table$taxonomy, get_depth, integer(1))

  # --- Split rows by resolution depth ---
  lowres <-
    table[table$depth < level_idx,]    # under-resolved, kept as-is
  highres <-
    table[table$depth >= level_idx,]  # resolved enough to collapse

  # --- Collapse highres rows down to the requested level ---
  highres$taxonomy <- vapply(highres$taxonomy, function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    paste(levels[seq_len(level_idx)], collapse = ";")
  }, character(1))

  # --- Group and sum, keeping one original ID ---
  highres <- highres %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(
      OTU_ID = dplyr::first(OTU_ID),
      # keep one of the original IDs
      dplyr::across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
      .groups = "drop"
    )

  # --- Combine lowres and highres ---
  table_final <- dplyr::bind_rows(lowres[, c("OTU_ID", "taxonomy", ordered_samples)],
                                  highres[, c("OTU_ID", "taxonomy", ordered_samples)])

  # --- Drop the depth column ---
  table_final <-
    table_final[, c("OTU_ID", "taxonomy", ordered_samples)]
  table_final <- tibble::column_to_rownames(table_final, "OTU_ID")

  # --- Convert to relative abundance ---
  if (rel_abun) {
    table_final[, ordered_samples] <- sweep(table_final[, ordered_samples, drop = FALSE],
                                            2,
                                            colSums(table_final[, ordered_samples, drop = FALSE], na.rm = TRUE),
                                            FUN = "/") * 100
  }

  # --- Long format ---
  table_long <- table_final %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(ordered_samples),
      names_to = "SAMPLEID",
      values_to = ifelse(rel_abun, "RelativeAbundance", "Counts")
    )

  # --- Save table to disk if requested ---
  if (save_table) {
    utils::write.table(
      table_final,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
    message(paste("Table saved as:", table_filename))
  }

  return(list(collapsed_table = table_final,
              long_format = table_long))
}
