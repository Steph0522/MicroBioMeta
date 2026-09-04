#' Merge feature table and taxonomy
#'
#' This function joins a feature table (counts per sample) with the corresponding taxonomy data.
#'
#' @param table A data frame or matrix with samples as columns and taxa (features) as rows.
#'              Row names must contain OTUIDs, ASVs, or species names.
#' @param taxonomy A data frame or matrix with taxonomy information.
#'                 Row names must match the identifiers in the table.
#' @param save_table Logical. If \code{TRUE}, saves the merged table to disk.
#'   Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"merged_feature_taxonomy.txt"}.
#'
#' @return A data frame with counts and taxonomy merged. If only one column in the taxonomy is present,
#'         it will be renamed to 'taxonomy'.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' full_table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' # Split the combined table into a counts-only feature table and a
#' # separate taxonomy data frame, as merge_feature_taxonomy expects them
#' feature_table <- full_table[, setdiff(colnames(full_table), "taxonomy")]
#' taxonomy_df <- full_table[, "taxonomy", drop = FALSE]
#'
#' merged <- merge_feature_taxonomy(
#'   table    = feature_table,
#'   taxonomy = taxonomy_df
#' )

merge_feature_taxonomy <- function(table, taxonomy,
                                   save_table = FALSE,
                                   table_filename = "merged_feature_taxonomy.txt") {
  # Ensure both inputs are data frames
  table <- as.data.frame(table)
  taxonomy <- as.data.frame(taxonomy)
  
  # Extract rownames (feature IDs)
  ids_table <- rownames(table)
  ids_taxonomy <- rownames(taxonomy)
  common_species <- intersect(ids_table, ids_taxonomy)
  
  # Validations
  if (length(common_species) == 0) {
    stop("No matching IDs found between the feature table and the taxonomy.")
  }
  
  if (!all(ids_taxonomy %in% ids_table)) {
    missing_from_table <- setdiff(ids_taxonomy, ids_table)
    warning(length(missing_from_table), 
            " taxonomy IDs are not present in the table and will be excluded from the result.")
  }
  
  # Subset to common features only
  table <- table[common_species, , drop = FALSE]
  taxonomy <- taxonomy[common_species, , drop = FALSE]
  
  # Prepare for join: convert rownames to a column
  table$OTUID <- rownames(table)
  taxonomy$OTUID <- rownames(taxonomy)
  
  # Join by OTUID
  joined <- dplyr::left_join(table, taxonomy, by = "OTUID")
  
  # Rename the taxonomy column to 'taxonomy' if there's only one column
  taxonomy_cols <- setdiff(colnames(taxonomy), "OTUID")
  if (length(taxonomy_cols) == 1) {
    colnames(joined)[colnames(joined) == taxonomy_cols] <- "taxonomy"
  }
  
  # Set OTUID as rownames again
  rownames(joined) <- joined$OTUID
  joined$OTUID <- NULL

  if (save_table) {
    utils::write.table(joined, file = table_filename, sep = "\t",
                       quote = FALSE, col.names = NA)
    message(paste("Table saved as:", table_filename))
  }

  return(joined)
}
