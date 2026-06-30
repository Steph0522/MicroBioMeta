#' Merge feature table and taxonomy
#'
#' This function joins a feature table (counts per sample) with the corresponding taxonomy data.
#'
#' @param table A data frame or matrix with samples as columns and taxa (features) as rows.
#'              Row names must contain OTUIDs, ASVs, or species names.
#' @param taxonomy A data frame or matrix with taxonomy information.
#'                 Row names must match the identifiers in the table.
#'
#' @return A data frame with counts and taxonomy merged. If only one column in the taxonomy is present,
#'         it will be renamed to 'taxonomy'.
#' @export
#'
#' @examples
#' \dontrun{
#' merged <- merge_feature_taxonomy(
#'   table    = feature_table,
#'   taxonomy = taxonomy_df
#' )
#' }

merge_feature_taxonomy <- function(table, taxonomy) {
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
  
  return(joined)
}
