#' Merge feature table and taxonomy
#'
#' This function joins a feature table (counts per sample) with the corresponding taxonomy data.
#'
#' @param table A data frame with features in rows and samples in columns. Row
#'   names must hold the feature IDs (OTUs, ASVs, species...). If the table has
#'   no row names (e.g. read without \code{row.names = 1}), its first column is
#'   used as feature IDs.
#' @param taxonomy A data frame or matrix with the taxonomy. Row names (or, if
#'   absent, the first column) must hold the feature IDs of \code{table}.
#'   Features are matched by ID, so both inputs can be in any order.
#' @param save_table Logical. If \code{TRUE}, saves the merged table as a
#'   tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"merged_feature_taxonomy.txt"}.
#'
#' @return A data frame with the counts of \code{table} and the taxonomy
#'   columns. A single taxonomy column is renamed to \code{taxonomy}.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "subset_table_bacteria.txt", package = "MicroBioMeta")
#' full_table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' feature_table <- full_table[, setdiff(colnames(full_table), "taxonomy")]
#' taxonomy_df <- full_table[, "taxonomy", drop = FALSE]
#'
#' merged <- merge_feature_taxonomy(
#'     table    = feature_table,
#'     taxonomy = taxonomy_df
#' )
merge_feature_taxonomy <- function(table, taxonomy,
                                   save_table = FALSE,
                                   table_filename = "merged_feature_taxonomy.txt") {
    table <- ids_to_rownames(as.data.frame(table), "table")
    taxonomy <- ids_to_rownames(as.data.frame(taxonomy), "taxonomy")

    ids_table <- rownames(table)
    ids_taxonomy <- rownames(taxonomy)
    common_species <- intersect(ids_table, ids_taxonomy)

    if (length(common_species) == 0) {
        stop("No matching IDs found between the feature table and the taxonomy.")
    }

    if (!all(ids_taxonomy %in% ids_table)) {
        missing_from_table <- setdiff(ids_taxonomy, ids_table)
        warning(
            length(missing_from_table),
            " taxonomy IDs are not present in the table and will be excluded from the result."
        )
    }

    table <- table[common_species, , drop = FALSE]
    taxonomy <- taxonomy[common_species, , drop = FALSE]

    table$OTUID <- rownames(table)
    taxonomy$OTUID <- rownames(taxonomy)

    joined <- dplyr::left_join(table, taxonomy, by = "OTUID")

    taxonomy_cols <- setdiff(colnames(taxonomy), "OTUID")
    if (length(taxonomy_cols) == 1) {
        colnames(joined)[colnames(joined) == taxonomy_cols] <- "taxonomy"
    }

    rownames(joined) <- joined$OTUID
    joined$OTUID <- NULL

    if (save_table) {
        utils::write.table(joined,
            file = table_filename, sep = "\t",
            quote = FALSE, col.names = NA
        )
        message("Table saved as: ", table_filename)
    }

    return(joined)
}

ids_to_rownames <- function(df, arg) {
    if (.row_names_info(df) > 0) {
        return(df)
    }
    first <- df[[1]]
    if (ncol(df) > 1 && (is.character(first) || is.factor(first)) &&
        !anyNA(first) && !anyDuplicated(first)) {
        message(
            "`", arg, "` has no feature IDs as row names; using its first column (\"",
            colnames(df)[1], "\") as feature IDs."
        )
        rownames(df) <- as.character(first)
        return(df[, -1, drop = FALSE])
    }
    stop(
        "`", arg, "` has no feature IDs as row names. Read it with ",
        "`row.names = 1` or set its row names to the feature IDs."
    )
}
