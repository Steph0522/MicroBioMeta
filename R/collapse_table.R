#' Collapse table
#'
#' Collapses an abundance table to a taxonomic level (e.g. genus, family,
#' phylum), summing the counts of the features that share the same taxonomy.
#' Features not resolved to that level are kept as they are. Optionally
#' converts the counts to relative abundance and saves the collapsed table.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' All the other columns are taken as samples, in the order of the table.
#' @param level Character. Taxonomic level: \code{"kingdom"}, \code{"phylum"},
#'   \code{"class"}, \code{"order"}, \code{"family"}, \code{"genus"} (default) or
#'   \code{"species"}. Case-insensitive.
#' @param rel_abun Logical. If \code{TRUE}, converts the counts to relative
#'   abundance (%) per sample. Default \code{FALSE}.
#' @param save_table Logical. If \code{TRUE}, saves the collapsed table as a
#'   tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"collapsed_table.txt"}.
#' @return A list with two elements: \code{collapsed_table} (wide format,
#'   one row per collapsed taxon) and \code{long_format} (the same data
#'   pivoted to one row per taxon/sample combination).
#' @export
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' collapse_table(
#'     table      = table,
#'     level      = "genus",
#'     rel_abun   = FALSE,
#'     save_table = FALSE
#' )
collapse_table <- function(table,
                           level = "genus",
                           rel_abun = FALSE,
                           save_table = FALSE,
                           table_filename = "collapsed_table.txt") {
    if (is.data.frame(level)) {
        stop("collapse_table() no longer takes `metadata`: call it as ",
            "collapse_table(table, level = ...).",
            call. = FALSE
        )
    }
    level <- .mbm_check_level(level)

    ordered_samples <- setdiff(names(table), "taxonomy")
    table <- table[, c("taxonomy", ordered_samples)]

    if (!is.null(rownames(table))) {
        table <- tibble::rownames_to_column(table, var = "OTU_ID")
    } else {
        table$OTU_ID <- paste0("OTU_", seq_len(nrow(table)))
    }

    table$taxonomy <- gsub("(;__)+$", "", table$taxonomy)

    level_idx <- switch(level,
        kingdom = 1,
        phylum = 2,
        class = 3,
        order = 4,
        family = 5,
        genus = 6,
        species = 7
    )

    get_depth <- function(tax) {
        levels <- unlist(strsplit(tax, ";"))
        sum(grepl("__", levels))
    }
    table$depth <- vapply(table$taxonomy, get_depth, integer(1))

    lowres <-
        table[table$depth < level_idx, ]
    highres <-
        table[table$depth >= level_idx, ]

    highres$taxonomy <- vapply(highres$taxonomy, function(tax) {
        levels <- unlist(strsplit(tax, ";"))
        paste(levels[seq_len(level_idx)], collapse = ";")
    }, character(1))

    highres <- highres %>%
        dplyr::group_by(taxonomy) %>%
        dplyr::summarise(
            OTU_ID = dplyr::first(OTU_ID),
            dplyr::across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
            .groups = "drop"
        )

    table_final <- dplyr::bind_rows(
        lowres[, c("OTU_ID", "taxonomy", ordered_samples)],
        highres[, c("OTU_ID", "taxonomy", ordered_samples)]
    )

    table_final <-
        table_final[, c("OTU_ID", "taxonomy", ordered_samples)]
    table_final <- tibble::column_to_rownames(table_final, "OTU_ID")

    if (rel_abun) {
        table_final[, ordered_samples] <- sweep(table_final[, ordered_samples, drop = FALSE],
            2,
            colSums(table_final[, ordered_samples, drop = FALSE], na.rm = TRUE),
            FUN = "/"
        ) * 100
    }

    table_long <- table_final %>%
        tidyr::pivot_longer(
            cols = dplyr::all_of(ordered_samples),
            names_to = "SAMPLEID",
            values_to = ifelse(rel_abun, "RelativeAbundance", "Counts")
        )

    if (save_table) {
        utils::write.table(
            table_final,
            file = table_filename,
            sep = "\t",
            quote = FALSE,
            col.names = NA
        )
        message("Table saved as: ", table_filename)
    }

    return(list(
        collapsed_table = table_final,
        long_format = table_long
    ))
}
