#' PERMANOVA or betadisper results table
#'
#' Runs a PERMANOVA (vegan::adonis2()) or a test of homogeneity of dispersions
#' (vegan::betadisper()) on the distances between samples and returns the
#' results as a table figure.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata.
#' Its first column must hold the sample IDs (the column names of `table`).
#' @param formula_str Character. Right-hand side of the model formula, with
#'   columns of \code{metadata} (e.g. \code{"Type_of_soil*Treatment"}).
#' @param distance Character. Distance metric: \code{"compositional"} (CLR with
#'   ALDEx2), \code{"aitchison"}, \code{"robust.aitchison"} or any other method
#'   of \code{vegan::vegdist()} (e.g. \code{"bray"}, \code{"jaccard"}). Default
#'   \code{"euclidean"}. Case-insensitive. Ignored when \code{table} is already a
#'   distance.
#' @param test Statistical test to run: one of \code{"permanova"} (default) or
#'   \code{"betadisper"}. Case-insensitive.
#' @param permutations Integer. Number of permutations for the test. Default
#'   \code{999}.
#' @param strata_var Character. Name of a column in \code{metadata} within
#'   which the permutations are restricted (strata). Optional; \code{NULL}
#'   (default) for none.
#' @param digits Integer. Decimal places of the numeric columns (except the
#'   p-value). Default \code{3}.
#' @param save_table Logical. If \code{TRUE}, saves the results table as a
#'   tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"beta_test_results.txt"}.
#' @param mc_samples Number of ALDEx2 Monte Carlo instances used when
#'   \code{distance = "compositional"}. With \code{1} (default) the clr values
#'   of one random instance are used: fast, but the result changes
#'   between runs (use \code{set.seed()}). With more, the clr values are
#'   averaged across instances, which gives an almost identical result in every
#'   run; \code{128} (ALDEx2's default) is suggested for final analyses, and
#'   takes longer. Ignored for other distances.
#'
#' @details The first column of \code{metadata} must hold the sample IDs;
#'   metadata rows are matched to the samples by ID, so their order doesn't
#'   matter. A precomputed distance is used as-is (\code{distance} is ignored).
#'   The p-values come from permutations; call \code{set.seed()} before the
#'   function to make the result reproducible.
#'
#' @return A ggplot object (\code{ggpubr::ggtexttable()}) with the results
#'   table; significant p-values are in bold. Use \code{save_table = TRUE}
#'   to get the results as a tab-delimited file.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' # Example using a data frame
#' beta_test_table(
#'     table = table,
#'     metadata = metadata,
#'     formula_str = "Location*Treatment",
#'     distance = "bray",
#'     test = "permanova",
#'     permutations = 999,
#'     strata_var = "Plot"
#' )
#'
#' # Example using a distance matrix
#' dist_matrix <- vegan::vegdist(
#'     t(table[, setdiff(colnames(table), "taxonomy")]),
#'     method = "bray"
#' )
#' beta_test_table(
#'     table       = dist_matrix,
#'     metadata    = metadata,
#'     formula_str = "Location",
#'     test        = "betadisper"
#' )
#'
#' # Compositional PERMANOVA (CLR via ALDEx2, then Euclidean)
#' beta_test_table(
#'     table = table,
#'     metadata = metadata,
#'     formula_str = "Location",
#'     distance = "compositional",
#'     test = "permanova",
#'     permutations = 999
#' )
beta_test_table <- function(table,
                            metadata,
                            formula_str,
                            distance = "euclidean",
                            test = c("permanova", "betadisper"),
                            permutations = 999,
                            mc_samples = 1,
                            strata_var = NULL,
                            digits = 3,
                            save_table = FALSE,
                            table_filename = "beta_test_results.txt") {
    test <- match.arg(tolower(test), c("permanova", "betadisper"))
    distance <- tolower(distance)
    is_square_dist <- is.matrix(table) && nrow(table) == ncol(table) &&
        isSymmetric(unname(table)) && all(diag(table) == 0)
    precomputed <- inherits(table, "dist") || is_square_dist
    raw_input <- is.data.frame(table)

    if (precomputed) {
        if (distance == "compositional") {
            stop(
                "distance = 'compositional' requires a raw abundance table ",
                "(data frame with a taxonomy column), not a precomputed distance matrix."
            )
        }
        dist_matrix <- stats::as.dist(table)
        sample_ids <- labels(dist_matrix)
    } else if (is.data.frame(table)) {
        tax_cols <- grep("taxonomy|taxon|Taxonomy|Taxa", names(table))
        if (length(tax_cols) > 0) {
            table <- table[, -tax_cols[1], drop = FALSE]
        }

        num_cols <- vapply(table, is.numeric, logical(1))
        if (!all(num_cols)) {}
        table <- as.matrix(table[, num_cols, drop = FALSE])
    } else if (!is.matrix(table)) {
        stop("'table' must be a matrix, dataframe, or dist object")
    }

    meta_ids <- trimws(as.character(metadata[[1]]))
    if (!precomputed) {
        n_cols_in_meta <- sum(colnames(table) %in% meta_ids)
        n_rows_in_meta <- sum(rownames(table) %in% meta_ids)
        if (n_cols_in_meta == 0 && n_rows_in_meta == 0) {
            stop("None of the sample names in table are in the first column of metadata.")
        }
        if (n_cols_in_meta >= n_rows_in_meta) table <- t(table)
        sample_ids <- rownames(table)
    }

    if (is.null(sample_ids)) {
        stop(
            "table has no sample names; they are needed to match it with ",
            "the sample IDs in the first column of metadata."
        )
    }
    metadata <- .mbm_align_metadata(sample_ids, metadata)
    keep <- metadata[[1]]
    if (precomputed) {
        if (length(keep) < length(sample_ids)) {
            dist_matrix <- stats::as.dist(as.matrix(dist_matrix)[keep, keep])
        }
    } else {
        table <- table[keep, , drop = FALSE]
    }

    if (test == "permanova") {
        vars <- all.vars(as.formula(paste("~", formula_str)))
        vars_in_metadata <- vars %in% colnames(metadata)
        if (!all(vars_in_metadata)) {
            stop("Variables ", paste(vars[!vars_in_metadata], collapse = ", "), " are not in metadata")
        }
    }

    strata <- NULL
    if (!is.null(strata_var)) {
        if (!strata_var %in% colnames(metadata)) {
            stop("Strata variable ", strata_var, " is not in metadata")
        }
        strata <- metadata[[strata_var]]
    }

    if (precomputed) {} else if (distance == "compositional") {
        if (!raw_input) {
            stop(
                "distance = 'compositional' requires a raw abundance table ",
                "(data frame with a taxonomy column), not a precomputed distance matrix."
            )
        }
        clr_samples <- .mbm_aldex_clr(t(table), mc_samples)
        dist_matrix <- stats::dist(clr_samples, method = "euclidean")
    } else if (distance %in% c("aitchison", "robust.aitchison")) {
        dist_matrix <- vegan::vegdist(table, method = distance, pseudocount = 0.5)
    } else {
        dist_matrix <- vegan::vegdist(table, method = distance)
    }

    if (test == "permanova") {
        resultado <- vegan::adonis2(as.formula(paste("dist_matrix ~", formula_str)),
            data = metadata,
            permutations = permutations,
            strata = strata,
            by = "terms"
        )
        tabla <- as.data.frame(resultado)
        tabla$Term <- rownames(tabla)
        rownames(tabla) <- NULL
    }

    if (test == "betadisper") {
        var_group <- all.vars(as.formula(paste("~", formula_str)))[1]
        disp <- vegan::betadisper(dist_matrix, metadata[[var_group]])
        perm <- vegan::permutest(disp, permutations = permutations)

        tabla <- as.data.frame(perm$tab)
        tabla$Term <- ifelse(rownames(tabla) == "Groups", var_group, rownames(tabla))
        rownames(tabla) <- NULL
    }

    col_p_name <- grep("Pr", names(tabla), ignore.case = TRUE, value = TRUE)
    tabla <- tabla %>%
        dplyr::mutate(across(where(is.numeric) & !dplyr::any_of(col_p_name), ~ round(., digits))) %>%
        dplyr::mutate(across(dplyr::any_of(col_p_name), ~ .mbm_format_pval(.))) %>%
        dplyr::mutate(across(everything(), as.character)) %>%
        dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "-", .)))

    if (save_table) {
        utils::write.table(tabla,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    tabla_display <- tabla
    names(tabla_display)[names(tabla_display) == "R2"] <- "R\u00b2"

    tab <- ggpubr::ggtexttable(tabla_display,
        rows = NULL,
        theme = ggpubr::ttheme(
            colnames = ggpubr::colnames_style(
                fill = "gray",
                color = "black",
                face = "bold",
                size = 12,
                fontface = "plain",
                fontfamily = "serif"
            ),
            tbody.style = ggpubr::tbody_style(
                fill = "white",
                color = "black",
                size = 12,
                fontfamily = "serif"
            )
        )
    )

    tab <- tab %>%
        ggpubr::tab_add_hline(at.row = 1, row.side = "top", linewidth = 4) %>%
        ggpubr::tab_add_hline(at.row = 2, row.side = "top", linewidth = 4)

    col_p <- grep("Pr", names(tabla), ignore.case = TRUE)
    if (length(col_p) > 0) {
        p_text <- tabla[[col_p]]
        is_below_threshold <- startsWith(p_text, "<")
        p_values <- suppressWarnings(as.numeric(p_text))
        filas_signif <- which(is_below_threshold | (!is.na(p_values) & p_values < 0.05))
        if (length(filas_signif) > 0) {
            for (fila in filas_signif) {
                tab <- ggpubr::table_cell_font(tab, row = fila + 1, column = col_p, face = "bold")
            }
        }
    }

    tab
}
