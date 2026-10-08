#' Correlation between environmental variables and taxa
#'
#' Computes the relative abundance of the taxa at a taxonomic level, correlates
#' it with environmental variables, and shows the correlations as a heatmap
#' (tile) or a bubble plot (circle).
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata.
#' Its first column must hold the sample IDs (the column names of `table`).
#' The environmental variables can be columns of `metadata` (see `env_vars`).
#' @param env_data Optional data frame of environmental variables, with row
#'   names matching the sample names in \code{table}. Default \code{NULL}:
#'   the variables named in \code{env_vars} are taken from \code{metadata}.
#'   Use it only when the variables are in a separate table.
#' @param env_vars Character vector of environmental variables to include in
#'   the correlation analysis (columns of \code{metadata}, or of
#'   \code{env_data} if given). Required when \code{env_data} is \code{NULL};
#'   with a separate \code{env_data}, \code{NULL} uses all its variables.
#' @param method Character. Correlation method: \code{"spearman"} (default),
#'   \code{"pearson"} or \code{"kendall"}.
#' @param hc.order Logical. If \code{TRUE} (default), taxa and variables are
#'   reordered by hierarchical clustering.
#' @param geom Character. Type of visualization to generate: `"tile"`
#'   (heatmap, default) or `"circle"` (bubble plot). Case-insensitive.
#' @param show_labels Logical. If \code{TRUE} (default), the correlation values
#'   are shown.
#' @param col_palette Optional character vector of colors for the correlation
#'   scale. If \code{NULL} (default), \code{diverging_palette} is used.
#' @param diverging_palette Character. Colorblind-friendly diverging palette,
#'   used when \code{col_palette} is \code{NULL}: \code{"BuOr"} (blue-orange,
#'   default), \code{"BuVm"} (blue-vermillion), \code{"BuPk"} (blue-pink),
#'   \code{"GnPk"} (green-pink) or \code{"PuYl"} (purple-yellow); or
#'   \code{"viridis"}, a sequential scale without a neutral midpoint (not
#'   recommended for correlations).
#' @param invert_axes Logical. If \code{TRUE} (default), taxa and variables
#'   swap axes.
#' @param taxonomy_db Character. Database the taxonomy strings come from:
#'   \code{"silva"} (default), \code{"gg2"} (Greengenes2, also \code{"gg"}),
#'   \code{"unite"} or \code{"Kraken2"} (also \code{"kraken"}). Case-insensitive.
#' @param level Character. Taxonomic level: \code{"kingdom"}, \code{"phylum"},
#'   \code{"class"}, \code{"order"}, \code{"family"}, \code{"genus"} (default) or
#'   \code{"species"}. Case-insensitive.
#' @param pval_threshold Numeric or \code{NULL}. P-value cutoff (after
#'   \code{p_adjust_method}) to keep only taxa with a significant correlation
#'   with at least one variable. Default \code{NULL} (no filtering).
#' @param p_adjust_method Character. Multiple-testing correction for the
#'   p-values of all taxon x variable correlations, any method of
#'   \code{stats::p.adjust()}. Default \code{"BH"}; \code{"none"} uses the raw
#'   p-values.
#' @param x_label_angle Numeric. Rotation (degrees) of the x-axis labels.
#'   Default \code{45}; \code{0} for horizontal labels.
#' @param save_table Logical. If \code{TRUE}, saves the plotted correlations
#'   (one row per taxon and variable, with the raw and adjusted p-values) as a
#'   tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"corr.txt"}.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' corr_env_abund_plot(
#'     table = table,
#'     metadata = metadata,
#'     env_vars = c("pH", "TOC", "FW", "Root_FW", "DW", "Root_L", "Stem_L"),
#'     method = "pearson",
#'     geom = "tile",
#'     hc.order = FALSE,
#'     invert_axes = TRUE,
#'     show_labels = FALSE,
#'     level = "phylum",
#'     taxonomy_db = "silva"
#' )
corr_env_abund_plot <- function(table,
                                env_data = NULL,
                                metadata = NULL,
                                env_vars = NULL,
                                method = "spearman",
                                hc.order = TRUE,
                                geom = c("tile", "circle"),
                                show_labels = TRUE,
                                col_palette = NULL,
                                diverging_palette = "BuOr",
                                invert_axes = TRUE,
                                taxonomy_db = "silva",
                                level = "genus",
                                pval_threshold = NULL,
                                p_adjust_method = "BH",
                                x_label_angle = 45,
                                save_table = FALSE,
                                table_filename = "corr.txt") {
    geom <- match.arg(tolower(geom), c("tile", "circle"))

    taxonomy_db <- .mbm_taxonomy_db(taxonomy_db)
    level <- .mbm_check_level(level)
    rownames(table) <- NULL


    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")

    table <- table %>%
        dplyr::filter(taxonomy != "d__Bacteria;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "d__Bacteria") %>%
        dplyr::filter(taxonomy != "d__Archaea;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "d__Archaea") %>%
        dplyr::filter(taxonomy != "d__Bacteria;p__;c__;o__;f__;g__;s__") %>%
        dplyr::filter(taxonomy != "d__Archaea;p__;c__;o__;f__;g__;s__") %>%
        dplyr::filter(taxonomy != "k__Bacteria;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "k__Fungi;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "k__Fungi;p__;c__;o__;f__;g__") %>%
        dplyr::filter(taxonomy != "k__Fungi") %>%
        dplyr::filter(taxonomy != "Unassigned") %>%
        dplyr::filter(taxonomy != "d__Eukaryota")


    env_data <- .mbm_env_data(env_data, metadata, env_vars)

    if (!is.null(metadata)) {
        metadata <- as.data.frame(metadata)
        common_samples <- Reduce(intersect, list(colnames(table), rownames(env_data), metadata[[1]]))
        table <- table[, c(tax_col, match(common_samples, colnames(table))), drop = FALSE]
        env_data <- env_data[common_samples, , drop = FALSE]
        metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
        rownames(metadata) <- metadata[[1]]
        table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
        ordered_samples <- intersect(metadata[, 1], colnames(table)[-1])
        table <- table[, c("taxonomy", ordered_samples)]
    } else {
        common_samples <- intersect(colnames(table), rownames(env_data))
        table <- table[, c(tax_col, match(common_samples, colnames(table))), drop = FALSE]
        env_data <- env_data[common_samples, , drop = FALSE]
        table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
        ordered_samples <- common_samples
        table <- table[, c("taxonomy", ordered_samples)]
    }

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

    lowres <- table[table$depth < level_idx, ]
    highres <- table[table$depth >= level_idx, ]

    highres$taxonomy <- vapply(highres$taxonomy, function(tax) {
        levels <- unlist(strsplit(tax, ";"))
        paste(levels[seq_len(level_idx)], collapse = ";")
    }, character(1))

    highres <- highres %>%
        dplyr::group_by(taxonomy) %>%
        dplyr::summarise(
            OTU_ID = paste(unique(OTU_ID), collapse = ";"),
            dplyr::across(dplyr::where(is.numeric), \(x) sum(x, na.rm = TRUE)),
            .groups = "drop"
        )


    table_final <- dplyr::bind_rows(
        lowres[, c("OTU_ID", "taxonomy", ordered_samples)],
        highres[, c("OTU_ID", "taxonomy", ordered_samples)]
    )

    table_final <- table_final[, c("OTU_ID", ordered_samples, "taxonomy")]
    table_final <- tibble::column_to_rownames(table_final, "OTU_ID")

    table <- table_final

    table <- .mbm_label_taxa(table, level, taxonomy_db)

    table <- table %>%
        tibble::rownames_to_column(var = "OTU_ID")

    taxon <- table %>%
        dplyr::select(OTU_ID, taxonomy)

    table <- table %>%
        dplyr::select(-taxonomy)

    table <- table %>%
        tibble::column_to_rownames(var = "OTU_ID")


    if (is.null(col_palette)) {
        if (identical(diverging_palette, "viridis")) {
            col_palette <- scales::viridis_pal(option = "C", direction = -1)(200)
        } else {
            preset <- .mbm_div_palettes[[diverging_palette]]
            if (is.null(preset)) {
                warning(
                    "Unknown diverging_palette '", diverging_palette,
                    "'. Using 'PuYl'. Valid options: 'viridis', ",
                    paste(names(.mbm_div_palettes), collapse = ", ")
                )
                preset <- .mbm_div_palettes[["PuYl"]]
            }
            col_palette <- grDevices::colorRampPalette(preset)(200)
        }
    }
    common_samples <- base::intersect(colnames(table), rownames(env_data))
    counts <- table[, common_samples, drop = FALSE]
    env <- env_data[common_samples, , drop = FALSE]


    if (!is.null(env_vars)) {
        env_vars <- env_vars[env_vars %in% colnames(env)]
        if (length(env_vars) == 0) stop("No matching variables found in env_data")
        env <- env[, env_vars, drop = FALSE]
    }


    is_num <- vapply(env, is.numeric, logical(1))
    if (!all(is_num)) {
        warning(
            "Dropping non-numeric columns from `env_data`: ",
            paste(names(env)[!is_num], collapse = ", ")
        )
        env <- env[, is_num, drop = FALSE]
    }

    env <- env[, apply(env, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
    abund <- sweep(counts, 2, colSums(counts, na.rm = TRUE), FUN = "/") * 100
    abund <- abund[apply(abund, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

    corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")


    pval_mat <- NULL
    if (!is.null(pval_threshold) || save_table) {
        pval_mat <- .mbm_cor_pvalues(env, abund, method)
        pval_raw <- pval_mat
        pval_mat[] <- stats::p.adjust(pval_mat, method = p_adjust_method)
    }
    if (!is.null(pval_threshold)) {
        signif_taxa <- rownames(abund)[apply(pval_mat, 2, function(x) any(x < pval_threshold, na.rm = TRUE))]
        if (length(signif_taxa) == 0) {
            stop("No taxa have a significant correlation (p < ", pval_threshold,
                " after p_adjust_method = '", p_adjust_method, "').\n",
                "Try a higher pval_threshold, a coarser level, or p_adjust_method = 'none'.",
                call. = FALSE
            )
        }
        abund <- abund[signif_taxa, , drop = FALSE]
        corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
    }

    if (hc.order) {
        d_row <- stats::dist(1 - corr_mat)
        hc_row <- stats::hclust(d_row)
        row_ord <- hc_row$labels[hc_row$order]
        d_col <- stats::dist(1 - t(corr_mat))
        hc_col <- stats::hclust(d_col)
        col_ord <- hc_col$labels[hc_col$order]
        corr_mat <- corr_mat[row_ord, col_ord]
    }

    corr_df <- as.data.frame(as.table(corr_mat))
    names(corr_df) <- c("Environmental", "Group", "Correlation")
    corr_df <- dplyr::left_join(
        corr_df,
        taxon,
        by = c("Group" = "OTU_ID")
    )
    corr_df$Taxon <- corr_df$taxonomy
    corr_df$taxonomy <- NULL

    if (save_table) {
        ids <- as.character(corr_df$Group)
        vars <- as.character(corr_df$Environmental)
        out <- data.frame(
            Taxon       = corr_df$Taxon,
            Variable    = vars,
            Correlation = corr_df$Correlation,
            p_value     = pval_raw[cbind(vars, ids)],
            p_adj       = pval_mat[cbind(vars, ids)],
            check.names = FALSE
        )
        utils::write.table(out,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }
    if (invert_axes) {
        p <- ggplot2::ggplot(
            corr_df,
            ggplot2::aes(x = Environmental, y = Taxon, fill = Correlation)
        )
    } else {
        p <- ggplot2::ggplot(
            corr_df,
            ggplot2::aes(x = Taxon, y = Environmental, fill = Correlation)
        )
    }

    if (geom == "tile") {
        p <- p + ggplot2::geom_tile(color = "gray80")
        if (show_labels) {
            p <-
                p + ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                    size = 3,
                    color = "black"
                )
        }
    } else if (geom == "circle") {
        p <- p + ggplot2::geom_point(ggplot2::aes(size = abs(Correlation)),
            shape = 21,
            color = "gray"
        ) +
            ggplot2::scale_size(range = c(2, 10))
        if (show_labels) {
            p <-
                p + ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                    size = 3,
                    vjust = 0.5
                )
        }
    }

    p <- p +
        ggplot2::scale_fill_gradientn(
            colours = col_palette,
            limits = c(-1, 1),
            name = "Correlation"
        ) +
        .mbm_theme(
            legend_position = "right",
            extra = ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    angle = x_label_angle,
                    vjust = if (x_label_angle == 0) 0.5 else 1,
                    hjust = if (x_label_angle == 0) 0.5 else 1,
                    size = 12, color = "black",
                    face = if (invert_axes) "plain" else "italic"
                ),
                axis.text.y = ggplot2::element_text(
                    size = 12, color = "black",
                    face = if (invert_axes) "italic" else "plain"
                ),
                panel.grid.major = ggplot2::element_blank()
            )
        ) +
        ggplot2::coord_fixed()

    return(p)
}
.mbm_cor_pvalues <- function(env, abund, method) {
    x <- as.matrix(env)
    y <- t(abund)
    pval <- matrix(NA_real_, ncol(x), ncol(y), dimnames = list(colnames(x), colnames(y)))

    loop <- matrix(TRUE, ncol(x), ncol(y))
    if (method %in% c("pearson", "spearman")) {
        x_na <- colSums(is.na(x)) > 0
        y_na <- colSums(is.na(y)) > 0
        loop <- outer(x_na, y_na, "|")
        n <- nrow(x)
        df <- n - 2
        if (method == "spearman") {
            x_ties <- apply(x, 2, anyDuplicated) > 0
            y_ties <- apply(y, 2, anyDuplicated) > 0
            exact <- n <= 1290 & !outer(x_ties, y_ties, "|")
            loop <- loop | exact
        }
        fast <- !loop
        if (any(fast)) {
            if (method == "pearson") {
                r <- suppressWarnings(stats::cor(x, y))
                t_stat <- sqrt(df) * r / sqrt(1 - r^2)
                p_fast <- 2 * pmin(stats::pt(t_stat, df), stats::pt(t_stat, df, lower.tail = FALSE))
            } else {
                r <- suppressWarnings(stats::cor(apply(x, 2, rank), apply(y, 2, rank)))
                q <- (n^3 - n) * (1 - r) / 6
                r <- 1 - q / ((n * (n^2 - 1)) / 6)
                t_stat <- r / sqrt((1 - r^2) / df)
                p_fast <- ifelse(q > (n^3 - n) / 6,
                    stats::pt(t_stat, df),
                    stats::pt(t_stat, df, lower.tail = FALSE)
                )
                p_fast <- pmin(2 * p_fast, 1)
            }
            pval[fast] <- p_fast[fast]
        }
    }

    for (k in which(loop)) {
        i <- (k - 1) %% nrow(pval) + 1
        j <- (k - 1) %/% nrow(pval) + 1
        pval[k] <- suppressWarnings(stats::cor.test(x[, i], y[, j], method = method))$p.value
    }
    pval
}
