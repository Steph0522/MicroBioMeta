#' ALDEx2 volcano or effect size plot
#'
#' Runs ALDEx2 between two groups and plots the results as a volcano plot or an
#' effect size plot.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata.
#' Its first column must hold the sample IDs (the column names of `table`).
#' @param group_col Character. Name of the column in \code{metadata} that
#'   defines the groups. It must have exactly two groups.
#' @param type Type of plot to generate: `"volcano"` (default, volcano plot)
#'   or `"effect"` (effect-size plot). Case-insensitive.
#' @param mc_samples Integer. Number of ALDEx2 Monte Carlo instances. Default
#'   \code{128} (ALDEx2's default), recommended for final analyses. Fewer
#'   instances are faster but the p-values and effect sizes are less stable (e.g.
#'   \code{16} for a quick look); call \code{set.seed()} before the function to
#'   make the result reproducible.
#' @param col_inf Character. Color of the taxa lower in \code{cond}. Default
#'   \code{"#56B4E9"} (Okabe-Ito blue).
#' @param col_sup Character. Color of the taxa higher in \code{cond}. Default
#'   \code{"#E69F00"} (Okabe-Ito orange).
#' @param threshold_lower Numeric. Lower cutoff on the x-axis (effect size or
#'   difference between groups), drawn as a dashed line. Default \code{-1.5}.
#' @param threshold_upper Numeric. Upper cutoff on the x-axis (effect size or
#'   difference between groups), drawn as a dashed line. Default \code{1.5}.
#' @param cond Name of the reference condition for the plot labels: positive
#'   values on the x-axis mean higher in \code{cond}. Default \code{NULL}, the
#'   condition of the first sample in \code{table}.
#' @param pval_threshold Numeric. P-value cutoff (after
#'   \code{p_adjust_method}), drawn as the dashed horizontal line. Default
#'   \code{0.05}.
#' @param p_adjust_method Character. \code{"BH"} (default) or \code{"none"}:
#'   use the Benjamini-Hochberg adjusted p-values of ALDEx2 (\code{wi.eBH}) or
#'   the raw ones (\code{wi.ep}) for the y-axis and \code{pval_threshold}. ALDEx2
#'   only computes the BH correction.
#' @param show_labels Logical. If \code{TRUE} (default), shows the "Higher
#'   in"/"Lower in" \code{cond} labels (only for \code{type = "effect"}).
#' @param taxa Optional data frame with columns \code{Feature.ID} and
#'   \code{Taxon}, used to label the taxa. If \code{NULL} (default), the
#'   \code{taxonomy} column of \code{table} is used.
#' @param x_axis_title,y_axis_title Titles for the x- and y-axis. Default
#'   \code{NULL}: \code{"Effect size"} (\code{"effect"}) or log2 fold change
#'   (\code{"volcano"}) on x, and the (adjusted) p-value on y.
#' @param label_size Numeric. Font size of the taxon labels drawn on
#'   significant points in the "volcano" plot. Default \code{3.5}.
#' @param filter_uncultured Logical. If \code{TRUE}, taxa whose name contains
#'   "uncultured"/"unculture" are still plotted as points but not labelled,
#'   keeping the volcano plot's text annotations to named taxa. Default
#'   \code{FALSE}.
#' @param save_table Logical. If \code{TRUE}, saves the ALDEx2 results table as
#'   a tab-delimited file. Default \code{FALSE}. \code{effect} and
#'   \code{diff.btw} are saved as ALDEx2 returns them (second group in
#'   alphabetical order minus the first); in the plot, positive means higher in
#'   \code{cond}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"aldex_pval_effect.txt"}.
#'
#' @return A ggplot object with the selected plot.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' # group_col must have exactly two groups; Location has two
#' # (Rhizosphere and Roots) in the bundled example data
#' aldex_volcano_plot(
#'     table = table,
#'     metadata = metadata,
#'     group_col = "Location",
#'     mc_samples = 16,
#'     type = "effect",
#'     col_inf = "#56B4E9",
#'     col_sup = "#E69F00",
#'     threshold_lower = -0.5,
#'     threshold_upper = 0.5,
#'     cond = "Rhizosphere",
#'     show_labels = TRUE
#' )
aldex_volcano_plot <- function(table,
                               metadata,
                               group_col,
                               type = "volcano",
                               mc_samples = 128,
                               col_inf = "#56B4E9",
                               col_sup = "#E69F00",
                               threshold_lower = -1.5,
                               threshold_upper = 1.5,
                               cond = NULL,
                               pval_threshold = 0.05,
                               p_adjust_method = "BH",
                               show_labels = TRUE,
                               taxa = NULL,
                               x_axis_title = NULL,
                               y_axis_title = NULL,
                               label_size = 3.5,
                               filter_uncultured = FALSE,
                               save_table = FALSE,
                               table_filename = "aldex_pval_effect.txt") {
    if (is.logical(p_adjust_method)) p_adjust_method <- if (isTRUE(p_adjust_method)) "BH" else "none"

    type <- tolower(type)
    if (!type %in% c("effect", "volcano")) {
        stop("`type` must be either \"effect\" or \"volcano\".", call. = FALSE)
    }

    if (!requireNamespace("ALDEx2", quietly = TRUE)) {
        stop("Package 'ALDEx2' needed for this function to work. Please install it.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' needed for this function to work. Please install it.")
    }

    if (!group_col %in% colnames(metadata)) {
        stop(
            "The column ", group_col,
            " does not exist in the object 'metadata'. Check the name is written correctly."
        )
    }

    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")

    if (is.null(taxa)) {
        last_col <- ncol(table)
        taxa_colname <- colnames(table)[last_col]
        taxa <- data.frame(
            Feature.ID = rownames(table),
            Taxon = table[[taxa_colname]],
            stringsAsFactors = FALSE
        )
        table <- table[, -last_col, drop = FALSE]
    }

    metadata <- .mbm_align_metadata(colnames(table), metadata)
    table <- table[, metadata[[1]], drop = FALSE]

    conditions <- as.character(metadata[[group_col]])
    groups <- unique(conditions)
    if (is.null(cond)) cond <- groups[1]
    if (!cond %in% groups) {
        stop("`cond` must be one of: ", paste(groups, collapse = ", "), call. = FALSE)
    }
    other_cond <- setdiff(groups, cond)[1]

    aldex_clr <- ALDEx2::aldex(table, conditions, mc.samples = mc_samples, denom = "all")

    processed_data <- aldex_clr %>%
        tibble::rownames_to_column(var = "Feature.ID") %>%
        dplyr::left_join(taxa, by = "Feature.ID") %>%
        dplyr::mutate(
            taxa = dplyr::case_when(
                stringr::str_detect(Taxon, "g__") ~ stringr::str_extract(Taxon, "(?<=g__)[^_;]+"),
                stringr::str_detect(Taxon, "f__") ~ stringr::str_extract(Taxon, "(?<=f__)[^_;]+"),
                stringr::str_detect(Taxon, "c__") ~ stringr::str_extract(Taxon, "(?<=c__)[^_;]+"),
                stringr::str_detect(Taxon, "o__") ~ stringr::str_extract(Taxon, "(?<=o__)[^_;]+"),
                TRUE ~ Feature.ID
            ),
            taxa = .mbm_composite_genus(taxa)
        )

    if (save_table) {
        utils::write.table(
            processed_data,
            file = table_filename,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    if (cond == sort(groups)[1]) {
        processed_data$diff.btw <- -processed_data$diff.btw
        processed_data$effect <- -processed_data$effect
    }

    p_adjust_method <- match.arg(p_adjust_method, c("BH", "none"))
    adjusted <- p_adjust_method == "BH"
    processed_data$.p <- if (adjusted) processed_data$wi.eBH else processed_data$wi.ep
    y_lab <- if (adjusted) {
        expression("-Log"[10] ~ "adjusted p-value (BH)")
    } else {
        expression("-Log"[10] ~ "p-value")
    }

    if (type == "effect") {
        plot_data <- processed_data %>%
            dplyr::mutate(
                grupo = dplyr::case_when(
                    effect <= threshold_lower ~ paste("Lower in", cond),
                    effect >= threshold_upper ~ paste("Higher in", cond),
                    TRUE ~ "Not significant"
                ),
                log_pvalue = -log10(.p + min(.p[.p > 0]) / 10)
            )

        top_taxa <- plot_data %>%
            dplyr::filter(grupo != "Not significant") %>%
            dplyr::group_by(grupo) %>%
            dplyr::slice_max(order_by = abs(effect), n = 3)

        lim_x <- max(abs(plot_data$effect), na.rm = TRUE)

        color_values <- c(col_sup, col_inf, "gray")
        names(color_values) <- c(paste("Higher in", cond), paste("Lower in", cond), "Not significant")

        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect, y = log_pvalue, color = grupo)) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_color_manual(values = color_values) +
            ggplot2::geom_vline(
                xintercept = c(threshold_lower, threshold_upper),
                linetype = 2,
                color = "black"
            ) +
            ggplot2::geom_hline(
                yintercept = -log10(pval_threshold),
                linetype = 2,
                color = "black"
            ) +
            ggplot2::labs(
                x = if (is.null(x_axis_title)) "Effect size" else x_axis_title,
                y = if (is.null(y_axis_title)) y_lab else y_axis_title,
                color = NULL
            ) +
            .mbm_theme(legend_position = "none") +
            ggplot2::scale_x_continuous(limits = c(-lim_x, lim_x)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30)))

        if (nrow(top_taxa) > 0) {
            plot_taxa <- if (filter_uncultured) {
                top_taxa[!grepl("uncultured|unculture", top_taxa$taxa, ignore.case = TRUE), ]
            } else {
                top_taxa
            }
            if (nrow(plot_taxa) > 0) {
                p <- p +
                    ggrepel::geom_text_repel(
                        data = plot_taxa,
                        ggplot2::aes(label = taxa),
                        color = "black",
                        family = "serif",
                        size = label_size,
                        fontface = "italic",
                        box.padding = 0.5,
                        min.segment.length = 0.3,
                        segment.color = "grey50",
                        max.overlaps = Inf,
                        seed = 1
                    )
            }
        }

        if (show_labels) {
            p <- p +
                ggplot2::annotate(
                    "text",
                    x = -Inf,
                    y = Inf,
                    label = paste0("Lower in\n", cond),
                    color = col_inf,
                    family = "serif",
                    fontface = "bold",
                    size = 5,
                    hjust = 0,
                    vjust = 1.3
                ) +
                ggplot2::annotate(
                    "text",
                    x = Inf,
                    y = Inf,
                    label = paste0("Higher in\n", cond),
                    color = col_sup,
                    family = "serif",
                    fontface = "bold",
                    size = 5,
                    hjust = 1,
                    vjust = 1.3
                )
        }
    } else if (type == "volcano") {
        plot_data <- processed_data %>%
            dplyr::mutate(
                log_pvalue = -log10(.p + min(.p[.p > 0]) / 10),
                significant = .p <= pval_threshold,
                direction = ifelse(diff.btw < 0,
                    paste("Lower in", cond),
                    paste("Higher in", cond)
                )
            )

        top_taxa <- plot_data %>%
            dplyr::filter(significant) %>%
            dplyr::group_by(direction) %>%
            dplyr::slice_max(order_by = abs(diff.btw), n = 3)

        color_values <- c(col_sup, col_inf)
        names(color_values) <- c(paste("Higher in", cond), paste("Lower in", cond))

        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = diff.btw, y = log_pvalue)) +
            ggplot2::geom_point(
                data = dplyr::filter(plot_data, !significant),
                color = "gray",
                size = 3,
                alpha = 0.7
            ) +
            ggplot2::geom_point(
                data = dplyr::filter(plot_data, significant),
                ggplot2::aes(color = direction),
                size = 3
            ) +
            ggplot2::scale_color_manual(values = color_values) +
            ggplot2::geom_vline(
                xintercept = c(threshold_lower, threshold_upper),
                color = "black",
                linetype = "dashed"
            ) +
            ggplot2::geom_hline(
                yintercept = -log10(pval_threshold),
                color = "black",
                linetype = "dashed"
            ) +
            ggplot2::labs(
                x = if (is.null(x_axis_title)) expression("Log"[2] ~ "Fold Change") else x_axis_title,
                y = if (is.null(y_axis_title)) y_lab else y_axis_title,
                color = NULL
            ) +
            .mbm_theme(legend_position = "none") +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.1)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30)))

        if (nrow(top_taxa) > 0) {
            plot_taxa <- if (filter_uncultured) {
                top_taxa[!grepl("uncultured|unculture", top_taxa$taxa, ignore.case = TRUE), ]
            } else {
                top_taxa
            }
            if (nrow(plot_taxa) > 0) {
                p <- p +
                    ggrepel::geom_text_repel(
                        data = plot_taxa,
                        ggplot2::aes(label = taxa),
                        color = "black",
                        family = "serif",
                        size = label_size,
                        fontface = "italic",
                        box.padding = 0.5,
                        min.segment.length = 0.3,
                        segment.color = "grey50",
                        max.overlaps = Inf,
                        seed = 1
                    )
            }
        }

        p <- p +
            ggplot2::annotate(
                "text",
                x = -Inf,
                y = Inf,
                label = paste0("Lower in\n", cond),
                color = col_inf,
                family = "serif",
                fontface = "bold",
                size = 5,
                hjust = 0,
                vjust = 1.3
            ) +
            ggplot2::annotate(
                "text",
                x = Inf,
                y = Inf,
                label = paste0("Higher in\n", cond),
                color = col_sup,
                family = "serif",
                fontface = "bold",
                size = 5,
                hjust = 1,
                vjust = 1.3
            )
    }

    return(p)
}
