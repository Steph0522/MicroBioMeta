#' Box plot of pairwise feature turnover (Hill numbers)
#'
#' Computes the pairwise beta diversity between samples with Hill numbers
#' (q = 0, 1, 2; same values as \code{hillR::hill_taxa_parti_pairwise()}) and
#' plots the proportion of feature turnover (\eqn{\beta - 1}, from 0 = identical
#' to 1 = no shared features) as box plots, one facet row per Hill order.
#' q = 0 weighs all features equally (presence/absence), q = 1 weighs them by
#' their abundance and q = 2 gives more weight to dominant features.
#'
#' Each box is one group pair from \code{condition1_col}, so a between-group
#' comparison (e.g. \code{"Rhizosphere_vs_Roots"}) can be contrasted with a
#' within-group baseline (e.g. \code{"Rhizosphere_vs_Rhizosphere"}): if the
#' between-group turnover is higher, the groups host distinct communities
#' beyond the variability among replicates of the same group. An optional
#' \code{condition2_col} keeps only pairs from the same level (e.g. the same
#' treatment) and facets by it, so the comparison is not confounded by it.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata.
#' Its first column must hold the sample IDs (the column names of `table`).
#' @param comparison_condition1 Character vector with the pairs of groups of
#'   \code{condition1_col} to show (e.g. \code{"Rhizosphere_vs_Roots"},
#'   \code{"Rhizosphere_vs_Rhizosphere"}; the order of the two groups does not
#'   matter). Each pair is one box, labeled as written.
#' @param comparison_condition2 Optional character vector with the pairs of
#'   groups of \code{condition2_col} to keep (e.g. \code{"TC_vs_TC"}); the plot
#'   is faceted by them. Requires \code{condition2_col}. If \code{NULL}
#'   (default), no pairs are filtered by \code{condition2_col}.
#' @param condition1_col Character. Name of the column in \code{metadata} whose
#'   groups are compared (the pairs of \code{comparison_condition1}).
#' @param condition2_col Optional. Metadata column used to build
#'   \code{comparison_condition2}'s group pairs, and to facet the plot.
#'   Required when \code{comparison_condition2} is set. If given without
#'   \code{comparison_condition2}, the plot is faceted by it keeping only
#'   pairs of samples with the same value (e.g. \code{"TC_vs_TC"},
#'   \code{"TD_vs_TD"}).
#' @param facet_colors Optional character vector of background colors for the
#'   facet strips. If \code{NULL} (default), a neutral \code{"grey85"} is used.
#' @param group_colors Optional character vector of colors, one per comparison
#'   on the x-axis, named after them or in their order. If \code{NULL} (default),
#'   the colorblind-friendly Okabe-Ito palette is used.
#' @param x_axis_title Character. Title of the x-axis. Default
#'   \code{"Section"}.
#' @param y_axis_title Character. Title of the y-axis. Default
#'   \code{"Proportion of feature turnover"}.
#' @param show_x_labels Logical. If \code{TRUE}, the x-axis labels are shown.
#'   Default \code{FALSE}, since the groups are already named in the legend.
#' @param x_label_angle Numeric. Rotation (degrees) of the x-axis labels (when
#'   \code{show_x_labels = TRUE}). Default \code{0} (horizontal); use \code{45}
#'   or \code{90} if they overlap.
#' @param strip_text_bold Logical. If \code{TRUE}, the facet strip labels are
#'   bold. Default \code{FALSE}.
#' @param strip_text_color Color of the top (x) facet strip labels, which sit on
#'   the \code{facet_colors} backgrounds. Default \code{"white"}, unless
#'   \code{facet_colors} is left at its own light \code{"grey85"} default,
#'   in which case this defaults to \code{"black"} instead so it stays legible.
#' @param aspect_ratio Numeric. Aspect ratio (height/width) of each panel.
#'   Default \code{NULL} (automatic).
#' @param stat Character or \code{NULL}. Statistical test to compare the boxes
#'   within each panel. \code{"wilcox.test"} or \code{"t.test"} compare every
#'   pair of boxes, each with its own bracket and p-value
#'   (\code{ggpubr::stat_pwc()}); \code{"kruskal.test"} or \code{"anova"} give
#'   one global p-value per panel (\code{ggpubr::stat_compare_means()}). Default
#'   \code{NULL} (no test shown).
#' @param p_adjust_method Character. Multiple-testing correction for the
#'   pairwise tests (\code{stat = "wilcox.test"} or \code{"t.test"}), applied
#'   within each panel, any method of \code{stats::p.adjust()}. Default
#'   \code{"holm"}; \code{"none"} uses the raw p-values.
#' @param save_table Logical. If \code{TRUE}, saves the turnover table as a
#'   tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"betadiv_turnover.txt"}.
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
#' colnames(metadata)[1] <- "OTUID"
#'
#' # comparison_condition1/2 match group pairs regardless of order, so listing
#' # "Rhizosphere_vs_Roots" also catches pairs the self-join recorded the
#' # other way around ("Roots_vs_Rhizosphere") - no need to list both.
#' beta_turnover_plot(
#'     table = table,
#'     metadata = metadata,
#'     comparison_condition1 = c("Rhizosphere_vs_Roots", "Rhizosphere_vs_Rhizosphere"),
#'     comparison_condition2 = c("Control_vs_Control", "Moderate_drought_vs_Moderate_drought"),
#'     condition1_col = "Location",
#'     condition2_col = "Treatment",
#'     facet_colors = c("#5D478B", "#8B668B"),
#'     group_colors = c(
#'         "Rhizosphere_vs_Roots" = "#56B4E9",
#'         "Rhizosphere_vs_Rhizosphere" = "#E69F00"
#'     )
#' )
beta_turnover_plot <- function(table,
                               metadata,
                               comparison_condition1,
                               comparison_condition2 = NULL,
                               condition1_col,
                               condition2_col = NULL,
                               facet_colors = NULL,
                               group_colors = NULL,
                               x_axis_title = "Section",
                               y_axis_title = "Proportion of feature turnover",
                               show_x_labels = FALSE,
                               x_label_angle = 0,
                               strip_text_bold = FALSE,
                               strip_text_color = "white",
                               aspect_ratio = NULL,
                               stat = NULL,
                               p_adjust_method = "holm",
                               save_table = FALSE,
                               table_filename = "betadiv_turnover.txt") {
    colnames(metadata)[1] <- "OTUID"

    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) == 1) table <- table[, -tax_col, drop = FALSE]

    otu_filter <- table %>%
        dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric))) != 0)

    if (nrow(otu_filter) == 0) stop("OTU table is empty after filtering rows with sum=0.")

    otu_filter_t <- as.data.frame(t(otu_filter))

    if (!condition1_col %in% names(metadata)) {
        stop("condition1_col '", condition1_col, "' is not a column of metadata.")
    }
    if (!is.null(condition2_col) && is.null(comparison_condition2)) {
        if (!condition2_col %in% names(metadata)) {
            stop("condition2_col '", condition2_col, "' is not a column of metadata.")
        }
        col2 <- metadata[[condition2_col]]
        lv2 <- if (is.factor(col2)) levels(droplevels(col2)) else sort(unique(as.character(col2)))
        comparison_condition2 <- paste0(lv2, "_vs_", lv2)
        facet2_labels <- lv2
    } else {
        facet2_labels <- NULL
    }
    if (!is.null(comparison_condition2)) {
        if (is.null(condition2_col)) {
            stop("condition2_col is required when comparison_condition2 is set.")
        }
        if (!condition2_col %in% names(metadata)) {
            stop("condition2_col '", condition2_col, "' is not a column of metadata.")
        }
    }

    groups1_needed <- unique(unlist(strsplit(comparison_condition1, "_vs_")))
    keep_samples <- metadata$OTUID[metadata[[condition1_col]] %in% groups1_needed]

    if (!is.null(comparison_condition2)) {
        groups2_needed <- unique(unlist(strsplit(comparison_condition2, "_vs_")))
        keep_samples <- intersect(keep_samples, metadata$OTUID[metadata[[condition2_col]] %in% groups2_needed])
    }

    otu_filter_t <- otu_filter_t[rownames(otu_filter_t) %in% keep_samples, , drop = FALSE]
    if (nrow(otu_filter_t) == 0) {
        stop("No samples left after restricting to the groups named in comparison_condition1/2.")
    }

    beta_q_list <- list()
    for (q in c(0, 1, 2)) {
        beta_res <- tryCatch(
            {
                .mbm_hill_pairwise(otu_filter_t, q) %>%
                    dplyr::mutate(Recambio = TD_beta - 1, q = q)
            },
            error = function(e) {
                warning("Could not compute hill_taxa_parti_pairwise with q = ", q, ": ", e$message)
                NULL
            }
        )
        if (!is.null(beta_res)) beta_q_list[[as.character(q)]] <- beta_res
    }

    beta_total <- dplyr::bind_rows(beta_q_list)

    if (nrow(beta_total) == 0) stop("Could not compute beta partitions (empty table).")

    beta_formato <- beta_total %>%
        dplyr::inner_join(metadata, by = c("site1" = "OTUID")) %>%
        dplyr::inner_join(metadata, by = c("site2" = "OTUID"))

    if (nrow(beta_formato) == 0) stop("Could not join beta_total with metadata (empty table).")

    c1_x <- as.character(beta_formato[[paste0(condition1_col, ".x")]])
    c1_y <- as.character(beta_formato[[paste0(condition1_col, ".y")]])
    beta_formato$compar_condition1 <- paste0(pmin(c1_x, c1_y), "_vs_", pmax(c1_x, c1_y))

    has_condition2 <- !is.null(comparison_condition2)
    if (has_condition2) {
        c2_x <- as.character(beta_formato[[paste0(condition2_col, ".x")]])
        c2_y <- as.character(beta_formato[[paste0(condition2_col, ".y")]])
        beta_formato$compar_condition2 <- paste0(pmin(c2_x, c2_y), "_vs_", pmax(c2_x, c2_y))
    }

    pair_norm <- .mbm_normalize_pair(comparison_condition1)
    label_by_norm <- stats::setNames(comparison_condition1, pair_norm)
    beta_formato$.compar_label <- factor(
        unname(label_by_norm[beta_formato$compar_condition1]),
        levels = unique(comparison_condition1)
    )

    if (has_condition2) {
        pair_norm2 <- .mbm_normalize_pair(comparison_condition2)
        label_by_norm2 <- stats::setNames(comparison_condition2, pair_norm2)
        beta_formato$.facet2_col <- factor(
            unname(label_by_norm2[beta_formato$compar_condition2]),
            levels = unique(comparison_condition2)
        )
        if (!is.null(facet2_labels)) levels(beta_formato$.facet2_col) <- facet2_labels
    }
    has_facet_by <- has_condition2

    beta_final <- beta_formato %>%
        dplyr::filter(compar_condition1 %in% pair_norm)

    if (has_condition2) {
        beta_final <- beta_final %>%
            dplyr::filter(compar_condition2 %in% .mbm_normalize_pair(comparison_condition2))
    }

    beta_final <- beta_final %>%
        dplyr::mutate(orden = dplyr::case_when(
            q == 0 ~ "q0", q == 1 ~ "q1", q == 2 ~ "q2", TRUE ~ as.character(q)
        ))

    if (nrow(beta_final) == 0) stop("After filtering comparisons, the table is empty.")

    if (is.null(group_colors)) {
        n_groups <- length(unique(beta_final$.compar_label))
        group_colors <- if (n_groups == 2) .mbm_colors_2group else rep_len(.mbm_colors, n_groups)
        names(group_colors) <- unique(beta_final$.compar_label)
    }
    .mbm_check_color_names(group_colors, as.character(unique(beta_final$.compar_label)))
    if (has_facet_by && is.null(facet_colors)) {
        if (missing(strip_text_color)) strip_text_color <- "black"
        n_facets <- length(unique(beta_final$.facet2_col))
        facet_colors <- rep("grey85", n_facets)
    }

    if (save_table) {
        utils::write.table(beta_final,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    q_labeller <- .mbm_q_labeller(strip_text_bold)

    x_text <- if (show_x_labels) .mbm_x_text(x_label_angle) else ggplot2::element_blank()
    x_ticks <- if (show_x_labels) ggplot2::element_line(colour = "black") else ggplot2::element_blank()

    facet_spec <- if (has_facet_by) {
        ggh4x::facet_grid2(orden ~ .facet2_col,
            scales = "free_x",
            labeller = ggplot2::labeller(orden = q_labeller),
            strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = facet_colors))
        )
    } else {
        ggh4x::facet_grid2(orden ~ .,
            scales = "free_x",
            labeller = ggplot2::labeller(orden = q_labeller)
        )
    }

    figura <- beta_final %>%
        ggpubr::ggboxplot(x = ".compar_label", y = "Recambio", fill = ".compar_label") +
        ggplot2::ylab(y_axis_title) +
        ggplot2::scale_fill_manual(values = group_colors) +
        ggplot2::labs(fill = "Comparison") +
        facet_spec +
        .mbm_theme(
            legend_position = "right",
            extra = ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                axis.text.x = x_text,
                axis.ticks.x = x_ticks,
                axis.title.y = ggplot2::element_text(
                    size = 14, color = "black",
                    margin = ggplot2::margin(t = 0, r = 0.5, b = 0, l = 0, "cm")
                ),
                strip.background.y = ggplot2::element_rect(fill = "grey", color = "black"),
                strip.text.x = .mbm_strip_text(strip_text_bold, colour = strip_text_color),
                strip.text.y = .mbm_strip_text(strip_text_bold),
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                panel.spacing.y = grid::unit(0.8, "lines")
            )
        ) +
        ggplot2::xlab(x_axis_title)

    if (!is.null(aspect_ratio)) {
        figura <- figura + ggplot2::theme(aspect.ratio = aspect_ratio)
    }

    top_expand <- if (is.null(stat)) 0.05 else 0.15
    figura <- figura +
        ggplot2::scale_y_continuous(
            breaks = seq(0, 1, by = 0.25),
            expand = ggplot2::expansion(mult = c(0.08, top_expand))
        )

    if (!is.null(stat)) {
        figura <- figura + .mbm_stat_layer(stat, p_adjust_method)
    }

    figura
}
