#' Alpha diversity plot
#'
#' This function generates a boxplot or barplot to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows.
#'              The first column must contain the OTUID, ASV, or species name.
#' @param metadata A data frame with metadata. The first column must match sample names in `table`.
#' @param type Type of plot: either "boxplot" or "barplot". Default is "boxplot".
#' @param stat Optional. A string indicating the test used for comparing means (e.g., "wilcox.test").
#' @param x_col Column in `metadata` to be used on the x-axis.
#' @param fill_col Column in `metadata` to define fill color.
#' @param facet_by Optional. A metadata column to facet alongside `q` (e.g., Treatment, Site).
#' @param facet_orientation Whether `facet_by` appears in columns ("horizontal", default) or rows ("vertical").
#' @param fill_palette Color palette to use: "colorb", "grey", "viridis", or "brewer". Default: "colorb".
#' @param custom_palette A vector of custom colors. Overrides `fill_palette` if provided.
#' @param n_cols Number of columns in facet wrap (optional).
#' @param n_rows Number of rows in facet wrap (optional).
#' @param strip_color Background color of facet strips. Default: "grey".
#' @param show_legend Logical. Show legend? Default: TRUE.
#' @param legend_title Title for the legend.
#' @param figure_title Title for the entire plot.
#' @param axis_x_title Title for the x-axis.
#' @param axis_y_title Title for the y-axis.
#' @param free_y Logical. Wether if scales in y are free or not.

#'
#' @return A ggplot object showing alpha diversity with Hill numbers.
#' @export
#' @examples
#' library(vegan)
#' data(dune)
#' data(dune.env)
#' alpha_hill_plot(
#'     table = t(dune),
#'     metadata = dune.env %>% tibble::rownames_to_column("SampleID"),
#'     x_col = "Management",
#'     fill_col = "Management",
#'     facet_by = "Use",
#'     facet_orientation = "vertical"
#' )
alpha_hill_plot <- function(
        table,
        metadata,
        type = "boxplot",
        stat = NULL,
        x_col,
        fill_col,
        facet_by = NULL,
        facet_orientation = "horizontal",
        fill_palette = "colorb",
        custom_palette = NULL,
        n_cols = NULL,
        n_rows = NULL,
        strip_color = "grey",
        show_legend = TRUE,
        figure_title = NULL,
        legend_title = NULL,
        axis_x_title = NULL,
        axis_y_title = "Effective number of ASVs",
        free_y = FALSE) {
    sample_order <- metadata[[1]]
    common_samples <- intersect(colnames(table), sample_order)
    if (length(common_samples) == 0) stop("No matching sample names between table and metadata.")

    table <- table[, common_samples, drop = FALSE]
    table <- table[, match(sample_order, colnames(table))]
    table <- data.frame(t(table))

    results <- data.frame(
        q0 = hillR::hill_taxa(comm = table, q = 0),
        q1 = hillR::hill_taxa(comm = table, q = 1),
        q2 = hillR::hill_taxa(comm = table, q = 2),
        metadata
    )

    results[[fill_col]] <- factor(results[[fill_col]], levels = unique(results[[fill_col]]))
    results[[x_col]] <- factor(results[[x_col]], levels = unique(results[[x_col]]))

    results_largo <- tidyr::pivot_longer(
        results,
        cols = tidyselect::starts_with("q"),
        names_to = "q",
        values_to = "value"
    )

    fill_scale <- if (!is.null(custom_palette)) {
        ggplot2::scale_fill_manual(values = custom_palette)
    } else {
        switch(fill_palette,
            "colorb" = ggplot2::scale_fill_manual(values = c(
                "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032",
                "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
                "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26"
            )),
            "grey" = ggplot2::scale_fill_grey(start = 0.9, end = 0.3),
            "viridis" = ggplot2::scale_fill_viridis_d(option = "plasma"),
            "brewer" = ggplot2::scale_fill_brewer(palette = "Set2"),
            ggplot2::scale_fill_grey(start = 0.9, end = 0.3)
        )
    }

    facet_scales <- if (free_y) "free_y" else "fixed"
    facet_config <- if (!is.null(facet_by)) {
        # Cuando se define orientación, usamos facet_wrap con ~ q + facet_by
        formula_facet <- stats::as.formula(paste("~", if (facet_orientation == "horizontal") {
            paste("q", facet_by, sep = " + ")
        } else {
            paste(facet_by, "q", sep = " + ")
        }))

        ggplot2::facet_wrap(
            formula_facet,
            ncol = if (!is.null(n_cols)) n_cols else NULL,
            nrow = if (!is.null(n_rows)) n_rows else NULL,
            scales = if (free_y) "free_y" else "fixed"
        )
    } else {
        # Solo hay q
        if (facet_orientation == "horizontal") {
            ggplot2::facet_wrap(~q,
                ncol = if (!is.null(n_cols)) n_cols else length(unique(results_largo$q)),
                scales = if (free_y) "free_y" else "fixed"
            )
        } else {
            ggplot2::facet_wrap(~q,
                nrow = if (!is.null(n_rows)) n_rows else length(unique(results_largo$q)),
                scales = if (free_y) "free_y" else "fixed"
            )
        }
    }




    capa_geom <- if (type == "boxplot") {
        ggplot2::geom_boxplot(
            width = 0.4,
            position = ggplot2::position_dodge(width = 0.9)
        )
    } else if (type == "barplot") {
        ggplot2::geom_bar(
            stat = "summary",
            fun = "mean",
            position = ggplot2::position_dodge(width = 0.9),
            width = 0.7,
            color = "black"
        )
    } else {
        stop("Invalid plot type. Use 'boxplot' or 'barplot'.")
    }

    capa_error <- if (type == "barplot") {
        ggplot2::stat_summary(
            fun.data = "mean_sdl",
            fun.args = list(mult = 1),
            geom = "errorbar",
            width = 0.2,
            position = ggplot2::position_dodge(width = 0.9)
        )
    } else if (type == "boxplot") {
        ggplot2::stat_summary(
            fun.data = "mean_sdl",
            fun.args = list(mult = 1),
            geom = "crossbar",
            width = 0.5,
            position = ggplot2::position_dodge(width = 0.9),
            color = "black"
        )
    } else {
        NULL
    }

    # Calcular label.y ajustado por q
    # Si hay facet_by, incluimos esa variable en el agrupamiento
    p <- ggplot2::ggplot(
        results_largo,
        ggplot2::aes(x = .data[[x_col]], y = value, fill = .data[[fill_col]])
    ) +
        capa_geom +
        capa_error +
        fill_scale +
        facet_config +
        ggplot2::labs(
            title = figure_title,
            x = axis_x_title,
            y = axis_y_title,
            fill = legend_title
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.spacing = grid::unit(1, "lines"),
            aspect.ratio = 0.4,
            strip.text = ggplot2::element_text(face = "bold", color = "black", size = 15),
            strip.background = ggplot2::element_rect(fill = strip_color),
            axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.text.x = ggplot2::element_text(size = 10),
            legend.title = ggplot2::element_text(size = 12, face = "bold"),
            legend.text = ggplot2::element_text(size = 11),
            legend.position = if (show_legend) "bottom" else "none"
        )
    if (!is.null(stat)) {
        split_vars <- if (!is.null(facet_by)) c("q", facet_by) else "q"
        p_vals_layers <- results_largo %>%
            dplyr::group_split(across(all_of(split_vars))) %>%
            purrr::map(~ {
                y_val <- max(.x$value, na.rm = TRUE) * 0.98
                ggpubr::stat_compare_means(
                    data = .x,
                    method = stat,
                    label = "p.format",
                    size = 3,
                    hide.ns = TRUE,
                    mapping = ggplot2::aes(label = ..p.format..),
                    label.y = y_val
                )
            })


        p <- p + p_vals_layers
    }

    return(p)
}
