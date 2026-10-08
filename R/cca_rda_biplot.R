#' CCA/RDA biplot
#'
#' Performs a Canonical Correspondence Analysis (CCA) or a Redundancy Analysis
#' (RDA) of the abundance table constrained by environmental variables, and
#' plots the samples and the environmental vectors as a biplot.
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
#' @param env_vars A character vector with the names of environmental variables
#'   (columns of \code{metadata}, or of \code{env_data} if given) to include in
#'   the analysis.
#' @param analysis Constrained ordination method: `"CCA"` (default, Canonical
#'   Correspondence Analysis) or `"RDA"` (Redundancy Analysis). Case-insensitive.
#' @param method Character. Transformation of the abundances, passed to
#'   \code{vegan::decostand()}. Default \code{"hell"} (Hellinger).
#' @param group_col Character. Name of the column in \code{metadata} that
#'   defines the groups. Used to color the samples. Optional; \code{NULL}
#'   (default) for no groups.
#' @param group_colors Optional character vector of colors, one per group,
#'   named after the groups or in their order. If \code{NULL} (default), the
#'   colorblind-friendly Okabe-Ito palette is used.
#' @param legend_title Character. Title of the legend. If \code{NULL}
#'   (default), the name of \code{group_col} is used.
#' @param scale_env Logical. If \code{TRUE} (default), the environmental
#'   variables are scaled (without centering) before the analysis.
#' @param scale_arrows Numeric. Multiplies the length of the environmental
#'   arrows; it only changes how they are drawn. Default \code{1}.
#' @param pval_threshold Numeric. P-value cutoff (after \code{p_adjust_method})
#'   to draw only the significant environmental variables. Default \code{0.05}.
#' @param p_adjust_method Character. Multiple-testing correction for the
#'   \code{vegan::envfit()} p-values of the environmental variables, any method
#'   of \code{stats::p.adjust()}. Default \code{"none"} (raw p-values), the usual
#'   choice with few variables.
#' @param show_all_env_vectors Logical. If \code{TRUE}, all the environmental
#'   variables are drawn, significant or not. Default \code{FALSE}.
#' @param title Character. Plot title. \code{"auto"} (default) shows \code{"CCA
#'   Biplot"} or \code{"RDA Biplot"}; \code{NULL} shows no title; any other text
#'   is used as the title.
#' @param save_table Logical. If \code{TRUE}, saves the sample scores and
#'   environmental vector loadings (one table, with a \code{type} column:
#'   \code{"site"} or \code{"vector"}) as a tab-delimited file. Default
#'   \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"cca_rda_scores.txt"}.
#'
#' @details The p-values of the environmental vectors come from
#'   \code{vegan::envfit()} permutations; call \code{set.seed()} before the
#'   function to make them reproducible.
#'
#' @return A ggplot object with the biplot. The axis titles show the percentage
#'   of the variance (inertia) explained by each constrained axis.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' # The environmental variables are columns of metadata, chosen with env_vars
#' cca_rda_biplot(
#'     table                = table,
#'     metadata             = metadata,
#'     env_vars             = c("pH", "TOC", "FW", "Root_FW", "DW"),
#'     analysis             = "RDA",
#'     show_all_env_vectors = TRUE,
#'     group_col            = "Location"
#' )
cca_rda_biplot <- function(table,
                           env_data = NULL,
                           env_vars,
                           method = "hell",
                           metadata,
                           group_col = NULL,
                           group_colors = NULL,
                           legend_title = NULL,
                           scale_env = TRUE,
                           pval_threshold = 0.05,
                           p_adjust_method = "none",
                           show_all_env_vectors = FALSE,
                           analysis = "CCA",
                           scale_arrows = 1,
                           title = "auto",
                           save_table = FALSE,
                           table_filename = "cca_rda_scores.txt") {
    analysis <- match.arg(toupper(analysis), c("CCA", "RDA"))

    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")

    taxonomy_col <- ncol(table)
    spp_table <- table[, -taxonomy_col, drop = FALSE]

    spp_table <- t(spp_table)

    colnames(metadata)[1] <- "SampleID"
    metadata <- as.data.frame(metadata)
    metadata <- .mbm_align_metadata(rownames(spp_table), metadata)
    spp_table <- spp_table[metadata$SampleID, , drop = FALSE]
    rownames(metadata) <- metadata$SampleID

    env_data <- .mbm_env_data(env_data, metadata, env_vars)
    missing_env <- setdiff(rownames(spp_table), rownames(env_data))
    if (length(missing_env) > 0) {
        stop(
            "Samples missing from the row names of env_data: ",
            paste(head(missing_env), collapse = ", ")
        )
    }
    env_data <- env_data[rownames(spp_table), , drop = FALSE]

    spp_hell <- vegan::decostand(spp_table, method = method)

    if (scale_env) {
        env_scaled <- scale(env_data[, env_vars, drop = FALSE], scale = TRUE, center = FALSE) %>% as.data.frame()
    } else {
        env_scaled <- env_data[, env_vars, drop = FALSE]
    }

    stopifnot(identical(rownames(spp_table), rownames(env_data)))

    if (toupper(analysis) == "RDA") {
        ord_result <- vegan::rda(spp_hell ~ ., data = env_scaled)
        axis_names <- c("RDA1", "RDA2")
    } else {
        ord_result <- vegan::cca(spp_hell ~ ., data = env_scaled)
        axis_names <- c("CCA1", "CCA2")
    }

    fit <- vegan::envfit(ord_result, env_scaled)

    if (show_all_env_vectors) {
        vars_to_plot <- rownames(vegan::scores(fit, display = "vectors"))
    } else {
        env_p <- stats::p.adjust(fit$vectors$pvals, method = p_adjust_method)
        sig_vars <- names(which(env_p < pval_threshold))
        if (length(sig_vars) == 0) {
            warning(
                "No significant environmental variables (p < ", pval_threshold,
                ", p_adjust_method = '", p_adjust_method, "')."
            )
            return(ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::ggtitle("No significant environmental variables"))
        }
        vars_to_plot <- sig_vars
    }

    vectors_scores <- vegan::scores(fit, display = "vectors")[vars_to_plot, , drop = FALSE] %>%
        as.data.frame()
    vectors_scores$Variable <- rownames(vectors_scores)

    site_scores <- vegan::scores(ord_result, display = "sites") %>% as.data.frame()
    colnames(site_scores)[seq_len(2)] <- axis_names
    site_scores$SampleID <- rownames(site_scores)

    if (save_table) {
        site_out <- site_scores
        colnames(site_out)[colnames(site_out) %in% axis_names] <- c("Axis1", "Axis2")
        site_out$type <- "site"
        vec_out <- vectors_scores
        colnames(vec_out)[colnames(vec_out) %in% axis_names] <- c("Axis1", "Axis2")
        vec_out$type <- "vector"
        combined_table <- dplyr::bind_rows(site_out, vec_out)
        utils::write.table(combined_table,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    if (!is.null(metadata) && !is.null(group_col) && group_col %in% colnames(metadata)) {
        if (!"SampleID" %in% colnames(metadata)) {
            metadata$SampleID <- rownames(metadata)
        }
        site_scores <- merge(site_scores, metadata[, c("SampleID", group_col)], by = "SampleID", all.x = TRUE)
        colnames(site_scores)[colnames(site_scores) == group_col] <- "Group"

        groups_present <- unique(site_scores$Group)
        default_colors <- if (length(groups_present) == 2) .mbm_colors_2group else .mbm_colors

        if (is.null(group_colors)) {
            color_values <- rep(default_colors, length.out = length(groups_present))
            names(color_values) <- groups_present
            group_colors <- color_values
        }

        legend_name <- ifelse(is.null(legend_title), group_col, legend_title)

        plot <- ggplot2::ggplot(site_scores, ggplot2::aes(x = .data[[axis_names[1]]], y = .data[[axis_names[2]]], fill = Group)) +
            ggplot2::geom_point(size = 4, shape = 21) +
            ggplot2::scale_fill_manual(name = legend_name, values = group_colors)
    } else {
        plot <- ggplot2::ggplot(site_scores, ggplot2::aes(x = .data[[axis_names[1]]], y = .data[[axis_names[2]]])) +
            ggplot2::geom_point(size = 4, shape = 21)
    }

    plot <- plot +
        ggplot2::geom_segment(
            data = vectors_scores,
            ggplot2::aes(
                x = 0, y = 0,
                xend = .data[[axis_names[1]]] * scale_arrows,
                yend = .data[[axis_names[2]]] * scale_arrows
            ),
            arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),
            color = "black",
            inherit.aes = FALSE
        ) +
        ggrepel::geom_label_repel(
            data = vectors_scores,
            ggplot2::aes(
                x = .data[[axis_names[1]]] * scale_arrows * 1.08,
                y = .data[[axis_names[2]]] * scale_arrows * 1.08,
                label = Variable
            ),
            color = "black",
            size = 5,
            family = "serif",
            fontface = "bold",
            fill = grDevices::adjustcolor("white", alpha.f = 0.75),
            label.size = NA,
            label.padding = grid::unit(0.1, "lines"),
            box.padding = 0.3,
            point.padding = 0,
            min.segment.length = Inf,
            max.overlaps = Inf,
            seed = 1,
            inherit.aes = FALSE
        ) +
        ggplot2::coord_fixed(ratio = 1)

    max_range <- max(abs(c(
        site_scores[[axis_names[1]]],
        vectors_scores[[axis_names[1]]] * scale_arrows,
        site_scores[[axis_names[2]]],
        vectors_scores[[axis_names[2]]] * scale_arrows
    )))
    buffer <- 1.2
    plot <- plot +
        ggplot2::scale_x_continuous(limits = c(-max_range, max_range) * buffer) +
        ggplot2::scale_y_continuous(limits = c(-max_range, max_range) * buffer) +
        ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        .mbm_theme(
            legend_position = "right",
            extra = ggplot2::theme(
                aspect.ratio = 1,
                legend.box = "vertical",
                panel.grid.major = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(
                    fill = NA, colour = "black",
                    linewidth = 0.5
                )
            )
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = legend_title))

    axis_pct <- 100 * vegan::eigenvals(ord_result)[axis_names] / ord_result$tot.chi
    axis_labels <- sprintf("%s (%.1f%%)", axis_names, axis_pct)
    auto_title <- paste(toupper(analysis), "Biplot")
    plot <- plot + ggplot2::labs(
        title = if (identical(title, "auto")) auto_title else title,
        x = axis_labels[1],
        y = axis_labels[2]
    )

    return(plot)
}
