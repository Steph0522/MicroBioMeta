#' Beta diversity ordination plot
#'
#' Computes the beta diversity between samples with several distance metrics
#' (including the compositional CLR/Aitchison distance via ALDEx2) and plots a
#' PCA, PCoA or NMDS ordination.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata.
#' Its first column must hold the sample IDs (the column names of `table`).
#' @param distance Distance method: one of `"euclidean"`, `"bray"`,
#'   `"jaccard"`, `"sorensen"`, `"compositional"` (default; CLR/Aitchison via
#'   ALDEx2), `"aitchison"`, or `"robust.aitchison"`. Case-insensitive.
#'   Note: `ordination = "PCA"` requires `distance = "compositional"`.
#' @param mc_samples Number of ALDEx2 Monte Carlo instances used when
#'   \code{distance = "compositional"}. With \code{1} (default) the clr values
#'   of one random instance are used: fast, but the result changes
#'   between runs (use \code{set.seed()}). With more, the clr values are
#'   averaged across instances, which gives an almost identical result in every
#'   run; \code{128} (ALDEx2's default) is suggested for final analyses, and
#'   takes longer. Ignored for other distances.
#' @param ordination Ordination method: one of `"PCA"` (default), `"PCoA"`, or
#'   `"NMDS"`. Case-insensitive.
#' @param group_col Character. Name of the column in \code{metadata} used to
#'   color the points; if it is numeric, a continuous color scale is used.
#'   Optional; \code{NULL} (default) for no groups.
#' @param palette Either a palette \strong{name} or a \strong{vector of fixed
#'   colors}; which scale it produces depends on whether \code{group_col} is
#'   discrete or continuous.
#'   \itemize{
#'     \item Named, discrete \code{group_col}: one of \code{"colorb"}
#'       (default; qualitative colorblind-friendly palette), \code{"grey"},
#'       \code{"viridis"}, or \code{"brewer"} (\code{"Set2"}).
#'     \item Named, continuous \code{group_col}: \code{"viridis"} (default;
#'       \code{option = "cividis"})
#'       or \code{"gradient"}
#'     \item Vector of colors, discrete \code{group_col}: used as-is, one
#'       color per level (\code{scale_*_manual}).
#'     \item Vector of colors, continuous \code{group_col}: used as gradient
#'       stops (\code{scale_*_gradientn}).
#'   }
#' @param shape_col Character. Name of the column in \code{metadata} used for
#'   the point shapes. Optional; \code{NULL} (default) for one shape.
#' @param legend_title Character. Title of the legend. If \code{NULL}
#'   (default), the name of \code{group_col} is used.
#' @param taxonomy_db Character. Database the taxonomy strings come from:
#'   \code{"silva"} (default), \code{"gg2"} (Greengenes2, also \code{"gg"}),
#'   \code{"unite"} or \code{"Kraken2"} (also \code{"kraken"}). Case-insensitive.
#'   Only used when \code{ordination = "PCA"}.
#' @param top_n Integer. Number of taxa that contribute most to the PCA, drawn
#'   as arrows. Default \code{5}.
#' @param arrows_size Numeric. Size/length scaling factor for biplot arrows. Default \code{10}.
#' @param title Character. Plot title. \code{"auto"} (default) shows
#'   \code{"Ordination - <distance>"}; \code{NULL} shows no title; any other text
#'   is used as the title.
#' @param save_table Logical. If \code{TRUE}, saves the ordination scores and
#'   loadings (one table) as a tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"ordination_scores.txt"}.
#' @details With \code{distance = "compositional"} and \code{mc_samples = 1},
#'   the clr values come from one random Monte Carlo instance of
#'   \code{ALDEx2::aldex.clr()}; call \code{set.seed()} before the function
#'   to make the result reproducible, or use \code{mc_samples = 128}.
#'
#' @return A ggplot object.
#' @export
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' beta_ord_plot(
#'     table      = table,
#'     metadata   = metadata,
#'     distance   = "aitchison",
#'     ordination = "NMDS",
#'     group_col  = "Location",
#'     top_n      = 5
#' )
beta_ord_plot <- function(table, metadata,
                          distance = "compositional",
                          mc_samples = 1,
                          ordination = "PCA",
                          group_col = NULL,
                          palette = "colorb",
                          shape_col = NULL,
                          legend_title = NULL,
                          taxonomy_db = "silva",
                          arrows_size = 10,
                          top_n = 5,
                          title = "auto",
                          save_table = FALSE,
                          table_filename = "ordination_scores.txt") {
    requireNamespace("vegan")
    requireNamespace("ggplot2")
    requireNamespace("ggrepel")
    requireNamespace("ALDEx2")
    requireNamespace("dplyr")
    requireNamespace("stringr")
    distance <- tolower(distance)

    ordination <- switch(toupper(ordination),
        "PCA" = "PCA",
        "PCOA" = "PCoA",
        "NMDS" = "NMDS",
        stop("Invalid `ordination`: '", ordination,
            "'. Choose one of: \"PCA\", \"PCoA\", \"NMDS\" (case-insensitive).",
            call. = FALSE
        )
    )

    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")

    if (ordination == "PCA" && distance != "compositional") {
        stop("PCA is only available with 'compositional' (Aitchison). Use PCoA or NMDS instead.")
    }

    tax_col <- names(table)[ncol(table)]
    taxonomy <- table[[tax_col]]
    abund_table <- table[, -ncol(table)]

    if ("Feature.ID" %in% names(abund_table)) {
        feature_ids <- abund_table$Feature.ID
        abund_table <- abund_table[, !(names(abund_table) == "Feature.ID")]
    } else {
        feature_ids <- rownames(abund_table)
    }
    rownames(abund_table) <- feature_ids
    otu_table <- as.data.frame(lapply(abund_table, as.numeric), check.names = FALSE)
    rownames(otu_table) <- feature_ids


    metadata_ids <- trimws(as.character(metadata[[1]]))
    sample_ids <- trimws(colnames(otu_table))
    colnames(otu_table) <- sample_ids
    metadata[[1]] <- metadata_ids

    metadata <- .mbm_align_metadata(sample_ids, metadata)
    otu_table <- otu_table[, metadata[[1]], drop = FALSE]

    if (distance == "compositional") {
        otu_trans <- .mbm_aldex_clr(otu_table, mc_samples)
        dist_matrix <- dist(otu_trans, method = "euclidean")
    } else if (distance %in% c("aitchison", "robust.aitchison")) {
        dist_matrix <- vegan::vegdist(t(otu_table), method = distance, pseudocount = 0.5)
        otu_trans <- NULL
    } else if (distance == "sorensen") {
        dist_matrix <- vegan::vegdist(t(otu_table), method = "bray", binary = TRUE)
        otu_trans <- NULL
    } else {
        dist_matrix <- vegan::vegdist(t(otu_table), method = distance)
        otu_trans <- NULL
    }

    expl_var <- NULL
    res <- switch(ordination,
        "PCA" = {
            pca_input <- if (!is.null(otu_trans)) otu_trans else t(otu_table)
            pca <- prcomp(pca_input)
            list(ord = pca, expl_var = round(100 * summary(pca)$importance[2, seq_len(2)], 1))
        },
        "PCoA" = {
            pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
            eigs <- pcoa$eig
            list(ord = pcoa, expl_var = round(100 * eigs[seq_len(2)] / sum(eigs[eigs > 0]), 1))
        },
        "NMDS" = list(ord = vegan::metaMDS(dist_matrix, k = 2, trymax = 100), expl_var = NULL),
        stop("Invalid ordination method.")
    )

    ord_res <- res$ord
    expl_var <- res$expl_var


    ord_df <- switch(ordination,
        "PCA" = as.data.frame(ord_res$x),
        "PCoA" = {
            df <- as.data.frame(ord_res$points)
            colnames(df) <- c("PCoA1", "PCoA2")
            df
        },
        "NMDS" = as.data.frame(ord_res$points)
    )

    ord_df$SampleID <- rownames(ord_df)
    colnames(metadata)[1] <- "SampleID"
    merged <- dplyr::inner_join(ord_df, metadata, by = "SampleID")
    if (nrow(merged) == 0) stop("No common samples between table and metadata.")

    x_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[1], " (", expl_var[1], "%)") else names(ord_df)[1]
    y_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[2], " (", expl_var[2], "%)") else names(ord_df)[2]

    is_continuous <- !is.null(group_col) && is.numeric(merged[[group_col]])
    legend_name <- if (!is.null(legend_title)) legend_title else group_col
    colorb_default <- if (!is_continuous && !is.null(group_col) &&
        length(unique(merged[[group_col]])) == 2) {
        .mbm_colors_2group
    } else {
        .mbm_colors
    }


    is_named_palette <- is.character(palette) && length(palette) == 1

    .group_scale <- function(aesthetic) {
        fill <- aesthetic == "fill"
        if (is_continuous) {
            if (!is_named_palette) {
                if (fill) {
                    ggplot2::scale_fill_gradientn(name = legend_name, colours = palette)
                } else {
                    ggplot2::scale_color_gradientn(name = legend_name, colours = palette)
                }
            } else if (palette == "gradient") {
                if (fill) {
                    ggplot2::scale_fill_gradient(name = legend_name, low = "#0072B2", high = "#E69F00")
                } else {
                    ggplot2::scale_color_gradient(name = legend_name, low = "#0072B2", high = "#E69F00")
                }
            } else {
                if (fill) {
                    ggplot2::scale_fill_viridis_c(name = legend_name, option = "cividis")
                } else {
                    ggplot2::scale_color_viridis_c(name = legend_name, option = "cividis")
                }
            }
        } else {
            if (!is_named_palette) {
                if (fill) {
                    ggplot2::scale_fill_manual(name = legend_name, values = palette)
                } else {
                    ggplot2::scale_color_manual(name = legend_name, values = palette)
                }
            } else {
                switch(palette,
                    "grey" = if (fill) {
                        ggplot2::scale_fill_grey(name = legend_name, start = 0.9, end = 0.3)
                    } else {
                        ggplot2::scale_color_grey(name = legend_name, start = 0.9, end = 0.3)
                    },
                    "viridis" = if (fill) {
                        ggplot2::scale_fill_viridis_d(name = legend_name)
                    } else {
                        ggplot2::scale_color_viridis_d(name = legend_name)
                    },
                    "brewer" = if (fill) {
                        ggplot2::scale_fill_brewer(name = legend_name, palette = "Set2")
                    } else {
                        ggplot2::scale_color_brewer(name = legend_name, palette = "Set2")
                    },
                    if (fill) {
                        ggplot2::scale_fill_manual(name = legend_name, values = colorb_default)
                    } else {
                        ggplot2::scale_color_manual(name = legend_name, values = colorb_default)
                    }
                )
            }
        }
    }

    if (is.null(shape_col)) {
        p <- ggplot2::ggplot(merged, ggplot2::aes(
            x = .data[[names(ord_df)[1]]],
            y = .data[[names(ord_df)[2]]],
            fill = if (!is.null(group_col)) .data[[group_col]]
        )) +
            (if (is.null(group_col)) {
                ggplot2::geom_point(size = 4, shape = 21, fill = "grey60")
            } else {
                ggplot2::geom_point(size = 4, shape = 21)
            }) +
            .group_scale("fill")
    } else {
        p <- ggplot2::ggplot(merged, ggplot2::aes(
            x = .data[[names(ord_df)[1]]],
            y = .data[[names(ord_df)[2]]],
            color = if (!is.null(group_col)) .data[[group_col]],
            shape = .data[[shape_col]]
        )) +
            ggplot2::geom_point(size = 4) +
            .group_scale("color")
    }

    p <- p +
        ggplot2::geom_vline(xintercept = 0, linetype = 2) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2) +
        ggplot2::labs(
            x = x_lab,
            y = y_lab,
            title = if (identical(title, "auto")) {
                paste(ordination, "-", distance)
            } else {
                title
            },
            caption = if (ordination == "NMDS") sprintf("Stress = %.3f", ord_res$stress) else NULL
        ) +
        .mbm_theme(
            legend_position = "right",
            extra = ggplot2::theme(
                legend.box = "vertical",
                panel.grid.major = ggplot2::element_blank(),
                plot.caption = ggplot2::element_text(
                    size = 11, color = "black",
                    hjust = 1
                )
            )
        )

    rot_df_out <- NULL
    if (ordination == "PCA") {
        rot_df <- as.data.frame(ord_res$rotation)
        rot_df$Feature.ID <- rownames(rot_df)
        rot_df$mag <- sqrt(rot_df$PC1^2 + rot_df$PC2^2)
        rot_df <- rot_df[order(rot_df$mag, decreasing = TRUE), ][seq_len(top_n), ]
        rot_df$PC1 <- rot_df$PC1 * arrows_size
        rot_df$PC2 <- rot_df$PC2 * arrows_size
        rot_df$Taxon <- taxonomy[match(rot_df$Feature.ID, feature_ids)]


        extract_clean_label <- function(taxon_string, taxonomy_db = "silva") {
            if (is.na(taxon_string) || taxon_string == "" || taxon_string == "Other") {
                return("Other")
            }

            levels <- unlist(strsplit(taxon_string, ";"))
            levels <- trimws(levels)

            clean_by_db <- list(
                silva   = function(x) sub("^[a-zA-Z]__", "", x),
                gg2     = function(x) sub("^[a-zA-Z]__", "", x),
                unite   = function(x) sub("^[a-zA-Z]__", "", x),
                kraken2 = function(x) sub("^[a-zA-Z]__", "", x, )
            )

            taxonomy_db_norm <- tolower(taxonomy_db)
            cleaner <- clean_by_db[[taxonomy_db_norm]]

            if (is.null(cleaner)) {
                cleaner <- function(x) sub(".*__", "", x)
            }

            levels_clean <- vapply(levels, cleaner, FUN.VALUE = character(1))
            levels_clean <- trimws(levels_clean)

            invalid_literals <- c(
                "", " ", "NA", "na", "unclassified", "Unassigned",
                "uncultured", "uncultured_soil", "metagenome", "__"
            )

            invalid_regex <- c("bacteriap[0-9]+")

            for (i in length(levels_clean):1) {
                lvl <- levels_clean[i]

                if (!(lvl %in% invalid_literals) &&
                    !any(grepl(invalid_regex, lvl, ignore.case = TRUE))) {
                    if (taxonomy_db_norm == "kraken2" && grepl("s__", levels[i])) {
                        genus_full <- stringr::str_extract(taxon_string, "g__[^;]*")

                        if (!is.na(genus_full)) {
                            genus <- sub("g__", "", genus_full)
                            species <- lvl

                            if (genus != "" && species != "") {
                                species <- gsub("_", " ", species)
                                return(paste(genus, species))
                            }
                        }
                    }

                    return(lvl)
                }
            }

            return("Unclassified")
        }


        rot_df$label <- vapply(rot_df$Taxon, extract_clean_label, character(1), taxonomy_db = taxonomy_db)
        rot_df$label <- gsub("_", " ", rot_df$label)
        rot_df$label <- gsub(" ", "\n", rot_df$label)
        rot_df_out <- rot_df

        p <- p +
            ggplot2::geom_segment(
                data = rot_df,
                ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
                arrow = ggplot2::arrow(length = grid::unit(0.7, "cm")),
                color = "gray30",
                inherit.aes = FALSE
            ) +
            ggrepel::geom_label_repel(
                data = rot_df,
                ggplot2::aes(x = PC1, y = PC2, label = label),
                fill = "white", color = "black",
                fontface = "italic", size = 4, family = "serif",
                inherit.aes = FALSE
            )
    }

    x_limits <- range(merged[[names(ord_df)[1]]], na.rm = TRUE)
    y_limits <- range(merged[[names(ord_df)[2]]], na.rm = TRUE)
    max_range <- max(abs(x_limits), abs(y_limits))
    p <- p + ggplot2::coord_cartesian(
        xlim = c(-max_range, max_range),
        ylim = c(-max_range, max_range)
    ) +
        ggplot2::coord_fixed()

    if (save_table) {
        site_out <- merged
        colnames(site_out)[colnames(site_out) %in% names(ord_df)[seq_len(2)]] <- c("Axis1", "Axis2")
        site_out$type <- "site"
        if (!is.null(rot_df_out)) {
            loading_out <- rot_df_out
            colnames(loading_out)[colnames(loading_out) %in% c("PC1", "PC2")] <- c("Axis1", "Axis2")
            loading_out$type <- "loading"
            combined_table <- dplyr::bind_rows(site_out, loading_out)
        } else {
            combined_table <- site_out
        }
        utils::write.table(combined_table,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    return(p)
}
