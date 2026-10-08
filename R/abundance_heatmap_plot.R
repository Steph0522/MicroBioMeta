#' Heatmap of relative abundance
#'
#' Creates a heatmap (with ComplexHeatmap) of the relative abundance of the
#' most abundant features (ASVs, OTUs or taxa).
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. Its first column
#'   must hold the sample IDs (the column names of `table`). Optional: if
#'   `NULL` (default), the heatmap has no sample annotations.
#' @param condition1 Character. Name of the column in \code{metadata} for the
#'   first column annotation. Optional; \code{NULL} (default) for none.
#' @param condition2 Character. Name of the column in \code{metadata} for the
#'   second column annotation. Optional; \code{NULL} (default) for none.
#' @param condition3 Character. Name of the column in \code{metadata} for the
#'   third column annotation. Optional; \code{NULL} (default) for none.
#' @param colors_condition1 Optional character vector of colors for
#'   \code{condition1}. If named, colors are matched to the values by name (e.g.
#'   \code{c(Roots = "#009E73", Rhizosphere = "#56B4E9")}); otherwise they are
#'   assigned in alphabetical order of the values. If \code{NULL} (default), the
#'   colorblind-friendly Okabe-Ito palette is used.
#' @param colors_condition2 Optional character vector of colors for
#'   \code{condition2} (same rules as \code{colors_condition1}). If \code{NULL}
#'   (default), colors of the colorblind-friendly "Safe" palette are used, which
#'   contrast with those of \code{condition1}.
#' @param colors_condition3 Optional character vector of colors for
#'   \code{condition3} (same rules as \code{colors_condition1}). If \code{NULL}
#'   (default), a third set of colorblind-friendly colors is used.
#' @param name_legend_condition1 Character. Title of the legend of
#'   \code{condition1}. If \code{NULL} (default), the name of \code{condition1}
#'   is used.
#' @param name_legend_condition2 Character. Title of the legend of
#'   \code{condition2}. If \code{NULL} (default), the name of \code{condition2}
#'   is used.
#' @param name_legend_condition3 Character. Title of the legend of
#'   \code{condition3}. If \code{NULL} (default), the name of \code{condition3}
#'   is used.
#' @param top_n Integer. Number of most abundant features to show. Default
#'   \code{15}.
#' @param exclude_unclassified Logical. If \code{TRUE} (default), taxa with no
#'   classification at any level ("Unclassified") are dropped before selecting
#'   the \code{top_n} most abundant features, so all the rows are identified
#'   taxa. Use \code{FALSE} to keep them.
#' @param cluster Logical. If \code{TRUE} (default), the rows are clustered; if
#'   \code{FALSE}, they are ordered by abundance.
#' @param show_column_names Logical. If \code{TRUE} (default), the sample names
#'   are shown.
#' @param cell_size Numeric. Side, in millimeters, of each (square) cell. If
#'   \code{NULL} (default), it is set automatically.
#' @param annotation_height Numeric. Height, in millimeters, of each column
#'   annotation bar (\code{condition1}, \code{condition2}, \code{condition3}). If
#'   \code{NULL} (default), the cell height is used.
#' @param save_table Logical. If \code{TRUE}, saves the abundance table shown
#'   in the heatmap as a tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"abundance_heatmap_table.txt"}.
#' @param feature_prefix Character. Prefix of the row labels, before the row
#'   number (e.g. \code{feature_prefix = "ASV"} labels the rows \code{"ASV1"},
#'   \code{"ASV2"}...). Default \code{""} (rows labeled \code{"1"},
#'   \code{"2"}...).
#' @param max_label_length Integer. Taxon names longer than this many
#'   characters are cut with an ellipsis in the row labels. \code{NULL} never
#'   cuts them. Default \code{35}.
#' @param composite_names Logical. If \code{TRUE} (default), SILVA composite
#'   names of three or more genera are labeled with the last genus plus
#'   "group" (e.g. "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
#'   becomes "Rhizobium group"); two-genus names (e.g.
#'   "Escherichia-Shigella") are kept whole. The full names are kept in the
#'   saved table.
#' @param draw Logical. If \code{TRUE} (default), the heatmap is drawn on the
#'   current device. Use \code{FALSE} to only build the returned grob without
#'   drawing it, e.g. to combine it with other plots; then, unless
#'   \code{cell_size} is given, cells stretch to fill the panel they are
#'   placed in instead of having a fixed size.
#' @return Invisibly, a \code{gTree} (grid grob) with the heatmap of the
#'   \code{top_n} most abundant features. Printing it (e.g. typing its name)
#'   draws the heatmap; it can also be combined with other plots (e.g.
#'   \code{cowplot::plot_grid()}, \code{patchwork::wrap_elements()}).
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' abundance_heatmap_plot(
#'     table                  = table,
#'     metadata               = metadata,
#'     condition1             = "Location",
#'     condition2             = "Treatment",
#'     condition3             = "Plot",
#'     top_n                  = 50,
#'     cluster                = TRUE,
#'     show_column_names      = FALSE,
#'     name_legend_condition1 = "Location",
#'     name_legend_condition2 = "Treatment",
#'     name_legend_condition3 = "Plot"
#' )
#' \donttest{
#' heat <- abundance_heatmap_plot(
#'     table = table,
#'     metadata = metadata,
#'     condition1 = "Location",
#'     condition2 = "Treatment",
#'     top_n = 20,
#'     show_column_names = FALSE,
#'     colors_condition1 = c(Rhizosphere = "#56B4E9", Roots = "#009E73"),
#'     colors_condition2 = c(
#'         Control = "#CC6677",
#'         Moderate_drought = "#332288",
#'         Severe_drought = "#999933"
#'     ),
#'     name_legend_condition1 = "Compartment",
#'     feature_prefix = "ASV",
#'     max_label_length = 30,
#'     draw = FALSE
#' )
#' heat # printing the returned object draws the heatmap
#' }
abundance_heatmap_plot <- function(table,
                                   metadata = NULL,
                                   condition1 = NULL,
                                   condition2 = NULL,
                                   condition3 = NULL,
                                   colors_condition1 = NULL,
                                   colors_condition2 = NULL,
                                   colors_condition3 = NULL,
                                   name_legend_condition1 = NULL,
                                   name_legend_condition2 = NULL,
                                   name_legend_condition3 = NULL,
                                   top_n = 15,
                                   exclude_unclassified = TRUE,
                                   cluster = TRUE,
                                   show_column_names = TRUE,
                                   cell_size = NULL,
                                   annotation_height = NULL,
                                   save_table = FALSE,
                                   table_filename = "abundance_heatmap_table.txt",
                                   feature_prefix = "",
                                   max_label_length = 35,
                                   composite_names = TRUE,
                                   draw = TRUE) {
    .mbm_check_table(table)
    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")
    if (is.null(metadata)) {
        if (!is.null(condition1) || !is.null(condition2) || !is.null(condition3)) {
            stop("`condition1`/`condition2`/`condition3` need `metadata`.", call. = FALSE)
        }
        metadata <- data.frame(SAMPLEID = names(table)[-tax_col])
    }

    taxonomy_col <- ncol(table)
    table_counts <- table[, -taxonomy_col, drop = FALSE]

    OTU_ID <- colnames(metadata)[1]

    table_tax <- table %>%
        tibble::rownames_to_column("OTUID") %>%
        dplyr::select(OTUID, taxonomy)

    phy.ra.complete <- t(t(table_counts) / colSums(table_counts) * 100) %>%
        as.data.frame()

    tax_strings <- unique(table_tax$taxonomy)
    tax_labels <- dplyr::case_when(
        grepl("g__[^;]*", tax_strings) & !grepl("g__uncultured|g__$", tax_strings) ~ sub(".*g__([^;]*).*", "\\1", tax_strings),
        grepl("f__[^;]*", tax_strings) & !grepl("f__uncultured|f__$", tax_strings) ~ paste0("other ", stringr::str_extract(tax_strings, "f__[^;]*") %>% sub("f__", "", .)),
        grepl("o__[^;]*", tax_strings) & !grepl("o__uncultured|o__$", tax_strings) ~ paste0("other ", stringr::str_extract(tax_strings, "o__[^;]*") %>% sub("o__", "", .)),
        grepl("c__[^;]*", tax_strings) & !grepl("c__uncultured|c__$", tax_strings) ~ paste0("other ", stringr::str_extract(tax_strings, "c__[^;]*") %>% sub("c__", "", .)),
        grepl("p__[^;]*", tax_strings) & !grepl("p__uncultured|p__$", tax_strings) ~ paste0("other ", stringr::str_extract(tax_strings, "p__[^;]*") %>% sub("p__", "", .)),
        TRUE ~ "Unclassified"
    )

    table_abundance <- phy.ra.complete %>%
        dplyr::mutate(abun = rowMeans(.)) %>%
        tibble::rownames_to_column(var = "OTUID") %>%
        dplyr::left_join(table_tax, by = "OTUID") %>%
        dplyr::mutate(
            taxonomy2 = taxonomy,
            taxonomy = tax_labels[match(taxonomy, tax_strings)]
        ) %>%
        {
            if (exclude_unclassified) dplyr::filter(., taxonomy != "Unclassified") else .
        } %>%
        dplyr::arrange(-abun) %>%
        dplyr::slice(seq_len(top_n)) %>%
        tidyr::separate(taxonomy2,
            into = c("dominio", "phylum", "clase", "orden", "familia", "genero", "especie"),
            sep = ";", fill = "right", extra = "merge"
        ) %>%
        dplyr::mutate(asv = paste0(feature_prefix, dplyr::row_number())) %>%
        tidyr::unite("taxa", asv, taxonomy, remove = FALSE) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), ~ trimws(.))) %>%
        dplyr::mutate(phylum = sub("^p__", "", phylum)) %>%
        dplyr::mutate(phylum = dplyr::case_when(
            phylum == "Proteobacteria" ~ "Pseudomonadata",
            phylum == "Firmicutes" ~ "Bacillota",
            phylum == "Actinobacteriota" ~ "Actinomycetota",
            phylum == "Bacteroidetes" ~ "Bacteroidota",
            TRUE ~ as.character(phylum)
        )) %>%
        dplyr::mutate(phylum = ifelse(is.na(phylum) | phylum == "", "Unclassified", phylum))
    warning(
        "Note: Some bacterial phylum names have been updated to match NCBI's revised taxonomy:\n",
        "\nReference: https://ncbiinsights.ncbi.nlm.nih.gov/2021/12/10/ncbi-taxonomy-prokaryote-phyla-added/"
    )

    if (save_table) {
        utils::write.table(table_abundance,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    ordered_taxa <- table_abundance$taxa

    heat <- table_abundance %>%
        tibble::column_to_rownames(var = "taxa") %>%
        t() %>%
        as.data.frame() %>%
        dplyr::slice(-1) %>%
        tibble::rownames_to_column(var = OTU_ID) %>%
        dplyr::inner_join(metadata, by = OTU_ID)

    heatmap <- heat %>%
        tibble::column_to_rownames(var = OTU_ID) %>%
        dplyr::select(all_of(ordered_taxa)) %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(.))) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(findInterval(
            ., c(0.001, 0.005, 0.01, 0.10, 0.20, 1, 2, 5, 10, 25, 50, 75),
            left.open = TRUE
        ))))

    heatmap <- as.matrix(heatmap)
    rownames(heatmap) <- .mbm_shorten_labels(ordered_taxa, max_label_length,
        composite = composite_names
    )

    if (!cluster) {
        row_order <- ordered_taxa
    } else {
        row_order <- NULL
    }

    annotation_rows <- table_abundance %>%
        dplyr::select(OTUID, phylum) %>%
        tibble::column_to_rownames(var = "OTUID")
    rownames(annotation_rows) <- rownames(heatmap)

    conditions <- Filter(Negate(is.null), list(condition1, condition2, condition3))

    if (length(conditions) > 0) {
        annotation_columns <- heat %>%
            dplyr::select(all_of(unlist(conditions)))
        annotation_columns[] <- lapply(annotation_columns, as.character)
        rownames(annotation_columns) <- colnames(heatmap)
    } else {
        annotation_columns <- data.frame()
    }

    gp_rows <- grid::gpar(fontsize = 12, fontface = "italic", fontfamily = "serif")
    gp_labels <- grid::gpar(fontsize = 12, fontfamily = "serif")
    gp_names <- grid::gpar(fontsize = 12, fontface = "bold", fontfamily = "serif")
    text_mm <- function(x, gp) {
        grid::convertWidth(ComplexHeatmap::max_text_width(x, gp = gp), "mm", valueOnly = TRUE)
    }
    min_row_mm <- 12 * 0.3528 * 1.15
    flexible <- !draw && is.null(cell_size)
    if (flexible) {
        cell_w <- cell_h <- 4
    } else if (is.null(cell_size)) {
        dev_mm <- grDevices::dev.size("in") * 25.4
        legend_mm <- 10 + max(
            text_mm(c(
                unique(table_abundance$phylum),
                unlist(lapply(annotation_columns, unique))
            ), gp_labels),
            text_mm(c(
                "Phylum", unlist(conditions),
                name_legend_condition1, name_legend_condition2,
                name_legend_condition3
            ), gp_names)
        )
        side_mm <- max(
            text_mm(rownames(heatmap), gp_rows),
            if (length(conditions) > 0) text_mm(unlist(conditions), gp_names) else 0
        )
        avail_w <- dev_mm[1] - side_mm - legend_mm - (if (cluster) 10 else 0) - 15
        avail_h <- dev_mm[2] - 25 - text_mm("Phylum", gp_names) - 8
        n_cond <- length(conditions)
        cell_w <- min(avail_w / (ncol(heatmap) + 1), avail_h / (nrow(heatmap) + n_cond))
        cell_h <- cell_w
        if (cell_h < min_row_mm) {
            cell_h <- min(min_row_mm, avail_h / (nrow(heatmap) + n_cond))
        }
        cell_w <- max(cell_w, 0.5)
    } else {
        cell_w <- cell_h <- cell_size
    }
    if (is.null(annotation_height)) annotation_height <- cell_h
    ann_name_size <- min(12, annotation_height / 0.3528 / 1.1)

    my_palette <- scales::viridis_pal(option = "C", direction = -1)(13)

    unique_phyla <- unique(annotation_rows$phylum)
    c5.phylum <- rep_len(scales::brewer_pal(palette = "Set2")(8), length(unique_phyla))
    cols_phyl <- list(Phylum = setNames(c5.phylum, unique_phyla))

    annphylum <- ComplexHeatmap::rowAnnotation(
        Phylum = annotation_rows$phylum,
        show_legend = FALSE,
        show_annotation_name = TRUE,
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold", fontfamily = "serif"),
        gp = grid::gpar(col = "white"),
        col = cols_phyl,
        simple_anno_size = grid::unit(cell_w, "mm")
    )

    heatmap_annotations <- list()
    legend_list <- list(lgd1 = ComplexHeatmap::Legend(
        at = unique(annotation_rows$phylum),
        legend_gp = grid::gpar(fill = c5.phylum),
        title = "Phylum",
        labels_gp = grid::gpar(fontsize = 12, fontfamily = "serif")
    ))

    if (!is.null(condition1)) {
        unique_vals <- sort(unique(annotation_columns[[condition1]]))

        if (is.null(colors_condition1)) {
            colors_condition1 <- if (length(unique_vals) == 2) {
                .mbm_colors_2group
            } else {
                rep_len(.mbm_colors, length(unique_vals))
            }
        } else if (is.null(names(colors_condition1)) && length(colors_condition1) < length(unique_vals)) {
            colors_condition1 <- rep_len(colors_condition1, length(unique_vals))
        }

        color_mapping <- .mbm_match_colors(colors_condition1, unique_vals)

        heatmap_annotations$ann1 <- ComplexHeatmap::HeatmapAnnotation(
            df = annotation_columns[condition1],
            which = "column",
            col = setNames(list(color_mapping), condition1),
            show_legend = FALSE,
            show_annotation_name = TRUE,
            annotation_name_gp = grid::gpar(fontsize = ann_name_size, fontface = "bold", fontfamily = "serif"),
            gp = grid::gpar(col = "white"),
            simple_anno_size = grid::unit(annotation_height, "mm")
        )

        legend_list$lgd2 <- ComplexHeatmap::Legend(
            at = names(color_mapping),
            legend_gp = grid::gpar(fill = color_mapping),
            title = if (!is.null(name_legend_condition1)) name_legend_condition1 else condition1,
            labels_gp = grid::gpar(fontsize = 12, fontfamily = "serif")
        )
    }
    if (!is.null(condition2)) {
        unique_vals <- sort(unique(annotation_columns[[condition2]]))

        if (is.null(colors_condition2)) {
            colors_condition2 <- rep_len(.mbm_colors_safe, length(unique_vals))
        } else if (is.null(names(colors_condition2)) && length(colors_condition2) < length(unique_vals)) {
            colors_condition2 <- rep_len(colors_condition2, length(unique_vals))
        }

        color_mapping <- .mbm_match_colors(colors_condition2, unique_vals)

        heatmap_annotations$ann2 <- ComplexHeatmap::HeatmapAnnotation(
            df = annotation_columns[condition2],
            which = "column",
            col = setNames(list(color_mapping), condition2),
            show_legend = FALSE,
            show_annotation_name = TRUE,
            annotation_name_gp = grid::gpar(fontsize = ann_name_size, fontface = "bold", fontfamily = "serif"),
            gp = grid::gpar(col = "white"),
            simple_anno_size = grid::unit(annotation_height, "mm")
        )

        legend_list$lgd3 <- ComplexHeatmap::Legend(
            at = names(color_mapping),
            legend_gp = grid::gpar(fill = color_mapping),
            title = if (!is.null(name_legend_condition2)) name_legend_condition2 else condition2,
            labels_gp = grid::gpar(fontsize = 12, fontfamily = "serif")
        )
    }

    if (!is.null(condition3)) {
        unique_vals <- sort(unique(annotation_columns[[condition3]]))

        if (is.null(colors_condition3)) {
            mbm_shift3 <- c(.mbm_colors[5:8], .mbm_colors[seq_len(4)])
            colors_condition3 <- rep_len(mbm_shift3, length(unique_vals))
        } else if (is.null(names(colors_condition3)) && length(colors_condition3) < length(unique_vals)) {
            colors_condition3 <- rep_len(colors_condition3, length(unique_vals))
        }

        color_mapping <- .mbm_match_colors(colors_condition3, unique_vals)

        heatmap_annotations$ann3 <- ComplexHeatmap::HeatmapAnnotation(
            df = annotation_columns[condition3],
            which = "column",
            col = setNames(list(color_mapping), condition3),
            show_legend = FALSE,
            show_annotation_name = TRUE,
            annotation_name_gp = grid::gpar(fontsize = ann_name_size, fontface = "bold", fontfamily = "serif"),
            gp = grid::gpar(col = "white"),
            simple_anno_size = grid::unit(annotation_height, "mm")
        )

        legend_list$lgd4 <- ComplexHeatmap::Legend(
            at = names(color_mapping),
            legend_gp = grid::gpar(fill = color_mapping),
            title = if (!is.null(name_legend_condition3)) name_legend_condition3 else condition3,
            labels_gp = grid::gpar(fontsize = 12, fontfamily = "serif")
        )
    }

    top_annotation <- do.call(c, unname(heatmap_annotations))

    pd.legends <- do.call(ComplexHeatmap::packLegend, legend_list)

    heats <- ComplexHeatmap::Heatmap(
        heatmap,
        col = my_palette,
        heatmap_legend_param = list(
            direction = "horizontal",
            labels_gp = grid::gpar(fontsize = 12, fontfamily = "serif"),
            title_gp = grid::gpar(fontsize = 14, fontface = "bold", fontfamily = "serif"),
            legend_gp = grid::gpar(fontsize = 12, fontfamily = "serif"),
            title = "Relative abundance (%)",
            title_position = "topcenter",
            at = 0:12
        ),
        rect_gp = grid::gpar(col = "black", lwd = 0.5),
        row_names_gp = gp_rows,
        row_names_max_width = ComplexHeatmap::max_text_width(rownames(heatmap), gp = gp_rows),
        column_names_gp = grid::gpar(fontsize = 12, fontfamily = "serif"),
        cluster_columns = FALSE,
        cluster_rows = cluster,
        width = if (!flexible) grid::unit(ncol(heatmap) * cell_w, "mm"),
        height = if (!flexible) grid::unit(nrow(heatmap) * cell_h, "mm"),
        row_order = row_order,
        show_column_names = show_column_names,
        show_heatmap_legend = TRUE,
        top_annotation = if (length(heatmap_annotations) > 0) top_annotation else NULL,
        left_annotation = annphylum
    )

    draw_heatmap <- function() {
        ComplexHeatmap::draw(
            heats,
            heatmap_legend_side = "top",
            annotation_legend_side = "right",
            merge_legend = FALSE,
            annotation_legend_list = pd.legends,
            background = "transparent"
        )
    }

    if (draw) {
        grid::grid.newpage()
        draw_heatmap()
        heatmap_output <- grid::grid.grab(wrap.grobs = TRUE)
    } else {
        heatmap_output <- grid::grid.grabExpr(draw_heatmap())
    }
    class(heatmap_output) <- c("mbm_heatmap", class(heatmap_output))
    return(invisible(heatmap_output))
}

#' @export
#' @noRd
print.mbm_heatmap <- function(x, ...) {
    grid::grid.newpage()
    grid::grid.draw(x)
    invisible(x)
}
