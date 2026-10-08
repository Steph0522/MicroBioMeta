#' Compare taxon abundance ratios between two conditions
#'
#' Computes relative abundances from a taxonomic abundance table and compares
#' two experimental conditions by calculating a directional abundance ratio
#' for each taxon.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. Its first column
#'   must hold the sample IDs (the column names of `table`).
#' @param group_col Character. Name of the column in \code{metadata} that
#'   defines the groups. \code{condition_A} and \code{condition_B} are two of its
#'   values.
#' @param condition_A Character. Name of the first condition to compare.
#' @param condition_B Character. Name of the second condition to compare.
#' @param taxonomy_db Character. Database the taxonomy strings come from:
#'   \code{"silva"} (default), \code{"gg2"} (Greengenes2, also \code{"gg"}),
#'   \code{"unite"} or \code{"Kraken2"} (also \code{"kraken"}). Case-insensitive.
#' @param top_n Integer. Number of taxa with the highest mean abundance to
#'   show. Default \code{30}.
#' @param level Character. Taxonomic level: \code{"kingdom"}, \code{"phylum"},
#'   \code{"class"}, \code{"order"}, \code{"family"}, \code{"genus"} (default) or
#'   \code{"species"}. Case-insensitive.
#' @param x_axis_title Character. Title of the x-axis (taxon names). Default
#'   \code{"Taxon"}.
#' @param legend_title Character. Title of the legend and of the axis showing
#'   the dominant condition. If \code{NULL} (default), the name of
#'   \code{group_col} is used.
#' @param group_colors Optional character vector of two colors, one per
#'   condition. If \code{NULL} (default), the package's orange/blue pair is used.
#' @param save_table Logical. If \code{TRUE}, saves the ratio table as a
#'   tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"ratios_bubble_table.txt"}.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' ratios_bubble_plot(
#'     table = table,
#'     metadata = metadata,
#'     group_col = "Location",
#'     condition_A = "Rhizosphere",
#'     condition_B = "Roots",
#'     taxonomy_db = "silva",
#'     level = "genus",
#'     top_n = 20
#' )
ratios_bubble_plot <- function(table,
                               metadata,
                               group_col,
                               condition_A,
                               condition_B,
                               taxonomy_db = "silva",
                               top_n = 30,
                               level = "genus",
                               x_axis_title = "Taxon",
                               legend_title = NULL,
                               group_colors = NULL,
                               save_table = FALSE,
                               table_filename = "ratios_bubble_table.txt") {
    taxonomy_db <- .mbm_taxonomy_db(taxonomy_db)
    level <- .mbm_check_level(level)

    metadata_sub <- metadata %>%
        dplyr::filter(.data[[group_col]] %in% c(condition_A, condition_B)) %>%
        dplyr::select(SampleID = 1, Condition = dplyr::all_of(group_col))

    samples <- metadata_sub$SampleID

    abundance_raw <- table %>%
        dplyr::select(taxonomy, dplyr::all_of(samples))

    abundance_raw <- .mbm_truncate_taxonomy(abundance_raw, level, taxonomy_db)
    abundance_raw <- .mbm_label_taxa(abundance_raw, level, taxonomy_db)
    abundance_raw$taxonomy <- .mbm_composite_genus(abundance_raw$taxonomy)

    abundance_raw <- abundance_raw %>%
        dplyr::group_by(taxonomy) %>%
        dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

    abundance_rel <- abundance_raw %>%
        tibble::column_to_rownames("taxonomy") %>%
        sweep(2, colSums(.), FUN = "/") %>%
        as.data.frame() %>%
        tibble::rownames_to_column("taxonomy")

    long_data <- abundance_rel %>%
        tidyr::pivot_longer(-taxonomy, names_to = "SampleID", values_to = "Abundance") %>%
        dplyr::left_join(metadata_sub, by = "SampleID")

    summary_data <- long_data %>%
        dplyr::group_by(taxonomy, Condition) %>%
        dplyr::summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
        tidyr::pivot_wider(
            names_from = Condition,
            values_from = MeanAbundance,
            values_fill = 0
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            Ratio = dplyr::if_else(
                .data[[condition_A]] > .data[[condition_B]],
                (.data[[condition_A]] - .data[[condition_B]]) / .data[[condition_B]],
                (.data[[condition_B]] - .data[[condition_A]]) / .data[[condition_A]]
            ),
            Dominant = dplyr::if_else(.data[[condition_A]] > .data[[condition_B]],
                condition_A,
                condition_B
            ),
            MeanAbund = mean(c(.data[[condition_A]], .data[[condition_B]]), na.rm = TRUE)
        ) %>%
        dplyr::ungroup()

    top_taxa <- summary_data %>%
        dplyr::slice_max(order_by = MeanAbund, n = top_n) %>%
        dplyr::arrange(desc(Ratio))

    if (save_table) {
        utils::write.table(top_taxa,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    legend_name <- if (!is.null(legend_title)) legend_title else group_col
    ggplot2::ggplot(
        top_taxa,
        ggplot2::aes(
            x = reorder(taxonomy, Ratio),
            y = Dominant,
            size = MeanAbund,
            fill = Dominant
        )
    ) +
        ggplot2::geom_point(shape = 21, color = "black") +
        ggplot2::scale_fill_manual(
            values = if (!is.null(group_colors)) group_colors else .mbm_colors_2group
        ) +
        ggplot2::scale_size(range = c(3, 10)) +
        ggplot2::labs(
            x    = x_axis_title,
            y    = legend_name,
            size = "Ratio",
            fill = legend_name
        ) +
        ggplot2::coord_flip() +
        .mbm_theme(
            legend_position = "right",
            extra = ggplot2::theme(
                axis.text.y = ggplot2::element_text(face = "italic")
            )
        )
}
