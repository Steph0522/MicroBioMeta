#' Venn diagram of taxa shared between sample groups
#'
#' Creates a Venn diagram of the taxa shared between groups of samples, with
#' ggvenn or ggVennDiagram.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. Its first column
#'   must hold the sample IDs (the column names of `table`).
#' @param merge_by Character. Name of the column in \code{metadata} that
#'   defines the groups (one set of taxa per group).
#' @param selected_samples Optional character vector of the sample IDs to use.
#'   If \code{NULL} (default), all samples are used.
#' @param min_prevalence Numeric (0-1). Minimum fraction of the samples of a
#'   group where a taxon must be present to count it in that group. Default
#'   \code{0} (present in at least one sample).
#' @param title Character. Plot title. \code{NULL} (default) shows no title.
#' @param method Character. Package used to draw the Venn diagram: `"ggvenn"`
#'   (default) or `"ggVennDiagram"`. Case-insensitive.
#' @param group_colors Optional character vector of colors, one per group,
#'   named after the groups or in their order. If \code{NULL} (default), the
#'   colorblind-friendly Okabe-Ito palette is used.
#' @param save_table Logical. If \code{TRUE}, saves the taxa of each group
#'   (long format) as a tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"venn_taxa_sets.txt"}.
#'
#' @return A ggplot object.
#' @importFrom ggvenn ggvenn
#' @importFrom ggVennDiagram ggVennDiagram
#' @export
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' venn_plot(
#'     table          = table,
#'     metadata       = metadata,
#'     merge_by       = "Location",
#'     min_prevalence = 0
#' )
#'
#' ## Filtering taxa by prevalence
#' venn_plot(
#'     table          = table,
#'     metadata       = metadata,
#'     merge_by       = "Location",
#'     min_prevalence = 0.2
#' )
#'
#' ## Custom colors
#' venn_plot(
#'     table          = table,
#'     metadata       = metadata,
#'     merge_by       = "Location",
#'     min_prevalence = 0,
#'     group_colors   = c("#1B9E77", "#D95F02")
#' )
venn_plot <- function(table, metadata, merge_by,
                      selected_samples = NULL, min_prevalence = 0,
                      title = NULL, method = "ggvenn",
                      group_colors = NULL,
                      save_table = FALSE,
                      table_filename = "venn_taxa_sets.txt") {
    table <- as.data.frame(table)
    metadata <- as.data.frame(metadata)

    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")

    if (!merge_by %in% colnames(metadata)) stop("Group or merge column is not in the metadata file.")

    common_samples <- intersect(colnames(table), metadata[[1]])
    table <- table[, c(common_samples), drop = FALSE]
    metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
    rownames(metadata) <- NULL

    if (!is.null(selected_samples)) {
        common_samples <- intersect(selected_samples, colnames(table))
        table <- table[, c(common_samples), drop = FALSE]
        metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
    }

    metadata_split <- split(metadata[[1]], metadata[[merge_by]])
    metadata_split <- metadata_split[vapply(metadata_split, length, integer(1)) > 0]
    num_groups <- length(metadata_split)

    use_manual <- !is.null(group_colors)
    if (!use_manual) {
        group_colors <- if (num_groups == 2) {
            .mbm_colors_2group
        } else {
            rep_len(.mbm_colors, num_groups)
        }
    } else if (!is.null(names(group_colors))) {
        group_colors <- unname(.mbm_match_colors(
            group_colors, names(metadata_split),
            "group_colors"
        ))
    } else {
        group_colors <- rep(group_colors, length.out = num_groups)
    }

    lista <- lapply(metadata_split, function(samps) {
        subset <- table[, samps, drop = FALSE]
        subset_core <- if (min_prevalence > 0) {
            subset_core <- subset[rowMeans(subset > 0) >= min_prevalence, ]
        } else {
            subset_core <- subset[rowSums(subset) != 0, ]
        }
        rownames(subset_core)
    })
    names(lista) <- names(metadata_split)

    if (save_table) {
        venn_table <- utils::stack(lista)
        colnames(venn_table) <- c("taxon_id", "group")
        utils::write.table(venn_table,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    if (tolower(method) == "ggvenndiagram") {
        venn_plot <- ggVennDiagram::ggVennDiagram(
            lista,
            label_alpha = 0,
            set_color   = group_colors,
            edge_size   = 1
        ) +
            ggplot2::scale_fill_gradient(
                low = "white", high = "grey60",
                na.value = NA, name = "Count"
            )
    } else if (tolower(method) == "ggvenn") {
        venn_plot <- ggvenn::ggvenn(lista,
            fill_color = group_colors,
            fill_alpha = 0.5,
            stroke_color = "grey30",
            set_name_size = 5,
            text_size = 4
        )
    } else {
        stop("Method must be 'ggvenn' or 'ggVennDiagram'.")
    }

    venn_plot <- venn_plot +
        ggplot2::theme_void(base_family = "serif") +
        ggplot2::theme(
            legend.position = "right",
            legend.text = ggplot2::element_text(size = 12, color = "black"),
            legend.title = ggplot2::element_text(
                size = 14, face = "bold",
                color = "black"
            )
        )

    if (!is.null(title)) {
        venn_plot <- venn_plot +
            ggplot2::labs(title = title) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    hjust = 0.5, face = "bold",
                    size = 14, color = "black"
                )
            )
    }
    return(venn_plot)
}
