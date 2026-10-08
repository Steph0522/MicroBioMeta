#' Random forest feature importance plot
#'
#' Fits a random forest (randomForest) that predicts a metadata variable from
#' the abundances, and plots the most important features as a lollipop plot
#' colored by phylum.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. Its first column
#'   must hold the sample IDs (the column names of `table`).
#' @param variable_to_predict Character. Name of the column in \code{metadata}
#'   to predict (a categorical variable).
#' @param top_n Integer. Number of most important features to show. Default
#'   \code{15}.
#' @param group_colors Optional character vector of colors, one per phylum
#'   (recycled if shorter). If \code{NULL} (default), the package palette is
#'   used.
#' @param title Character. Plot title. \code{NULL} (default) shows no title.
#' @param size Numeric. Size of the lollipop points. Default \code{8}.
#' @param save_table Logical. If \code{TRUE}, saves the feature-importance
#'   table as a tab-delimited file. Default \code{FALSE}.
#' @param x_axis_title Character. Title of the x-axis (feature importance).
#'   Default \code{"Feature importance (MeanDecreaseGini)"}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"randomforest_importance.txt"}.
#'
#' @details The random forest is random; call \code{set.seed()} before the
#'   function to make the result reproducible.
#'
#' @return A ggplot object.
#' @importFrom randomForest randomForest importance
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' random_forest_lollipop_plot(
#'     table               = table,
#'     metadata            = metadata,
#'     variable_to_predict = "Location",
#'     top_n               = 20,
#'     size                = 6
#' )
random_forest_lollipop_plot <- function(table,
                                        metadata,
                                        top_n = 15,
                                        size = 8,
                                        variable_to_predict,
                                        group_colors = NULL,
                                        title = NULL,
                                        save_table = FALSE,
                                        table_filename = "randomforest_importance.txt",
                                        x_axis_title = "Feature importance (MeanDecreaseGini)") {
    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")

    if (is.null(title)) {
        title <- sprintf("Top %d most important features (Random Forest)", top_n)
    }

    ncols <- ncol(table)
    table_numeric <- table[-ncols]
    taxonomy <- table[ncols]

    if (ncol(table_numeric) < nrow(table_numeric)) {
        table_numeric <- data.frame(t(table_numeric), check.names = FALSE)
    }

    metadata_filtered <- .mbm_align_metadata(rownames(table_numeric), metadata)
    otu_filtered <- table_numeric[metadata_filtered[[1]], , drop = FALSE]

    response <- as.factor(metadata_filtered[[variable_to_predict]])

    if (length(unique(response)) < 2) {
        stop("The variable to predict must contain at least two classes.")
    }

    rf_model <- randomForest::randomForest(
        x = otu_filtered,
        y = response,
        importance = TRUE,
        ntree = 500
    )

    importance_df <- randomForest::importance(rf_model)
    importance_df <- as.data.frame(importance_df)

    importance_df$ASV <- rownames(importance_df)
    top_asvs <- importance_df %>%
        dplyr::arrange(desc(MeanDecreaseGini)) %>%
        head(top_n)

    top_asvs <- dplyr::left_join(
        top_asvs,
        data.frame(ASV = colnames(as.data.frame(table_numeric)), taxonomy = taxonomy),
        by = "ASV"
    ) %>%
        dplyr::mutate(taxonomy_original = taxonomy)

    top_asvs <- top_asvs %>%
        tidyr::separate(
            col = taxonomy_original,
            into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
            sep = ";",
            fill = "right"
        ) %>%
        dplyr::mutate(across(everything(), ~ trimws(.)))

    top_asvs.modified <- top_asvs %>%
        dplyr::mutate(
            taxonomy = dplyr::case_when(
                grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~
                    sub(".*g__([^;]*).*", "\\1", taxonomy),
                grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~
                    paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
                grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~
                    paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
                grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~
                    paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
                grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~
                    paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
                TRUE ~ "Unclassified"
            ),
            Phylum = dplyr::case_when(
                Phylum == "p__Proteobacteria" ~ "Pseudomonadota",
                Phylum == "p__Actinobacteriota" ~ "Actinomycetota",
                Phylum == "p__Bacteroidetes" ~ "Bacteroidota",
                Phylum == "p__Firmicutes" ~ "Bacillota",
                TRUE ~ sub("^p__", "", as.character(Phylum))
            ),
            row_id = paste0("row", dplyr::row_number())
        )
    updated_phyla <- c(
        "p__Proteobacteria", "p__Actinobacteriota",
        "p__Bacteroidetes", "p__Firmicutes"
    )
    present_phyla <- intersect(updated_phyla, top_asvs$Phylum)

    if (length(present_phyla) > 0) {
        old_new <- data.frame(
            Old = c("Proteobacteria", "Actinobacteriota", "Bacteroidetes", "Firmicutes"),
            New = c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Bacillota")
        )

        changes <- old_new[old_new$Old %in% sub("^p__", "", present_phyla), ]

        warning(
            "Note: Some bacterial phylum names have been updated to match NCBI's revised taxonomy:\n",
            paste(sprintf("  - '%s' changed to '%s'", changes$Old, changes$New), collapse = "\n"),
            "\nReference: https://ncbiinsights.ncbi.nlm.nih.gov/2021/12/10/ncbi-taxonomy-prokaryote-phyla-added/"
        )
    }
    n_phyla <- length(unique(top_asvs.modified$Phylum))
    if (!is.null(group_colors)) {
        if (length(group_colors) < n_phyla) {
            warning(sprintf(
                "Custom palette has %d colors but %d are needed. Recycling palette.",
                length(group_colors), n_phyla
            ))
            group_colors <- rep_len(group_colors, n_phyla)
        }
        fill_colors <- group_colors
    } else {
        fill_colors <- rep_len(.mbm_colors, n_phyla)
    }

    top_asvs.modified$MeanDecreaseGini <- as.numeric(top_asvs.modified$MeanDecreaseGini)

    if (save_table) {
        utils::write.table(top_asvs.modified,
            file = table_filename, sep = "\t",
            quote = FALSE, row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    row_labels <- setNames(top_asvs.modified$taxonomy, top_asvs.modified$row_id)
    lollipop <- ggplot2::ggplot(
        top_asvs.modified,
        ggplot2::aes(
            x = reorder(row_id, MeanDecreaseGini),
            y = MeanDecreaseGini,
            fill = Phylum
        )
    ) +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = reorder(row_id, MeanDecreaseGini),
                xend = reorder(row_id, MeanDecreaseGini),
                y = 0,
                yend = MeanDecreaseGini
            ),
            color = "black",
            linewidth = 1
        ) +
        ggplot2::ylab(x_axis_title) +
        ggplot2::geom_point(size = size, shape = 21, color = "black") +
        ggplot2::scale_fill_manual(values = fill_colors) +
        ggplot2::scale_x_discrete(labels = row_labels) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = title) +
        .mbm_theme(
            legend_position = "right",
            extra = ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(
                    size = 12, color = "black",
                    face = "italic"
                )
            )
        )

    return(lollipop)
}
