#'Compare taxon abundance ratios between two conditions
#'
#' Computes relative abundances from a taxonomic abundance table and compares
#' two experimental conditions by calculating a directional abundance ratio
#' for each taxon. Taxa are ranked by mean abundance and visualized as a bubble
#' plot, where bubble size represents mean relative abundance and color indicates
#' the dominant condition.
#'
#' @param table A data frame containing a taxonomic abundance table with a
#'  taxonomy column and sample columns with numeric counts.
#' @param metadata A data frame containing sample metadata. The first column must
#'   correspond to sample IDs and include a column defining the experimental
#'   conditions.
#' @param condition_col Character. Name of the metadata column defining the
#'   experimental condition.
#' @param condition_A Character. Name of the first condition to compare.
#' @param condition_B Character. Name of the second condition to compare.
#' @param taxonomy_db Character. Taxonomic database used for annotation.
#'  ("silva","Kraken2").
#' @param top_n Integer. Number of taxa with the highest mean abundance to display.
#' @param level Character. Taxonomic level to use for comparison
#'   (e.g. "phylum", "genus", "species").
#' @param x_axis_title Character. Label for the x-axis (taxon names).
#' @param group_colors Character vector of colors used to represent the dominant
#'   condition.
#' @param save_table Logical. If \code{TRUE}, saves the underlying ratio table
#'   to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"ratios_bubble_table.txt"}.
#'
#' @return A ggplot2 object showing abundance ratios between the two
#'   conditions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ratio_plot(
#'   table         = table,
#'   metadata      = metadata,
#'   condition_col = "metodo",
#'   condition_A   = "kit",
#'   condition_B   = "fenol",
#'   level         = "genus",
#'   top_n         = 20
#' )
#' }

ratio_plot <- function(table,
                                 metadata,
                                 condition_col,
                                 condition_A,
                                 condition_B,
                                 taxonomy_db = "silva",
                                 top_n = 30,
                                 level = "genus",
                                 x_axis_title = "Taxon",
                                 group_colors = NULL,
                                 save_table = FALSE,
                                 table_filename = "ratios_bubble_table.txt") {
  # Filtrar metadatos a las condiciones deseadas
  metadata_sub <- metadata %>%
    dplyr::filter(.data[[condition_col]] %in% c(condition_A, condition_B)) %>%
    dplyr::select(SampleID = 1, Condition = dplyr::all_of(condition_col))

  samples <- metadata_sub$SampleID

  # Separar tabla de taxonomía y abundancias
  abundance_raw <- table %>%
    dplyr::select(taxonomy, dplyr::all_of(samples))

  # Colapsar por taxonomía si hay duplicados
  abundance_raw <- abundance_raw %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  # Corregir taxonomía según base y nivel
  abundance_raw <- abundance_raw %>%
    dplyr::mutate(
      taxonomy = dplyr::case_when(
        taxonomy_db == "Kraken2" & level == "specie" ~ dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*;.*s__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$|s__uncultured|s__$", taxonomy) ~
            paste0(
              stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .),
              " ",
              stringr::str_extract(taxonomy, "s__[^;]*") %>% sub("s__", "", .)
            ),
          grepl("g__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$", taxonomy) ~
            paste0(
              "other ",
              stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .)
            ),
          grepl("f__[^;]*", taxonomy) &
            !grepl("f__uncultured|f__$", taxonomy) ~
            paste0(
              "other ",
              stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)
            ),
          grepl("o__[^;]*", taxonomy) &
            !grepl("o__uncultured|o__$", taxonomy) ~
            paste0(
              "other ",
              stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)
            ),
          grepl("c__[^;]*", taxonomy) &
            !grepl("c__uncultured|c__$", taxonomy) ~
            paste0(
              "other ",
              stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)
            ),
          grepl("p__[^;]*", taxonomy) &
            !grepl("p__uncultured|p__$", taxonomy) ~
            paste0(
              "other ",
              stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)
            ),
          TRUE ~ "Unclassified"
        ),

        taxonomy_db %in% c("silva", "Kraken2") &
          level == "genus" ~ dplyr::case_when(
            taxonomy == "Other" ~ "Other",
            grepl("g__[^;]*", taxonomy) &
              !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
            grepl("f__[^;]*", taxonomy) &
              !grepl("f__uncultured|f__$", taxonomy) ~ paste0(
                "other ",
                stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)
              ),
            grepl("o__[^;]*", taxonomy) &
              !grepl("o__uncultured|o__$", taxonomy) ~ paste0(
                "other ",
                stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)
              ),
            grepl("c__[^;]*", taxonomy) &
              !grepl("c__uncultured|c__$", taxonomy) ~ paste0(
                "other ",
                stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)
              ),
            grepl("p__[^;]*", taxonomy) &
              !grepl("p__uncultured|p__$", taxonomy) ~ paste0(
                "other ",
                stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)
              ),
            TRUE ~ "Unclassified"
          ),

        taxonomy_db %in% c("silva", "Kraken2") &
          level == "phylum" ~ dplyr::case_when(
            taxonomy == "Other" ~ "Other",
            grepl("p__[^;]*", taxonomy) ~ stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
            TRUE ~ "Unclassified"
          ),

        TRUE ~ taxonomy
      )
    ) %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  # Calcular abundancias relativas por muestra
  abundance_rel <- abundance_raw %>%
    tibble::column_to_rownames("taxonomy") %>%
    sweep(2, colSums(.), FUN = "/") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("taxonomy")

  # Convertir a formato largo
  long_data <- abundance_rel %>%
    tidyr::pivot_longer(-taxonomy, names_to = "SampleID", values_to = "Abundance") %>%
    dplyr::left_join(metadata_sub, by = "SampleID")

  # Calcular abundancia media por condición
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
                         condition_B),
      MeanAbund = mean(c(.data[[condition_A]], .data[[condition_B]]), na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  # Seleccionar top_n taxones más abundantes
  top_taxa <- summary_data %>%
    dplyr::slice_max(order_by = MeanAbund, n = top_n) %>%
    dplyr::arrange(desc(Ratio))

  if (save_table) {
    utils::write.table(top_taxa, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # Bubble plot
  ggplot2::ggplot(top_taxa,
         ggplot2::aes(
           x = reorder(taxonomy, Ratio),
           y = Dominant,
           size = MeanAbund,
           fill = Dominant
         )) +
    ggplot2::geom_point(shape = 21, color = "black") +
    ggplot2::scale_fill_manual(
      values = if (!is.null(group_colors)) group_colors else .mbm_colors
    ) +
    ggplot2::scale_size(range = c(3, 10)) +
    ggplot2::labs(
      x    = x_axis_title,
      y    = "Condition",
      size = "Ratio",
      fill = "Condition"
    ) +
    ggplot2::coord_flip() +
    .mbm_theme(legend_position = "right")
  
}
