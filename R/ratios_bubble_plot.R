#' ratio_plot
#'
#' @param table 
#' @param metadata 
#' @param condition_col 
#' @param condition_A 
#' @param condition_B 
#' @param taxonomy_db 
#' @param top_n 
#' @param level 
#' @param x_axis_title 
#' @param fill_palette 
#'
#' @return
#' @export
#'
#' @examples
ratio_plot <- function(table,
                                 metadata,
                                 condition_col,
                                 condition_A,
                                 condition_B,
                                 taxonomy_db = "silva",
                                 top_n = 30,
                                 level = "genus",
                                 x_axis_title = "Taxon",
                                 fill_palette = c("#1f77b4", "#ff7f0e")) {
  library(tidyverse)
  
  # Filtrar metadatos a las condiciones deseadas
  metadata_sub <- metadata %>%
    filter(.data[[condition_col]] %in% c(condition_A, condition_B)) %>%
    select(SampleID = 1, Condition = all_of(condition_col))
  
  samples <- metadata_sub$SampleID
  
  # Separar tabla de taxonomía y abundancias
  abundance_raw <- table %>%
    select(taxonomy, all_of(samples))
  
  # Colapsar por taxonomía si hay duplicados
  abundance_raw <- abundance_raw %>%
    group_by(taxonomy) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  # Corregir taxonomía según base y nivel
  abundance_raw <- abundance_raw %>%
    mutate(
      taxonomy = case_when(
        taxonomy_db == "Kraken2" & level == "specie" ~ case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*;.*s__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$|s__uncultured|s__$", taxonomy) ~
            paste0(
              str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .),
              " ",
              str_extract(taxonomy, "s__[^;]*") %>% sub("s__", "", .)
            ),
          grepl("g__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$", taxonomy) ~
            paste0(
              "other ",
              str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .)
            ),
          grepl("f__[^;]*", taxonomy) &
            !grepl("f__uncultured|f__$", taxonomy) ~
            paste0(
              "other ",
              str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)
            ),
          grepl("o__[^;]*", taxonomy) &
            !grepl("o__uncultured|o__$", taxonomy) ~
            paste0(
              "other ",
              str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)
            ),
          grepl("c__[^;]*", taxonomy) &
            !grepl("c__uncultured|c__$", taxonomy) ~
            paste0(
              "other ",
              str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)
            ),
          grepl("p__[^;]*", taxonomy) &
            !grepl("p__uncultured|p__$", taxonomy) ~
            paste0(
              "other ",
              str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)
            ),
          TRUE ~ "Unclassified"
        ),
        
        taxonomy_db %in% c("silva", "Kraken2") &
          level == "genus" ~ case_when(
            taxonomy == "Other" ~ "Other",
            grepl("g__[^;]*", taxonomy) &
              !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
            grepl("f__[^;]*", taxonomy) &
              !grepl("f__uncultured|f__$", taxonomy) ~ paste0(
                "other ",
                str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)
              ),
            grepl("o__[^;]*", taxonomy) &
              !grepl("o__uncultured|o__$", taxonomy) ~ paste0(
                "other ",
                str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)
              ),
            grepl("c__[^;]*", taxonomy) &
              !grepl("c__uncultured|c__$", taxonomy) ~ paste0(
                "other ",
                str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)
              ),
            grepl("p__[^;]*", taxonomy) &
              !grepl("p__uncultured|p__$", taxonomy) ~ paste0(
                "other ",
                str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)
              ),
            TRUE ~ "Unclassified"
          ),
        
        taxonomy_db %in% c("silva", "Kraken2") &
          level == "phylum" ~ case_when(
            taxonomy == "Other" ~ "Other",
            grepl("p__[^;]*", taxonomy) ~ str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
            TRUE ~ "Unclassified"
          ),
        
        TRUE ~ taxonomy
      )
    ) %>%
    group_by(taxonomy) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  # Calcular abundancias relativas por muestra
  abundance_rel <- abundance_raw %>%
    column_to_rownames("taxonomy") %>%
    sweep(2, colSums(.), FUN = "/") %>%
    as.data.frame() %>%
    rownames_to_column("taxonomy")
  
  # Convertir a formato largo
  long_data <- abundance_rel %>%
    pivot_longer(-taxonomy, names_to = "SampleID", values_to = "Abundance") %>%
    left_join(metadata_sub, by = "SampleID")
  
  # Calcular abundancia media por condición
  summary_data <- long_data %>%
    group_by(taxonomy, Condition) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    pivot_wider(
      names_from = Condition,
      values_from = MeanAbundance,
      values_fill = 0
    ) %>%
    rowwise() %>%
    mutate(
      Ratio = if_else(
        .data[[condition_A]] > .data[[condition_B]],
        (.data[[condition_A]] - .data[[condition_B]]) / .data[[condition_B]],
        (.data[[condition_B]] - .data[[condition_A]]) / .data[[condition_A]]
      ),
      Dominant = if_else(.data[[condition_A]] > .data[[condition_B]],
                         condition_A,
                         condition_B),
      MeanAbund = mean(c(.data[[condition_A]], .data[[condition_B]]), na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Seleccionar top_n taxones más abundantes
  top_taxa <- summary_data %>%
    slice_max(order_by = MeanAbund, n = top_n) %>%
    arrange(desc(Ratio))
  
  # Bubble plot
  ggplot(top_taxa,
         aes(
           x = reorder(taxonomy, Ratio),
           y = Dominant,
           size = MeanAbund,
           fill = Dominant
         )) +
    geom_point(shape = 21, color = "black") +
    scale_fill_manual(values = fill_palette) +
    scale_size(range = c(3, 10)) +
    labs(
      x = x_axis_title,
      y = "Condition",
      size = "Ratio",
      fill = "Condition"
    ) +
    theme_classic(base_size = 12) +
    coord_flip() +
    theme(
      legend.title = ggplot2::element_text(size = 14, color = "black", family = "Times New Roman", face = "bold"),
      legend.text = ggplot2::element_text(size = 12, color = "black", family = "Times New Roman"),
      axis.title.x = ggplot2::element_text(size = 14, color = "black", family = "Times New Roman"),
      axis.title.y = ggplot2::element_text(size = 14, color = "black", family = "Times New Roman"),
      axis.text.x = ggplot2::element_text(size = 12, color = "black", family = "Times New Roman"),
      axis.text.y = ggplot2::element_text(size = 12, color = "black", family = "Times New Roman"))
  
}
