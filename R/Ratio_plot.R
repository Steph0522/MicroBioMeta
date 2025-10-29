ratio_plot2 <- function(table,
                       metadata,
                       condition_col,
                       condition_A,
                       condition_B,
                       taxonomy_db = "silva",
                       top_n = 30,
                       level = "genus",
                       x_axis_title = "Taxon",
                       fill_palette = c("#1f77b4", "#ff7f0e", "#999999"),  # A, B, Neutral
                       x_limits = NULL,
                       neutral_threshold = 1,
                       save_table = TRUE,
                       table_filename = "ratio.txt") {
  library(tidyverse)
  
  # Filtrar metadatos
  metadata_sub <- metadata %>%
    filter(.data[[condition_col]] %in% c(condition_A, condition_B)) %>%
    select(SampleID = 1, Condition = all_of(condition_col))
  
  samples <- metadata_sub$SampleID
  
  # Colapsar por taxón
  abundance_raw <- table %>%
    select(taxonomy, all_of(samples)) %>%
    group_by(taxonomy) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  # Corregir taxonomía
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
            !grepl("g__uncultured|g__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .)),
          grepl("f__[^;]*", taxonomy) &
            !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) &
            !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) &
            !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) &
            !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        ),
        taxonomy_db %in% c("silva", "Kraken2") & level == "genus" ~ case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        ),
        taxonomy_db %in% c("silva", "Kraken2") & level == "phylum" ~ case_when(
          taxonomy == "Other" ~ "Other",
          grepl("p__[^;]*", taxonomy) ~ str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
          TRUE ~ "Unclassified"
        ),
        TRUE ~ taxonomy
      )
    ) %>%
    group_by(taxonomy) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop")
  
  # Abundancia relativa
  abundance_rel <- abundance_raw %>%
    column_to_rownames("taxonomy") %>%
    sweep(2, colSums(.), FUN = "/") %>%
    as.data.frame() %>%
    rownames_to_column("taxonomy")
  
  # Formato largo
  long_data <- abundance_rel %>%
    pivot_longer(-taxonomy, names_to = "SampleID", values_to = "Abundance") %>%
    left_join(metadata_sub, by = "SampleID")
  
  # Cálculo de medias y ratios
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
      MeanAbund = mean(c(.data[[condition_A]], .data[[condition_B]]), na.rm = TRUE) * 100,
      SignedRatio = if_else(Dominant == condition_A, Ratio, -Ratio),
      RatioCategory = case_when(
        SignedRatio > neutral_threshold ~ condition_A,
        SignedRatio < -neutral_threshold ~ condition_B,
        TRUE ~ "Neutral"
      )
    ) %>%
    ungroup()
  
  #guardar tabla
  if (save_table) {
    utils::write.table(
      summary_data,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message(paste("Table saved as:", table_filename))
  }
  
  # Top N por abundancia
  top_taxa <- summary_data %>%
    slice_max(order_by = MeanAbund, n = top_n) %>%
    arrange(desc(MeanAbund))
  
  # Eje x automático
  if (identical(x_limits, "auto")) {
    max_ratio <- ceiling(max(abs(top_taxa$SignedRatio), na.rm = TRUE))
    x_limits <- c(-max_ratio, max_ratio)
  }
  
  # Gráfico
  p <- ggplot(top_taxa,
              aes(
                x = SignedRatio,
                y = reorder(taxonomy, MeanAbund),
                fill = RatioCategory
              )) +
    geom_point(shape = 21, color = "black", alpha = 0.85, size = 5) +  # tamaño fijo
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = setNames(fill_palette, c(condition_A, condition_B, "Neutral"))) +
    labs(
      x = "Ratio",
      y = x_axis_title,
      fill = "Dominant Condition") +
    coord_cartesian(xlim = x_limits) +
    theme_test(base_size = 12) + 
    theme(
      legend.title = ggplot2::element_text(size = 14, color = "black", family = "serif", face = "bold"),
      legend.text = ggplot2::element_text(size = 12, family = "serif"),
      axis.title.x = ggplot2::element_text(size = 14, color = "black", family = "serif"),
      axis.title.y = ggplot2::element_text(size = 14, color = "black", family = "serif"),
      axis.text.x = ggplot2::element_text(size = 12, color = "black", family = "serif"),
      axis.text.y = ggplot2::element_text(size = 12, color = "black", family = "serif"))
  
  return(p)
}
