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
  
  # --- Filtrar metadatos ---
  metadata_sub <- metadata %>%
    filter(.data[[condition_col]] %in% c(condition_A, condition_B)) %>%
    select(SampleID = 1, Condition = all_of(condition_col))
  
  samples <- metadata_sub$SampleID
  
  # --- Ordenar columnas y empatar con metadata ---
  table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
  ordered_samples <- intersect(metadata[, 1], colnames(table)[-1])
  table <- table[, c("taxonomy", ordered_samples)]
  
  # --- Asegurar identificadores únicos ---
  if (!is.null(rownames(table))) {
    table <- tibble::rownames_to_column(table, var = "OTU_ID")
  } else {
    table$OTU_ID <- paste0("OTU_", seq_len(nrow(table)))
  }
  
  # --- Limpiar terminaciones vacías en la taxonomía ---
  table$taxonomy <- gsub("(;__)+$", "", table$taxonomy)
  
  # --- Determinar índice del nivel taxonómico ---
  level_idx <- switch(level,
                      kingdom = 1, phylum = 2, class = 3, order = 4,
                      family = 5, genus = 6, species = 7)
  
  # --- Calcular profundidad taxonómica de cada OTU ---
  get_depth <- function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    sum(grepl("__", levels))
  }
  table$depth <- sapply(table$taxonomy, get_depth)
  
  # --- Separar filas según resolución taxonómica ---
  lowres <- table[table$depth < level_idx, ]    
  highres <- table[table$depth >= level_idx, ]  
  
  # --- Recortar y colapsar taxonomías con suficiente resolución ---
  highres$taxonomy <- sapply(highres$taxonomy, function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    paste(levels[1:level_idx], collapse = ";")
  })
  
  # --- Colapsar correctamente taxones repetidos ---
  highres <- highres %>%
    group_by(taxonomy) %>%
    summarise(
      OTU_ID = paste(unique(OTU_ID), collapse = ";"),
      across(where(is.numeric), sum, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # --- Combinar lowres y highres ---
  table_final <- dplyr::bind_rows(
    lowres[, c("OTU_ID", "taxonomy", ordered_samples)],
    highres[, c("OTU_ID", "taxonomy", ordered_samples)]
  )
  
  # --- Limpiar tabla final ---
  table_final <- table_final[, c("OTU_ID",  ordered_samples, "taxonomy")]
  table_final <- tibble::column_to_rownames(table_final, "OTU_ID")
  
  # --- Corregir taxonomía ---
  abundance_raw <- table_final %>%
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
    )
  
  # --- Abundancia relativa ---
  abundance_rel <- abundance_raw %>%
    mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE) * 100))
  
  # --- Formato largo ---
  long_data <- abundance_rel %>%
    rownames_to_column("OTU_ID") %>%
    pivot_longer(cols = all_of(ordered_samples), names_to = "SampleID", values_to = "Abundance") %>%
    left_join(metadata_sub, by = "SampleID")
  
  # --- Cálculo de medias y ratios (manteniendo ASVs) ---
  summary_data <- long_data %>%
    group_by(OTU_ID, taxonomy, Condition) %>%
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
  
  # --- Crear etiqueta combinada para ASVs repetidos ---
  # Crear etiquetas cortas ASV1, ASV2, ...
  summary_data <- summary_data %>%
    mutate(ASV_label = paste0("ASV", row_number()))
  
  #summary_data <- summary_data %>%
   # mutate(
    #  taxonomy_display = if_else(
     #   duplicated(taxonomy) | duplicated(taxonomy, fromLast = TRUE),
      #  paste0(taxonomy, " (", OTU_ID, ")"),
       # taxonomy
      #)
    #)
  
  summary_data <- summary_data %>%
    mutate(
      taxonomy_display = if_else(
        duplicated(taxonomy) | duplicated(taxonomy, fromLast = TRUE),
        paste0(taxonomy, " (", ASV_label, ")"),
        taxonomy
      )
    )
  
  
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
                y = reorder(taxonomy_display, MeanAbund),
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
