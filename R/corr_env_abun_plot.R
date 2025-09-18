#' Plot Correlation Between Environmental Variables and Taxonomic Groups
#'
#' This function calculates the relative abundances of taxa at a specified taxonomic level 
#' (phylum, genus, or species) from a count table, computes correlations between these 
#' abundances and environmental variables, and visualizes the results as either a heatmap 
#' (tile) or a bubble plot (circle).
#'
#
#'
#' @param table 
#' @param env_table 
#' @param metadata 
#' @param cond_vect 
#' @param method 
#' @param hc.order 
#' @param geom 
#' @param show_labels 
#' @param col_palette 
#' @param invert_axes 
#' @param taxonomy_db 
#' @param level 
#' @param pval_threshold 
#'
#' @return
#' @export
#'
#' @examples
#'   colores<- c("pink","white","purple")
#    corr_env_abund_plot(table = table, 
#     env_table = env_data,
#     metadata=metadata,
#     cond_vect= c("pH","OM", "NO3","NH4"),
#     method = "pearson", 
#     geom = "tile", 
#     hc.order = FALSE, 
#     col_palette = colores,
#     invert_axes = TRUE,
#     show_labels = FALSE,
#     level = "species",
#     taxonomy_db = "unite",
#     pval_threshold= 0.05)
#
corr_env_abund_plot <- function(table,
                                env_table,
                                metadata,
                                cond_vect= NULL,
                                method = "spearman",
                                hc.order = TRUE,
                                geom = c("tile", "circle"),
                                show_labels = TRUE,
                                col_palette = NULL,
                                invert_axes = TRUE,
                                taxonomy_db = "silva",
                                level = "genus",
                                pval_threshold= NULL) {
  rownames(table) <- NULL
  
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  # Remove uninformative taxonomy strings
  table <- table %>%
    dplyr::filter(taxonomy != "d__Bacteria;__;__;__;__;__") %>%
    dplyr::filter(taxonomy != "d__Bacteria") %>%
    dplyr::filter(taxonomy != "d__Archaea;__;__;__;__;__") %>%
    dplyr::filter(taxonomy != "d__Archaea") %>%
    dplyr::filter(taxonomy != "d__Bacteria;p__;c__;o__;f__;g__;s__") %>%
    dplyr::filter(taxonomy != "d__Archaea;p__;c__;o__;f__;g__;s__") %>%
    dplyr::filter(taxonomy != "k__Bacteria;__;__;__;__;__")%>%
    dplyr::filter(taxonomy != "k__Fungi;__;__;__;__;__")%>%
    dplyr::filter(taxonomy != "k__Fungi;p__;c__;o__;f__;g__")%>%
    dplyr::filter(taxonomy != "k__Fungi")%>%
    dplyr::filter(taxonomy != "Unassigned")%>%
    dplyr::filter(taxonomy != "d__Eukaryota")
  
  
  #muestras entre table, env_table y metadata
  common_samples <- Reduce(intersect, list(colnames(table), rownames(env_table), metadata$SAMPLEID))
  table <- table[, c(tax_col, match(common_samples, colnames(table))), drop = FALSE]
  env_table <- env_table[common_samples, , drop = FALSE]
  metadata <- metadata[metadata$SAMPLEID %in% common_samples, , drop = FALSE]
  rownames(metadata) <- metadata$SAMPLEID
  
  
  #colapsar la tabla al nivel taxonómico deseado
  if (taxonomy_db %in% c("silva", "Kraken2", "gg2")) {
    if (level == "kingdom") {
      table$taxonomy <- sub(";.*", "", table$taxonomy)
    }
    if (level == "phylum") {
      table$taxonomy <- sub(";\\s?c__.*", "", table$taxonomy)
    }
    if (level == "class") {
      table$taxonomy <- sub(";\\s?o__.*", "", table$taxonomy)
    }
    if (level == "order") {
      table$taxonomy <- sub(";\\s?f__.*", "", table$taxonomy)
    }
    if (level == "family") {
      table$taxonomy <- sub(";\\s?g__.*", "", table$taxonomy)
    }
    if (level == "genus") {
      table$taxonomy <- sub(";\\s?s__.*", "", table$taxonomy)
    }
    if (level == "species") {
      table$taxonomy <- table$taxonomy
    }
  }
  
  # ---- Colapsar taxonomía según nivel ----
  if (taxonomy_db == "unite") {
    if (level == "kingdom") {
      table$taxonomy <- sub(";.*", "", table$taxonomy)
    }
    if (level == "phylum") {
      table$taxonomy <- sub(";\\s?c__.*", "", table$taxonomy)
    }
    if (level == "class") {
      table$taxonomy <- sub(";\\s?o__.*", "", table$taxonomy)
    }
    if (level == "order") {
      table$taxonomy <- sub(";\\s?f__.*", "", table$taxonomy)
    }
    if (level == "family") {
      table$taxonomy <- sub(";\\s?g__.*", "", table$taxonomy)
    }
    if (level == "genus") {
      table$taxonomy <- sub(";\\s?s__.*", "", table$taxonomy)
    }
    if (level == "species") {
      table$taxonomy <- sub(";\\s?sh__.*", "", table$taxonomy)
    }
  }
  
  
  table <- table %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE))
  
  #modificar la columna taxonomy para solo conservar el nombre al nivel que colapsamos
  #esto hace que al graficar salga sólo ese nombre y no toda la taxonomía
  
  if (taxonomy_db %in% c("unite","silva", "gg2") && level == "species") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*;.*s__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$|s__uncultured|s__$", taxonomy) ~
            paste0(
              stringr::str_extract(taxonomy, "s__[^;]*") %>% sub("s__", "", .)
            ),
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .)),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  if (taxonomy_db == "Kraken2" && level == "species") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*;.*s__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$|s__uncultured|s__$", taxonomy) ~
            paste0(
              stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .), " ",
              stringr::str_extract(taxonomy, "s__[^;]*") %>% sub("s__", "", .)
            ),
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .)),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  
  # Simplify taxonomy for SILVA
  if (taxonomy_db %in% c("silva") && level == "genus") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$|g__Incertae_Sedis", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$|f__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$|o__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  
  if (taxonomy_db %in% c("unite", "Kraken2", "gg2") && level == "genus") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "family") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$|f__Incertae_Sedis", taxonomy) ~ sub(".*f__([^;]*).*", "\\1", taxonomy),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$|o__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "order") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$|o__Incertae_Sedis", taxonomy) ~ sub(".*o__([^;]*).*", "\\1", taxonomy),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "class") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ sub(".*c__([^;]*).*", "\\1", taxonomy),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "phylum") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("p__[^;]*", taxonomy) ~ stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  # Asegurar que la columna "taxonomy" sea rownames
  if ("taxonomy" %in% colnames(table)) {
    table <- tibble::column_to_rownames(table, "taxonomy")
  }
  
  # Paleta por defecto
  if (is.null(col_palette)) {
    col_palette <-
      grDevices::colorRampPalette(c("blue", "white", "red"))(200)
  }
  # --- Filas comunes
  common_samples <- base::intersect(colnames(table), rownames(env_table))
  counts <- table[, common_samples, drop = FALSE]
  env <- env_table[common_samples, , drop = FALSE]
  

  # Seleccionar solo las variables ambientales indicadas en cond_vect, verificando coincidencias
  if (!is.null(cond_vect)) {
    cond_vect <- cond_vect[cond_vect %in% colnames(env)]
    if(length(cond_vect) == 0) stop("No matching variables found in env_table")
    env <- env[, cond_vect, drop = FALSE]
  }
  
  
  # Filtrar variables constantes
  env <- env[, apply(env, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  abund <- sweep(counts, 2, colSums(counts, na.rm = TRUE), FUN = "/") * 100
  abund <- abund[apply(abund, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
  
  # Matriz de correlación general
  corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
  
  # --- Calcular p-values si se indica pval_threshold
  if (!is.null(pval_threshold)) {
    pval_mat <- matrix(NA, 
                       nrow = ncol(env), 
                       ncol = nrow(abund),
                       dimnames = list(colnames(env), rownames(abund)))
    
    for (env_var in colnames(env)) {
      for (taxon in rownames(abund)) {
        test <- suppressWarnings(
          cor.test(env[[env_var]], as.numeric(abund[taxon, ]), method = method)
        )
        pval_mat[env_var, taxon] <- test$p.value
      }
    }
    
    # Mantener solo taxones significativos en al menos una variable
    signif_taxa <- rownames(abund)[apply(pval_mat, 2, function(x) any(x < pval_threshold, na.rm = TRUE))]
    abund <- abund[signif_taxa, , drop = FALSE]
    corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
  }
  
  
  
  # Clustering jerárquico
  if (hc.order) {
    d_row <- stats::dist(1 - corr_mat)
    hc_row <- stats::hclust(d_row)
    row_ord <- hc_row$labels[hc_row$order]
    d_col <- stats::dist(1 - t(corr_mat))
    hc_col <- stats::hclust(d_col)
    col_ord <- hc_col$labels[hc_col$order]
    corr_mat <- corr_mat[row_ord, col_ord]
  }
  
  # Convertir a formato largo
  corr_df <- reshape2::melt(corr_mat,
                            varnames = c("Environmental", "Group"),
                            value.name = "Correlation")
  
  # Base del gráfico
  if (invert_axes) {
    p <- ggplot2::ggplot(corr_df,
                         ggplot2::aes(x = Environmental, y = Group, fill = Correlation))
  } else {
    p <- ggplot2::ggplot(corr_df,
                         ggplot2::aes(x = Group, y = Environmental, fill = Correlation))
  }
  
  # Elegir tipo de gráfico, "tile" es como heatmap y "circle" como bubbleplot
  if (geom == "tile") {
    p <- p + ggplot2::geom_tile(color = "gray80")
    if (show_labels) {
      p <-
        p + ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                               size = 3,
                               color = "black")
    }
  } else if (geom == "circle") {
    p <- p + ggplot2::geom_point(ggplot2::aes(size = abs(Correlation)),
                                 shape = 21,
                                 color = "gray") +
      ggplot2::scale_size(range = c(2, 10))
    if (show_labels) {
      p <-
        p + ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                               size = 3,
                               vjust = 0.5)
    }
  }
  
  # Colores y tema
  p <- p + ggplot2::scale_fill_gradientn(colours = col_palette,
                                         limits = c(-1, 1),
                                         name = "Correlation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 10)
    ) +
    ggplot2::coord_fixed()
  
  return(p)
}
