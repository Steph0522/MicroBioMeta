#' Plot Correlation Between Environmental Variables and Taxonomic Groups
#'
#' This function calculates the relative abundances of taxa at a specified taxonomic level 
#' (phylum, genus, or species) from a count table, computes correlations between these 
#' abundances and environmental variables, and visualizes the results as either a heatmap 
#' (tile) or a bubble plot (circle).
#'
#' @param table A data frame containing count data with a column named `taxonomy`. Each 
#'   row corresponds to a taxon, and each column (besides `taxonomy`) corresponds to a sample.
#' @param env_table A data frame of environmental variables, with samples as row names.
#' @param method Correlation method to use. Options include `"spearman"`, `"pearson"`, or `"kendall"`. Default is `"spearman"`.
#' @param cond_vect Vector with environmental variables to consider in the correlation
#' @param hc.order Logical. If `TRUE`, performs hierarchical clustering to reorder rows and columns based on correlation similarity. Default is `TRUE`.
#' @param geom Character. Type of plot to generate: `"tile"` for a heatmap or `"circle"` for a bubble plot. Default is `"tile"`.
#' @param show_labels Logical. Whether to display correlation values on the plot. Default is `TRUE`.
#' @param col_palette A vector of colors for the gradient scale. If `NULL`, a default blue-white-red palette is used.
#' @param invert_axes Logical. If `TRUE`, environmental variables are shown on the x-axis and taxonomic groups on the y-axis. Default is `TRUE`.
#' @param taxonomy_db Character. Database used for taxonomy annotation. Options are `"silva"` or `"Kraken2"`. Default is `"silva"`.
#' @param level Character. Taxonomic level for collapsing counts: `"phylum"`, `"genus"`, or `"specie"`. Default is `"genus"`.
#'
#' @ret
#' @export
#'
#' @examples
#' 
#' 

#colores<- c("pink","white","purple")
#corr_env_abund_plot(table = abund, 
#                    env_table = env,
#                   cond_vect= c("ph","OM")
#                    method = "pearson", 
#                    geom = "tile", 
#                    hc.order = FALSE, 
#                    col_palette = colores,
#                    invert_axes = TRUE,
#                    show_labels = FALSE,
#                    level = "genus")


corr_env_abund_plot <- function(table,
                                env_table,
                                cond_vect= NULL,
                                method = "spearman",
                                hc.order = TRUE,
                                geom = c("tile", "circle"),
                                show_labels = TRUE,
                                col_palette = NULL,
                                invert_axes = TRUE,
                                taxonomy_db = "silva",
                                level = "genus") {
  rownames(table) <- NULL
  
  #colapsar la tabla al nivel taxonómico deseado: filo, género o especie
  if (level == "genus") {
    table$taxonomy <- gsub(";\\s?s__.*", "", table$taxonomy)
  }
  if (level == "phylum") {
    table$taxonomy <- gsub(";\\s?c__.*", "", table$taxonomy)
  }
  if (level == "specie") {
    table$taxonomy <- table$taxonomy
  }
  
  table <- table %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE))
  
  #modificar la columna taxonomy para solo conservar el nombre al nivel que colapsamos
  #esto hace que al graficar salga sólo ese nombre y no toda la taxonomía
  
  if (taxonomy_db %in% c("silva", "Kraken2") && level == "phylum") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("p__[^;]*", taxonomy) ~ stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  if (taxonomy_db %in% c("silva", "Kraken2") && level == "genus") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
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
  
  # Filas comunes
  common_samples <- base::intersect(colnames(table), rownames(env_table))
  counts <- table[, common_samples, drop = FALSE]
  env <- env_table[common_samples, , drop = FALSE]
  
  # Seleccionar solo las variables ambientales indicadas en el vector cond_vect
  if (!is.null(cond_vect)) {
    env <- env[, cond_vect, drop = FALSE]
  }
  
  
  # Convertir a abundancias relativas (porcentaje)
  abund <-
    sweep(counts, 2, colSums(counts, na.rm = TRUE), FUN = "/") * 100
  
  # Matriz de correlación
  corr_mat <-
    stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
  
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
