#' This function generates either an effect size plot or a volcano plot based on ALDEx2 results.
#'
#' @param table Data frame with count data; columns represent samples, rows represent features.
#' @param metadata Data frame containing metadata for the samples.
#' @param col_cond Name of the column in `metadata` that contains the experimental conditions.
#' @param type Type of plot to generate: "effect" for effect size plot or "volcano" for volcano plot.
#' @param col_inf Color for points with effect size/difference lower than threshold (default for "effect" plot).
#' @param col_sup Color for points with effect size/difference higher than threshold (default for "effect" plot).
#' @param threshold_lower Lower threshold for effect size/difference (x-axis).
#' @param threshold_upper Upper threshold for effect size/difference (x-axis).
#' @param cond Name of the condition that appears first in `table` (used in plot labels).
#' @param cutoff.pval p-value cutoff for significance (default = 0.05).
#' @param show_labels Logical. Whether to display "Higher/Lower in cond" labels (for "effect" plot only, default is TRUE).
#' @param taxa Data frame with taxonomic information (required for "volcano" plot only).
#'
#' @return A `ggplot` object with the selected plot.
#' @export
#'
aldex_volcano_plot <- function(table,
                               metadata,
                               col_cond,
                               type = "volcano",
                               col_inf = "blue",
                               col_sup = "red",
                               threshold_lower = -1.5,
                               threshold_upper = 1.5,
                               cond = NULL,
                               cutoff.pval = 0.05,
                               show_labels = TRUE,
                               taxa = NULL) {
  
  # Verificar que type tiene un valor válido
  if (!type %in% c("effect", "volcano")) {
    stop("type must be either 'effect' or 'volcano'")
  }
  
  # Verificar paquetes requeridos
  # Cargar ggtext explícitamente
  if (!require(ggtext, quietly = TRUE)) {
    stop("Package 'ggtext' required for formatted text. Please install it.")
  }
  
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("Package 'ALDEx2' needed for this function to work. Please install it.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.")
  }
  
  # Verificar que la columna de condición existe
  if (!col_cond %in% colnames(metadata)) {
    stop(
      paste(
        "The column",
        col_cond,
        "does not exist in the object 'metadata'. Check the name is written correctly."
      )
    )
  }
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  # Identificar automáticamente la columna taxonómica (última columna)
  if (is.null(taxa)) {
    last_col <- ncol(table)
    taxa_colname <- colnames(table)[last_col]
    taxa <- data.frame(
      Feature.ID = rownames(table),
      Taxon = table[[taxa_colname]],
      stringsAsFactors = FALSE
    )
    # Eliminar la última columna de table para el análisis
    table <- table[, -last_col, drop = FALSE]
  }
  
  conditions <- metadata[[col_cond]]
  groups <- unique(conditions)
  if (is.null(cond)) cond <- groups[1]
  other_cond <- setdiff(groups, cond)[1]
  
  aldex_clr <- ALDEx2::aldex(table, conditions, mc.samples = 128, denom = "all")
  
  # Procesar datos taxonómicos para ambos tipos de gráficos
  processed_data <- aldex_clr %>%
    tibble::rownames_to_column(var = "Feature.ID") %>%
    dplyr::left_join(taxa, by = "Feature.ID") %>%
    dplyr::mutate(
      taxa = dplyr::case_when(
        stringr::str_detect(Taxon, "g__") ~ stringr::str_extract(Taxon, "(?<=g__)[^_;]+"),
        stringr::str_detect(Taxon, "f__") ~ stringr::str_extract(Taxon, "(?<=f__)[^_;]+"),
        stringr::str_detect(Taxon, "c__") ~ stringr::str_extract(Taxon, "(?<=c__)[^_;]+"),
        stringr::str_detect(Taxon, "o__") ~ stringr::str_extract(Taxon, "(?<=o__)[^_;]+"),
        TRUE ~ Feature.ID
      )
    )
  
  if (type == "effect") {
    # Preparar datos para effect plot
    plot_data <- processed_data %>%
      dplyr::mutate(
        grupo = dplyr::case_when(
          effect <= threshold_lower ~ paste("Lower in", cond),
          effect >= threshold_upper ~ paste("Higher in", cond),
          TRUE ~ "Not significant"
        ),
        log_pvalue = -log10(wi.ep + min(wi.ep[wi.ep > 0])/10)
      )
    
    # Obtener los top taxones para etiquetar
    top_taxa <- plot_data %>%
      dplyr::filter(grupo != "Not significant") %>%
      dplyr::group_by(grupo) %>%
      dplyr::slice_max(order_by = abs(effect), n = 3)
    
    lim_x <- max(abs(plot_data$effect), na.rm = TRUE)
    
    # Crear vector de colores con nombres dinámicos
    color_values <- c(col_sup, col_inf, "gray")
    names(color_values) <- c(paste("Higher in", cond), paste("Lower in", cond), "Not significant")
    
    # Crear el gráfico base
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect, y = log_pvalue, color = grupo)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::geom_vline(
        xintercept = c(threshold_lower, threshold_upper),
        linetype = 2,
        color = "black"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(cutoff.pval),
        linetype = 2,
        color = "black"
      ) +
      ggplot2::labs(
        x = "Effect size",
        y = expression("-Log"[10]~"p-value"),
        color = NULL
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 12, color = "black"),
        axis.title = ggplot2::element_text(size = 12),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(limits = c(-lim_x, lim_x))
    
    # Añadir etiquetas de taxones significativos
    if (nrow(top_taxa) > 0) {
      p <- p +
        ggplot2::geom_text(
          data = top_taxa,
          ggplot2::aes(label = taxa),
          color = "black",
          size = 3,
          vjust = -0.5,
          fontface = "italic"
        )
    }
    
    # Añadir etiquetas de condición con formato richtext
    if (show_labels) {
      p <- p +
        ggplot2::annotate(
          "richtext",
          x = threshold_lower,
          y = max(plot_data$log_pvalue) * 0.95,
          label = paste0("<b style='color:", col_inf, "'>Lower in ", cond, "</b>"),
          size = 5,
          hjust = 1,
          vjust = 1
        ) +
        ggplot2::annotate(
          "richtext",
          x = threshold_upper,
          y = max(plot_data$log_pvalue) * 0.95,
          label = paste0("<b style='color:", col_sup, "'>Higher in ", cond, "</b>"),
          size = 5,
          hjust = 0,
          vjust = 1
        )
    }
    
  } else if (type == "volcano") {
    # Preparar datos para volcano plot
    plot_data <- processed_data %>%
      dplyr::mutate(
        log_pvalue = -log10(wi.ep + min(wi.ep[wi.ep > 0])/10),
        significant = wi.ep <= cutoff.pval,
        direction = ifelse(diff.btw < 0, 
                           paste("Lower in", cond),
                           paste("Higher in", cond))
      )
    
    # Obtener los top taxones para etiquetar
    top_taxa <- plot_data %>%
      dplyr::filter(significant) %>%
      dplyr::group_by(direction) %>%
      dplyr::slice_max(order_by = abs(diff.btw), n = 3)
    
    # Crear vector de colores con nombres dinámicos
    color_values <- c(col_sup, col_inf)
    names(color_values) <- c(paste("Higher in", cond), paste("Lower in", cond))
    
    # Crear el gráfico base
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = diff.btw, y = log_pvalue)) +
      ggplot2::geom_point(
        data = dplyr::filter(plot_data, !significant),
        color = "gray",
        size = 3,
        alpha = 0.7
      ) +
      ggplot2::geom_point(
        data = dplyr::filter(plot_data, significant),
        ggplot2::aes(color = direction),
        size = 3
      ) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::geom_vline(
        xintercept = c(threshold_lower, threshold_upper),
        color = 'black',
        linetype = 'dashed'
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(cutoff.pval),
        color = 'black',
        linetype = 'dashed'
      ) +
      ggplot2::labs(
        x = expression("Log"[2]~"Fold Change"),
        y = expression("-Log"[10]~"p-value"),
        color = NULL
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 12, color = "black"),
        axis.title = ggplot2::element_text(size = 12),
        legend.position = "none"
      )
    
    # Añadir etiquetas de taxones
    if (nrow(top_taxa) > 0) {
      p <- p +
        ggplot2::geom_text(
          data = top_taxa,
          ggplot2::aes(label = taxa),
          color = "black",
          size = 3,
          vjust = -0.5,
          fontface = "italic"
        )
    }
    
    # Añadir etiquetas de condición con formato richtext
    p <- p +
      ggplot2::annotate(
        "richtext",
        x = min(plot_data$diff.btw) + 2,
        y = max(plot_data$log_pvalue) * 0.95,
        label = paste0("<b style='color:", col_inf, "'>Lower in ", cond, "</b>"),
        color = col_inf,
        size = 5,
        hjust = 1,
        vjust = 1
      ) +
      ggplot2::annotate(
        "richtext",
        x = max(plot_data$diff.btw) - 2,
        y = max(plot_data$log_pvalue) * 0.95,
        label = paste0("<b style='color:", col_sup, "'>Higher in ", cond, "</b>"),
        color = col_sup,
        size = 5,
        hjust = 0,
        vjust = 1
      ) 
  }
  
  return(p)
}