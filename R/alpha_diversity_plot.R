#' Alpha diversity plot
#'
#' This function generates a boxplot or barplot to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows.
#'              The first column must contain the OTUID, ASV, or species name.
#' @param metadata A data frame with metadata. The first column must match sample names in `table`.
#' @param type Type of plot: either "boxplot" or "barplot". Default is "boxplot".
#' @param stat Optional. A string indicating the test used for comparing means (e.g., "wilcox.test").
#' @param x_col Column in `metadata` to be used on the x-axis.
#' @param fill_col Column in `metadata` to define fill color.
#' @param facet_by Optional. A metadata column to facet (e.g., Treatment, Site).
#' @param facet_by2 Optional. A metadata column to double facet (e.g., Treatment, Site).
#' @param facet_orientation Whether `facet_by` appears in columns ("horizontal", default) or rows ("vertical").
#' @param palette Color palette to use: "colorb", "grey", "viridis", or "brewer". Default: "colorb".
#' @param group_colors A vector of custom colors. Overrides `palette` if provided.
#' @param n_cols Number of columns in facet wrap (optional).
#' @param n_rows Number of rows in facet wrap (optional).
#' @param strip_color Background color of facet strips. Default: "grey".
#' @param show_legend Logical. Show legend? Default: TRUE.
#' @param legend_title Title for the legend.
#' @param legend_position Position of the legend: "bottom", "top", "right", or "left". Default is "bottom".
#' @param title Title for the entire plot.
#' @param x_axis_title Title for the x-axis.
#' @param y_axis_title Title for the y-axis.
#' @param free_y Logical. Whether y-axis scales are free across facets. Default \code{FALSE}.
#' @param rarefy_depth Integer. If provided, rarefies samples to this depth before diversity estimation. Default \code{NULL}.
#' @param panel_label_case Character. Case of the auto-generated panel tags
#'   (A, B, C... added per panel when \code{stat} is used). One of
#'   \code{"upper"} (default, "A", "B", "C") or \code{"lower"} ("a", "b", "c").
#'   Ignored if \code{panel_labels} is supplied.
#' @param panel_labels Optional character vector of custom panel tags, one per
#'   panel, used as-is (e.g. \code{c("(a)", "(b)", "(c)")} or
#'   \code{c("a.", "b.", "c.")}) — for journal styles that
#'   \code{panel_label_case} alone can't produce. Overrides
#'   \code{panel_label_case} when provided.
#' @param panel_label_bold Logical. If \code{TRUE} (default), panel tags are
#'   bold. Set to \code{FALSE} for journals that require plain (non-bold)
#'   panel tags.
#' @param save_table Logical. If \code{TRUE}, saves the diversity table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table. Default \code{"diversity.txt"}.
#'
#' @return A ggplot object showing alpha diversity with Hill numbers.
#' @export
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "table_with_taxonomy.tsv", package = "MicroBioMeta")
#' table <- read.delim(table_path, skip = 1, comment.char = "", check.names = FALSE, row.names = 1)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE, comment.char = "")
#' colnames(metadata)[1] <- "SAMPLEID"
#'
#' alpha_diversity_plot(
#'   table           = table,
#'   metadata        = metadata,
#'   type            = "barplot",
#'   x_col           = "Type_of_soil",
#'   fill_col        = "Type_of_soil",
#'   facet_by2       = "Treatment",
#'   free_y          = TRUE,
#'   legend_position = "top",
#'   stat            = "t.test"
#' )
#' }

alpha_diversity_plot <- function(
    table,
    metadata,
    type = "boxplot",
    stat = NULL,
    x_col,
    fill_col,
    facet_by = NULL,
    facet_by2 = NULL,  # Nuevo parámetro
    facet_orientation = "horizontal",
    palette = "colorb",
    group_colors = NULL,
    n_cols = NULL,
    n_rows = NULL,
    strip_color = "grey",
    show_legend = TRUE,
    title = NULL,
    legend_title = NULL,
    legend_position = "bottom",
    x_axis_title = NULL,
    y_axis_title = "Diversity measure",
    free_y = FALSE,
    rarefy_depth=NULL,
    panel_label_case = "upper",
    panel_labels = NULL,
    panel_label_bold = TRUE,
    save_table = FALSE,
    table_filename = "diversity.txt") {

  # Treat metadata's first column as the sample ID regardless of its original name
  colnames(metadata)[1] <- "SAMPLEID"

  #considero que esta parte no debería ser necesaria, no sé si es buena idea remover asvs
  #que no tengan una identificación taxonómica no quiere decir que no se tengan que considerar
  #table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
  # Remove uninformative taxonomy strings
  #table <- table %>%
   # dplyr::filter(taxonomy != "d__Bacteria;__;__;__;__;__")
  
  # Keep only samples present in both table and metadata. Final row order
  # doesn't matter here: per-sample diversity indices are computed
  # independently, and results are joined back to metadata by SAMPLEID below.
  sample_order <- metadata[[1]]
  common_samples <- intersect(colnames(table), sample_order)
  if (length(common_samples) == 0) stop("No matching sample names between table and metadata.")
  table <- table[, common_samples, drop = FALSE]
  table <- data.frame(t(table))
  
  if (!is.null(rarefy_depth) && rarefy_depth > 0) {
    table <- vegan::rrarefy(table, sample = rarefy_depth)
  }
  
  # Estimar índices
  est_richness <- vegan::estimateR(table)  # muestras como filas
  chao1 <- est_richness["S.chao1", ]
  shannon <- vegan::diversity(table, index = "shannon")
  simpson <- vegan::diversity(table, index = "simpson")
  
  # Crear data frame unificado
  results <- data.frame(
    SAMPLEID = rownames(table),
    Chao1 = as.numeric(chao1),
    Shannon = shannon,
    Simpson = simpson
  )
  
  # Unir con metadata por SAMPLEID
  results <- dplyr::left_join(results, metadata, by = "SAMPLEID")
  
  # Guardar tabla si se solicita
  if (save_table) {
    utils::write.table(
      results,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message(paste("Table saved as:", table_filename))
  }
  
  results[[fill_col]] <- factor(results[[fill_col]], levels = unique(results[[fill_col]]))
  results[[x_col]] <- factor(results[[x_col]], levels = unique(results[[x_col]]))
  
  results_largo <- tidyr::pivot_longer(
    results,
    cols = c("Chao1", "Shannon", "Simpson"),
    names_to = "Index",
    values_to = "value"
  )
  
  fill_scale <- if (!is.null(group_colors)) {
    ggplot2::scale_fill_manual(values = group_colors)
  } else {
    switch(palette,
           "colorb" = ggplot2::scale_fill_manual(values = .mbm_colors),
           "grey" = ggplot2::scale_fill_grey(start = 0.9, end = 0.3),
           "viridis" = ggplot2::scale_fill_viridis_d(option = "plasma"),
           "brewer" = ggplot2::scale_fill_brewer(palette = "Set2"),
           ggplot2::scale_fill_grey(start = 0.9, end = 0.3)
    )
  }
  
  facet_scales <- if (free_y) "free_y" else "fixed"  # Esto ahora se usa correctamente
  
  facet_config <- if (!is.null(facet_by) && !is.null(facet_by2)) {
    formula_facet <- if (facet_orientation == "horizontal") {
      stats::as.formula(paste(facet_by2, "~",  "Index +", facet_by))
    } else {
      stats::as.formula(paste("Index +", facet_by, "~", facet_by2))
    }
    ggh4x::facet_nested(
      formula_facet,
      nest_line = ggplot2::element_line(colour = "black"),
      scales = facet_scales
    )
  } else if (!is.null(facet_by)) {
    formula_facet <- if (facet_orientation == "horizontal") {
      stats::as.formula(paste(facet_by, "~ Index"))
    } else {
      stats::as.formula(paste("Index ~", facet_by))
    }
    
    if (free_y) {
      ggh4x::facet_grid2(
        formula_facet,
        scales = "free_y",
        independent = "y"
      )
    } else {
      ggh4x::facet_grid2(
        formula_facet,
        scales = "fixed"
      )
    }
  } else {
    ggplot2::facet_wrap(
      ~Index,
      ncol = if (facet_orientation == "horizontal") {
        length(unique(results_largo$Index))
      } else {
        NULL
      },
      nrow = if (facet_orientation == "vertical") {
        length(unique(results_largo$Index))
      } else {
        NULL
      },
      scales = if (free_y) "free_y" else "fixed"
    )
  }
  
  
  capa_geom <- if (type == "boxplot") {
    ggplot2::geom_boxplot(
      width = 0.5,
      position = ggplot2::position_dodge(width = 0.9)
    )
  } else if (type == "barplot") {
    ggplot2::geom_bar(
      stat = "summary",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.5,
      color = "black"
    )
  } else {
    stop("Invalid plot type. Use 'boxplot' or 'barplot'.")
  }
  
  capa_error <- if (type == "barplot") {
    ggplot2::stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.9)
    )
  } else {
    NULL
  }
  
  aspect_ratio_theme <- if (facet_orientation == "horizontal") {
    ggplot2::theme(aspect.ratio = 1.5)
  } else {
    ggplot2::theme(aspect.ratio = 0.7)
  }
  
  # Calcular label.y ajustado por q
  p <- ggplot2::ggplot(
    results_largo,
    ggplot2::aes(x = .data[[x_col]], y = value, fill = .data[[fill_col]])
  ) +
    capa_geom +
    capa_error +
    fill_scale +
    facet_config +
    ggplot2::labs(
      title = title,
      x = x_axis_title,
      y = y_axis_title,
      fill = legend_title
    ) +
    .mbm_theme(
      legend_position = if (show_legend) legend_position else "none",
      extra = ggplot2::theme(
        panel.grid     = ggplot2::element_blank(),
        panel.spacing  = grid::unit(1, "lines"),
        strip.background = ggplot2::element_rect(fill = strip_color, color = "black")
      )
    ) + aspect_ratio_theme 
  
  if (!is.null(stat)) {
    split_vars <- if (!is.null(facet_by)) c("Index", facet_by) else "Index"
    p_vals_layers <- results_largo %>%
      dplyr::group_split(dplyr::across(dplyr::all_of(split_vars))) %>%
      purrr::map(~ {
        y_val <- max(.x$value, na.rm = TRUE) * 0.98
        
        # Obtener los niveles únicos del eje x y calcular el valor central
        x_levels <- levels(factor(.x[[x_col]]))
        x_numeric <- match(.x[[x_col]], x_levels)
        x_center <- mean(range(x_numeric, na.rm = TRUE))
        
        ggpubr::stat_compare_means(
          data = .x,
          method = stat,
          mapping = ggplot2::aes(
            label = scales::label_pvalue(
              accuracy = 0.001, prefix = c("p < ", "p = ", "p > ")
            )(ggplot2::after_stat(p))
          ),
          size = 3.5,
          family= "serif",
          hide.ns = TRUE,
          label.y = y_val,
          label.x = x_center
        )
      })
    
    p <- p + p_vals_layers
    # Añadir etiquetas tipo A, B, C... a los paneles
    panel_letters <- if (!is.null(panel_labels)) {
      panel_labels
    } else if (identical(panel_label_case, "lower")) letters else LETTERS
    gb <- ggplot2::ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = panel_letters[lay$PANEL], x = -Inf, y = Inf)

    p <- p + ggplot2::geom_text(
      data = tags,
      mapping = ggplot2::aes(x = x, y = y, label = label),
      hjust = -0.5,
      vjust = 1.5,
      fontface = if (panel_label_bold) "bold" else "plain",
      family= "serif",
      size= 6,
      inherit.aes = FALSE
    )
  }
  
  return(p)
}
