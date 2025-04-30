#' Alpha diversity correlation plot
#'
#' This function generates a boxplot or barplot to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows.
#'              The first column must contain the OTUID, ASV, or species name.
#' @param facet_orientation Whether `facet_by` appears in columns ("horizontal", default) or rows ("vertical").
#' @param n_cols Number of columns in facet wrap (optional).
#' @param n_rows Number of rows in facet wrap (optional).
#' @param strip_color Background color of facet strips. Default: "grey".
#' @param free_y Logical. Wether if scales in y are free or not.

#'
#' @return A ggplot object showing alpha diversity with Hill numbers.
#' @export
#' @examples
#' library(vegan)
#' data(dune)
#' data(dune.env)
#' alpha_hill_corplot(
#'     table = t(dune),
#'     metadata = dune.env %>% tibble::rownames_to_column("SampleID"),
#'     facet_orientation = "horizontal"
#' )
alpha_hill_plot <- function(table,
                            facet_orientation = "horizontal",
                            n_cols = NULL,
                            n_rows = NULL,
                            strip_color = "grey",
                            free_y = FALSE) {
  table <- data.frame(t(table))
  
  q_data <- data.frame(
    q0 = hillR::hill_taxa(comm = table, q = 0),
    q1 = hillR::hill_taxa(comm = table, q = 1),
    q2 = hillR::hill_taxa(comm = table, q = 2),
  )
  
  q0_vs_depth_fltr <- ggscatter(
    q_data ,
    x = "Frequency",
    y = "q0",
    xlab = "Sequencing depth (number of reads)",
    ylab = "q=0",
    #ylab="Alpha diversity q=0 (effective number of total ASVs)",
    add = "reg.line",
    # Add regression line
    add.params = list(color = "#B03A2E", fill = "#566573"),
    # Customize reg. line
    conf.int = TRUE,
    # Add confidence interval
    cor.coef = TRUE,
    # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n"
    )
  ) + theme_grey() +
    theme(legend.title = element_blank(), legend.position = "none") + labs(y =
                                                                             expression(paste(
                                                                               italic("q"), "=0", " (number of total features)"
                                                                             )))
  
  q1_vs_depth_fltr <- ggscatter(
    q_data,
    x = "Frequency",
    y = "q1",
    xlab = "Sequencing depth (number of reads)",
    ylab = "q=1",
    #ylab="Alpha diversity q=1 (effective number of total ASVs)",
    add = "reg.line",
    # Add regression line
    add.params = list(color = "#B03A2E", fill = "#566573"),
    # Customize reg. line
    conf.int = TRUE,
    # Add confidence interval
    cor.coef = TRUE,
    # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n"
    )
  ) + theme_grey() +
    theme(legend.title = element_blank(), legend.position = "none") + labs(y =
                                                                             expression(paste(
                                                                               italic("q"), "=1", " (number of frequent features)"
                                                                             )))
  
  
  
  q2_vs_depth_fltr <- ggscatter(
    q_data,
    x = "Frequency",
    y = "q2",
    xlab = "Sequencing depth (number of reads)",
    ylab = "q=2",
    #ylab="Alpha diversity q=2 (effective number of total ASVs)",
    add = "reg.line",
    # Add regression line
    add.params = list(color = "#B03A2E", fill = "#566573"),
    # Customize reg. line
    conf.int = TRUE,
    # Add confidence interval
    cor.coef = TRUE,
    # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n"
    )
  ) + theme_grey() +
    theme(legend.title = element_blank(), legend.position = "none") + labs(y =
                                                                             expression(paste(
                                                                               italic("q"), "=2", " (number of dominant features)"
                                                                             )))
  
  
  library(ggplot2)
  library(cowplot)
  
  #save plots
  
  #Combine plot
  title_corr_plot <- ggdraw() + draw_label("Alpha diversity depth correlation with samples")
  correlation_plot_q012_fltr <- plot_grid(
    q0_vs_depth_fltr,
    q1_vs_depth_fltr,
    q2_vs_depth_fltr,
    labels = c("A", "B", "C"),
    nrow = 1,
    rel_heights = c(1, 1, 1)
  )

  facet_scales <- if (free_y)
    "free_y"
  else
    "fixed"
  facet_config <- if (!is.null(facet_by)) {
    # Cuando se define orientación, usamos facet_wrap con ~ q + facet_by
    formula_facet <- stats::as.formula(paste("~", if (facet_orientation == "horizontal") {
      paste("q", facet_by, sep = " + ")
    } else {
      paste(facet_by, "q", sep = " + ")
    }))
    
    ggplot2::facet_wrap(
      formula_facet,
      ncol = if (!is.null(n_cols))
        n_cols
      else
        NULL,
      nrow = if (!is.null(n_rows))
        n_rows
      else
        NULL,
      scales = if (free_y)
        "free_y"
      else
        "fixed"
    )
  } else {
    # Solo hay q
    if (facet_orientation == "horizontal") {
      ggplot2::facet_wrap(
        ~ q,
        ncol = if (!is.null(n_cols))
          n_cols
        else
          length(unique(results_largo$q)),
        scales = if (free_y)
          "free_y"
        else
          "fixed"
      )
    } else {
      ggplot2::facet_wrap(
        ~ q,
        nrow = if (!is.null(n_rows))
          n_rows
        else
          length(unique(results_largo$q)),
        scales = if (free_y)
          "free_y"
        else
          "fixed"
      )
    }
  }
  
  
  
  
  capa_geom <- if (type == "boxplot") {
    ggplot2::geom_boxplot(width = 0.5,
                          position = ggplot2::position_dodge(width = 0.9))
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
  # Si hay facet_by, incluimos esa variable en el agrupamiento
  p <- ggplot2::ggplot(results_largo, ggplot2::aes(x = .data[[x_col]], y = value, fill = .data[[fill_col]])) +
    capa_geom +
    capa_error +
    fill_scale +
    facet_config +
    ggplot2::labs(
      title = figure_title,
      x = axis_x_title,
      y = axis_y_title,
      fill = legend_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(1, "lines"),
      # aspect.ratio = 0.6,
      strip.text = ggplot2::element_text(
        face = "bold",
        color = "black",
        size = 15
      ),
      strip.background = ggplot2::element_rect(fill = strip_color),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      axis.text.x = ggplot2::element_text(size = 10, color = "black"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = if (show_legend)
        legend_position
      else
        "none"
    ) + aspect_ratio_theme
  
  if (!is.null(stat)) {
    split_vars <- if (!is.null(facet_by))
      c("q", facet_by)
    else
      "q"
    p_vals_layers <- results_largo %>%
      dplyr::group_split(across(all_of(split_vars))) %>%
      purrr::map( ~ {
        y_val <- max(.x$value, na.rm = TRUE) * 0.98
        ggpubr::stat_compare_means(
          data = .x,
          method = stat,
          label = "p.format",
          size = 3,
          hide.ns = TRUE,
          mapping = ggplot2::aes(label = ..p.format..),
          label.y = y_val
        )
      })
    
    
    p <- p + p_vals_layers
  }
  
  return(p)
}
