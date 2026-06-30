#' Alpha diversity correlation plot
#'
#' This function generates a boxplot or barplot to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows.
#'              The first column must contain the OTUID, ASV, or species name.
#' @param facet_orientation Whether `facet_by` appears in columns ("horizontal", default) or rows ("vertical").
#' @param plot_title Character. Title for the plot. Default \code{"default"}.
#'
#' @return A ggplot object showing alpha diversity with Hill numbers.
#' @export
#' @examples
#' \dontrun{
#' alpha_hill_corrplot(
#'   table             = table,
#'   facet_orientation = "horizontal"
#' )
#' }

alpha_hill_corrplot <- function(table,
                                facet_orientation = "horizontal",
                                plot_title = "default") {
  
  # --- Preparación de datos ---
  table <- table[, !colnames(table) %in% "taxonomy"]
  table <- data.frame(t(table))
  
  # --- Calcular Hill numbers ---
  q_data <- data.frame(
    Frequency = rowSums(table),
    q0 = hillR::hill_taxa(comm = table, q = 0),
    q1 = hillR::hill_taxa(comm = table, q = 1),
    q2 = hillR::hill_taxa(comm = table, q = 2)
  )
  
  # --- Tema base común ---
  base_theme <- .mbm_theme(
    legend_position = "none",
    extra = ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(size = 8, color = "black")
    )
  )
  
  # --- Relación de aspecto ---
  aspect_ratio_theme <- if (facet_orientation == "horizontal") NULL else 0.5
  
  # --- Gráfico q0 ---
  q0_vs_depth <- ggpubr::ggscatter(
    q_data, x = "Frequency", y = "q0",
    xlab = "Sequencing depth (number of reads)",
    add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
    add.params = list(color = "#D55E00", fill = "#56B4E9"),
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n",
      size = 4.5,
      family = "serif",
      fontface = "bold",
      color = "black"
    )
  ) +
    base_theme +
    ggplot2::labs(y = expression(paste(italic("q"), "=0", " (number of total features)"))) +
    ggplot2::theme(aspect.ratio = aspect_ratio_theme)
  
  # --- Gráfico q1 ---
  q1_vs_depth <- ggpubr::ggscatter(
    q_data, x = "Frequency", y = "q1",
    xlab = "Sequencing depth (number of reads)",
    add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
    add.params = list(color = "#D55E00", fill = "#56B4E9"),
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n",
      size = 4.5,
      family = "serif",
      fontface = "bold",
      color = "black"
    )
  ) +
    base_theme +
    ggplot2::labs(y = expression(paste(italic("q"), "=1", " (number of frequent features)"))) +
    ggplot2::theme(aspect.ratio = aspect_ratio_theme)
  
  # --- Gráfico q2 ---
  q2_vs_depth <- ggpubr::ggscatter(
    q_data, x = "Frequency", y = "q2",
    xlab = "Sequencing depth (number of reads)",
    add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
    add.params = list(color = "#D55E00", fill = "#56B4E9"),
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n",
      size = 4.5,
      family = "serif",
      fontface = "bold",
      color = "black"
    )
  ) +
    base_theme +
    ggplot2::labs(y = expression(paste(italic("q"), "=2", " (number of dominant features)"))) +
    ggplot2::theme(aspect.ratio = aspect_ratio_theme)
  
  # --- Componer grid de gráficos ---
  grid_plot <- if (facet_orientation == "horizontal") {
    cowplot::plot_grid(
      q0_vs_depth, q1_vs_depth, q2_vs_depth,
      labels = c("A", "B", "C"),
      nrow = 1,
      label_fontfamily = "serif",
      label_fontface = "bold",
      label_size = 16
    )
  } else {
    cowplot::plot_grid(
      q0_vs_depth, q1_vs_depth, q2_vs_depth,
      labels = c("A", "B", "C"),
      ncol = 1,
      label_fontfamily = "serif",
      label_fontface = "bold",
      label_size = 16
    )
  }
  
  # --- Título del gráfico ---
  if (is.character(plot_title)) {
    if (tolower(plot_title) == "none") {
      return(grid_plot)
    } else if (tolower(plot_title) == "default") {
      title <- cowplot::ggdraw() +
        cowplot::draw_label("Alpha diversity vs sequencing depth",
                            fontface = "bold",
                            fontfamily = "serif")
    } else {
      title <- cowplot::ggdraw() +
        cowplot::draw_label(plot_title,
                            fontface = "bold",
                            fontfamily = "serif")
    }
    cowplot::plot_grid(title, grid_plot, ncol = 1, rel_heights = c(0.1, 1))
  } else {
    return(grid_plot)
  }
}




