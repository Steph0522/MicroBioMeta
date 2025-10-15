#' Alpha diversity correlation plot
#'
#' This function generates a boxplot or barplot to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows.
#'              The first column must contain the OTUID, ASV, or species name.
#' @param facet_orientation Whether `facet_by` appears in columns ("horizontal", default) or rows ("vertical").
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
  base_theme <- ggplot2::theme_test() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(), 
      legend.position = "none", 
      axis.title.x = ggplot2::element_text(color = "black", family = "Times New Roman", size = 12),
      axis.title.y = ggplot2::element_text(color = "black", family = "Times New Roman", size = 12),
      axis.text.y = ggplot2::element_text(colour = "black", family = "Times New Roman", size = 8),
      axis.text.x = ggplot2::element_text(colour = "black", family = "Times New Roman", size = 4)
    )
  
  # --- Relación de aspecto ---
  aspect_ratio_theme <- if (facet_orientation == "horizontal") NULL else 0.5
  
  # --- Gráfico q0 ---
  q0_vs_depth <- ggpubr::ggscatter(
    q_data, x = "Frequency", y = "q0",
    xlab = "Sequencing depth (number of reads)",
    add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
    add.params = list(color = "#B03A2E", fill = "#566573"),
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n",
      size = 4.5,
      family = "Times New Roman",
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
    add.params = list(color = "#B03A2E", fill = "#566573"),
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n",
      size = 4.5,
      family = "Times New Roman",
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
    add.params = list(color = "#B03A2E", fill = "#566573"),
    cor.coeff.args = list(
      method = "pearson",
      label.x = 3,
      label.sep = "\n",
      size = 4.5,
      family = "Times New Roman",
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
      label_fontfamily = "Times New Roman",
      label_fontface = "bold",
      label_size = 16
    )
  } else {
    cowplot::plot_grid(
      q0_vs_depth, q1_vs_depth, q2_vs_depth,
      labels = c("A", "B", "C"),
      ncol = 1,
      label_fontfamily = "Times New Roman",
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
                            fontfamily = "Times New Roman")
    } else {
      title <- cowplot::ggdraw() +
        cowplot::draw_label(plot_title,
                            fontface = "bold",
                            fontfamily = "Times New Roman")
    }
    cowplot::plot_grid(title, grid_plot, ncol = 1, rel_heights = c(0.1, 1))
  } else {
    return(grid_plot)
  }
}




