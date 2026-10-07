#' Alpha diversity correlation plot
#'
#' Computes Hill numbers (q = 0, 1, 2) per sample and plots each against
#' sequencing depth (total reads) as a scatter plot with a fitted regression
#' line and correlation coefficient, combining the three plots (q0,
#' q1, q2) into a single ggplot2 plot
#'
#' @param table A data frame with taxa in rows and samples in columns. 
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param method Character. Correlation method passed to
#'   \code{ggpubr::stat_cor}. One of \code{"spearman"} (default, rank-based
#'   and robust to non-linear/non-normal relationships), \code{"pearson"}, or
#'   \code{"kendall"}.
#' @param facet_orientation Whether the three q0/q1/q2 panels are arranged in
#'   a row ("horizontal", default) or a column ("vertical").
#' @param title Character. Title for the combined figure. \code{"default"}
#'   (default) shows "Alpha diversity vs sequencing depth"; \code{"none"}
#'   shows no title; any other string is used as-is.
#' @param panel_label_case Character. Case of the auto-generated A/B/C panel
#'   tags. One of \code{"upper"} (default, "A", "B", "C") or \code{"lower"}
#'   ("a", "b", "c"). Ignored if \code{panel_labels} is supplied.
#' @param panel_labels Optional character vector of 3 custom panel tags (one
#'   per q0/q1/q2 panel), used as-is (e.g. \code{c("(a)", "(b)", "(c)")} or
#'   \code{c("a.", "b.", "c.")}) — for journal styles that
#'   \code{panel_label_case} alone can't produce. Overrides
#'   \code{panel_label_case} when provided.
#' @param panel_label_bold Logical. If \code{TRUE} (default), panel tags are
#'   bold. Set to \code{FALSE} for journals that require plain (non-bold)
#'   panel tags.
#' @param save_table Logical. If \code{TRUE}, saves the Hill numbers table to
#'   disk. Default \code{FALSE}.
#' @param x_axis_title Title for the x-axis of the three panels. Default
#'   \code{"Sequencing depth (number of reads)"}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"hill.txt"}.
#'
#' @return A ggplot object showing alpha diversity with Hill numbers.
#' @export
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' # Always tagged A/B/C (see panel_label_case and panel_labels to customize)
#' alpha_hill_corr_plot(
#'   table             = table,
#'   facet_orientation = "horizontal"
#' )
#'
#' ## Using Pearson correlation instead of the default Spearman
#' alpha_hill_corr_plot(
#'   table  = table,
#'   method = "pearson"
#' )

alpha_hill_corr_plot <- function(table,
                                method = "spearman",
                                facet_orientation = "horizontal",
                                title = "default",
                                panel_label_case = "upper",
                                panel_labels = NULL,
                                panel_label_bold = TRUE,
                                save_table = FALSE,
                                table_filename = "hill.txt",
                                x_axis_title = "Sequencing depth (number of reads)") {

  table <- table[, !colnames(table) %in% "taxonomy"]
  table <- data.frame(t(table))

  q_data <- data.frame(
    SampleID = rownames(table),
    Frequency = rowSums(table),
    q0 = hillR::hill_taxa(comm = table, q = 0),
    q1 = hillR::hill_taxa(comm = table, q = 1),
    q2 = hillR::hill_taxa(comm = table, q = 2)
  )

  if (save_table) {
    utils::write.table(q_data, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message("Table saved as: ", table_filename)
  }

  base_theme <- .mbm_theme(
    legend_position = "none",
    extra = ggplot2::theme(
      legend.title = ggplot2::element_blank()
    )
  )
  
  aspect_ratio_theme <- NULL
  
  q0_vs_depth <- ggpubr::ggscatter(
    q_data, x = "Frequency", y = "q0",
    xlab = x_axis_title,
    add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
    add.params = list(color = "#D55E00", fill = "#56B4E9"),
    cor.coeff.args = list(
      method = method,
      label.x = 3,
      label.sep = "\n",
      p.accuracy = 0.001,
      r.accuracy = 0.001,
      size = 4.5,
      family = "serif",
      fontface = "bold",
      color = "black"
    )
  ) +
    base_theme +
    ggplot2::labs(y = expression(paste(italic("q"), "=0", " (number of total features)"))) +
    ggplot2::theme(aspect.ratio = aspect_ratio_theme)
  
  q1_vs_depth <- ggpubr::ggscatter(
    q_data, x = "Frequency", y = "q1",
    xlab = x_axis_title,
    add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
    add.params = list(color = "#D55E00", fill = "#56B4E9"),
    cor.coeff.args = list(
      method = method,
      label.x = 3,
      label.sep = "\n",
      p.accuracy = 0.001,
      r.accuracy = 0.001,
      size = 4.5,
      family = "serif",
      fontface = "bold",
      color = "black"
    )
  ) +
    base_theme +
    ggplot2::labs(y = expression(paste(italic("q"), "=1", " (number of frequent features)"))) +
    ggplot2::theme(aspect.ratio = aspect_ratio_theme)
  
  q2_vs_depth <- ggpubr::ggscatter(
    q_data, x = "Frequency", y = "q2",
    xlab = x_axis_title,
    add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
    add.params = list(color = "#D55E00", fill = "#56B4E9"),
    cor.coeff.args = list(
      method = method,
      label.x = 3,
      label.sep = "\n",
      p.accuracy = 0.001,
      r.accuracy = 0.001,
      size = 4.5,
      family = "serif",
      fontface = "bold",
      color = "black"
    )
  ) +
    base_theme +
    ggplot2::labs(y = expression(paste(italic("q"), "=2", " (number of dominant features)"))) +
    ggplot2::theme(aspect.ratio = aspect_ratio_theme)
  
  
  resolved_labels <- if (!is.null(panel_labels)) {
    panel_labels
  } else if (identical(panel_label_case, "lower")) c("a", "b", "c") else c("A", "B", "C")

  resolved_title <- if (!is.character(title) || tolower(title) == "none") {
    NULL
  } else if (tolower(title) == "default") {
    "Alpha diversity vs sequencing depth"
  } else {
    title
  }

  horizontal <- identical(facet_orientation, "horizontal")
  .mbm_patchwork_grid(
    list(q0_vs_depth, q1_vs_depth, q2_vs_depth),
    ncol = if (horizontal) 3 else 1,
    tags = resolved_labels,
    title = resolved_title,
    show_legend = FALSE,
    tag_bold = panel_label_bold
  )
}
