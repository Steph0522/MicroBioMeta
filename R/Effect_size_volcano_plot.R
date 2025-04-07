#' Volcano plot of effect size
#'
#' This function generates a volcano plot based on the results of ALDEx2,
#' displaying the effect size on the x-axis and the Kruskal-Wallis p-value on the y-axis.

#' @param table Data frame with count data; columns represent samples, rows represent features.
#' @param metadata Data frame containing metadata for the samples.
#' @param col_cond Name of the column in `metadata` that contains the experimental conditions.
#' @param col_inf Color for points with effect size lower than `threshold_lower`.
#' @param col_sup Color for points with effect size higher than `threshold_upper`.
#' @param threshold_lower Lower threshold for effect size (x-axis).
#' @param threshold_upper Upper threshold for effect size (x-axis).
#' @param cond Name of the condition that appears first in `table` (used in plot labels).
#' @param show_labels Logical. Whether to display "Higher/Lower in cond" labels (default is TRUE).
#' 
#' #' @return A `ggplot` object with the volcano plot.
#' @export
#'
#' @examples

#' 
#' 
#' 
effect_size_plot <- function(table,
                             metadata,
                             col_cond,
                             col_inf,
                             col_sup,
                             threshold_lower,
                             threshold_upper,
                             cond,
                             show_labels = TRUE) {
  library(ALDEx2)
  library(ggplot2)
  library(scales)
  library(ggtext)
  
  condiciones <- metadata[[col_cond]]
  
  aldex_clr <- aldex.clr(table, condiciones, mc.samples = 128, denom = "all")
  
  effect_size <- aldex.effect(
    aldex_clr,
    verbose = TRUE,
    include.sample.summary = FALSE,
    useMC = FALSE,
    CI = FALSE
  )
  
  KW <- aldex.kw(aldex_clr, useMC = FALSE, verbose = FALSE)
  
  resultado <- cbind(effect_size, KW)
  
  resultado$grupo <- ifelse(
    resultado$effect <= threshold_lower,
    "Menor en condicion 1",
    ifelse(
      resultado$effect >= threshold_upper,
      "Mayor en condicion 1",
      "Normal"
    )
  )
  
  p <- ggplot(resultado, aes(x = effect, y = kw.ep, color = grupo)) +
    geom_point(size = 3.5) +
    scale_color_manual(
      values = c(
        "Menor en condicion 1" = col_inf,
        "Normal" = "gray",
        "Mayor en condicion 1" = col_sup
      )
    ) +
    geom_vline(xintercept = c(threshold_lower, threshold_upper), linetype = 2, color = "black") +
    geom_hline(yintercept = 0.05, linetype = 2, color = "black") +
    labs(
      x = "Effect size",
      y = bquote(italic("p") ~ " value"),
      color = "Grupo"
    ) +
    ylim(c(0, 1)) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 12)
    ) +
    theme(legend.position = "none")
  
  if (show_labels) {
    p <- p +
      geom_richtext(
        aes(x = threshold_lower, y = 0.95,
            label = paste("<b>Menor en", cond, "</b>")),
        color = col_inf,
        size = 5,
        hjust = 1,
        vjust = 1,
        inherit.aes = FALSE
      ) +
      geom_richtext(
        aes(x = threshold_upper, y = 0.95,
            label = paste("<b>Mayor en", cond, "</b>")),
        color = col_sup,
        size = 5,
        hjust = 0,
        vjust = 1,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}