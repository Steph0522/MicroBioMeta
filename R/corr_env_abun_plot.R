



#' Title
#'
#' @param abund_table 
#' @param env_table 
#' @param method 
#' @param hc.order 
#' @param geom 
#' @param show_labels 
#' @param col_palette 
#' @param invert_axes 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#colores<- c("pink","white","purple")
#corr_env_abund_plot(abund_table = abund,
#                    env_table = env,
#                    method = "spearman",
#                    geom = "tile",
#                    hc.order = FALSE,
#                    col_palette = colores,
#                    invert_axes = TRUE,
#                    show_labels = FALSE)


corr_env_abund_plot <- function(abund_table,
                                env_table,
                                method = "spearman",
                                hc.order = TRUE,
                                geom = c("tile", "circle"),
                                show_labels = TRUE,
                                col_palette = NULL,
                                invert_axes = TRUE) {
  # Paleta por defecto
  if (is.null(col_palette)) {
    col_palette <-
      grDevices::colorRampPalette(c("blue", "white", "red"))(200)
  }
  
  # Filas comunes
  common_samples <-
    base::intersect(colnames(abund_table), rownames(env_table))
  abund <- abund_table[, common_samples, drop = FALSE]   
  env <- env_table[common_samples, , drop = FALSE]      
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
  corr_df <-
    reshape2::melt(corr_mat,
                   varnames = c("Environmental", "Group"),
                   value.name = "Correlation")
  
  # Base del gráfico
  if (invert_axes) {
    p <-
      ggplot2::ggplot(corr_df,
                      ggplot2::aes(x = Environmental, y = Group, fill = Correlation))
  } else {
    p <-
      ggplot2::ggplot(corr_df,
                      ggplot2::aes(x = Group, y = Environmental, fill = Correlation))
  }
  
  # Elegir tipo de gráfico
  if (geom == "tile") {
    p <- p + ggplot2::geom_tile(color = "gray80")
    if (show_labels) {
      p <-
        p + ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                               size = 3,
                               color = "black")
    }
  } else if (geom == "circle") {
    p <-
      p + ggplot2::geom_point(ggplot2::aes(size = abs(Correlation)),
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
  p <-
    p + ggplot2::scale_fill_gradientn(colours = col_palette,
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
