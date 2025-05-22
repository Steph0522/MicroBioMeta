#' CCA Biplot with ggplot2
#'
#' Performs Canonical Correspondence Analysis (CCA) based on a species abundance table
#' and selected environmental variables, returning a biplot with ggplot2 that visualizes
#' sample scores and significant environmental vectors.
#'
#' @param table A data frame or matrix of species abundances (samples as rows, species as columns).
#' @param env_data A data frame of environmental variables (rows must match `table`).
#' @param env_vars A character vector with the names of environmental variables to include in the CCA.
#' @param method Transformation method passed to `decostand` (default is `"hell"` for Hellinger).
#' @param metadata Optional data frame with sample metadata for grouping in the plot.
#' @param group_col Optional name of the column in `metadata` used to define sample groups.
#' @param group_colors Optional named vector of colors to use for each group.
#' @param legend_title Optional custom title for the group legend.
#' @param scale_env Logical; whether to scale environmental variables (default is `TRUE`).
#' @param pval_threshold P-value threshold for selecting significant environmental variables (default is `0.05`).
#' @param seed Random seed for reproducibility (default is `126`).
#' @param scale_arrows Numeric value to scale environmental vectors in the plot.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object displaying the CCA biplot with sample scores and significant environmental vectors.
#' @export
#'
#' @examples
#' # cca_biplot(table, env_data, env_vars = c("pH", "Temp"),
#' #            metadata, group_col = "Treatment")

cca_biplot <- function(table,
                       env_data,
                       env_vars,
                       method = "hell",
                       metadata = NULL,
                       group_col = NULL,
                       group_colors = NULL,
                       legend_title = NULL,
                       scale_env = TRUE,
                       pval_threshold = 0.05,
                       seed = 126,
                       scale_arrows = 1,
                       title = NULL) {
  # Cargar paquetes necesarios
  require(vegan)
  require(ggplot2)
  require(ggvegan)
  require(dplyr)
  
  # 1. Transformación
  spp_hell <- vegan::decostand(table, method = method)
  
  # 2. Escalar variables ambientales seleccionadas
  if (scale_env) {
    env_scaled <- scale(env_data[, env_vars], scale = TRUE, center = FALSE) %>%
      as.data.frame()
  } else {
    env_scaled <- env_data[, env_vars]
  }
  
  # 3. Asegurar que las filas coincidan entre species_table, env_data y metadata
  if (!is.null(metadata)) {
    metadata <- metadata[rownames(species_table), , drop = FALSE]
  }
  
  # 4. Análisis CCA
  set.seed(seed)
  cca_result <- vegan::cca(spp_hell ~ ., data = env_scaled)
  
  # 5. Ajuste de vectores ambientales
  fit <- vegan::envfit(cca_result, env_scaled)
  
  # 6. Selección de variables significativas
  sig_vars <- names(which(fit$vectors$pvals < pval_threshold))
  
  if (length(sig_vars) == 0) {
    warning("No hay variables ambientales significativas (p <", pval_threshold, ")")
    return(ggplot() + theme_void() + ggtitle("Sin variables significativas"))
  }
  
  vectors_scores <- vegan::scores(fit, display = "vectors")[sig_vars, , drop = FALSE] %>%
    as.data.frame()
  vectors_scores$Variable <- rownames(vectors_scores)
  
  # 7. Coordenadas de sitios
  site_scores <- vegan::scores(cca_result, display = "sites") %>%
    as.data.frame()
  site_scores$SampleID <- rownames(site_scores)
  
  if (!is.null(metadata) &&
      !is.null(group_col) && group_col %in% colnames(metadata)) {
    
    site_scores$Group <- metadata[[group_col]]
    
    # Paleta discreta por defecto (ColorBrewer Set2)
    default_colors <- c(
      "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
    )
    
    groups_present <- unique(site_scores$Group)
    
    if (is.null(group_colors)) {
      color_values <- rep(default_colors, length.out = length(groups_present))
      names(color_values) <- groups_present
      group_colors <- color_values
    }
    
    # Título de leyenda personalizado si se proporciona
    legend_name <- ifelse(is.null(legend_title), group_col, legend_title)
    
    plot <- ggplot(site_scores, aes(x = CCA1, y = CCA2, color = Group)) +
      geom_point(size = 3) +
      scale_color_manual(name = legend_name, values = group_colors)
    
  } else {
    plot <- ggplot(site_scores, aes(x = CCA1, y = CCA2)) +
      geom_point(size = 3)
  }
  
  # 8. Añadir vectores ambientales escalados
  plot <- plot +
    geom_segment(
      data = vectors_scores,
      aes(
        x = 0, y = 0,
        xend = CCA1 * scale_arrows,
        yend = CCA2 * scale_arrows
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = vectors_scores,
      aes(x = CCA1 * scale_arrows, y = CCA2 * scale_arrows, label = Variable),
      color = "black",
      hjust = 0.5,
      vjust = -0.5,
      inherit.aes = FALSE
    )
  
  # 9. Título del gráfico si se proporciona
  if (!is.null(title)) {
    plot <- plot + ggtitle(title)
  }
  
  plot <- plot + theme_minimal()
  
  return(plot)
}


