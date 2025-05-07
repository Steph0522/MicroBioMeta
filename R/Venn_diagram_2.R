#' Title
#'
#' @param table 
#' @param metadata 
#' @param merge_by 
#' @param selected_samples 
#' @param min_prevalence 
#' @param title 
#' @param show_legend 
#' @param group_colors 
#'
#' @return
#' @export
#'
#' @examples
venn_diagram <- function(table,
                         metadata,
                         merge_by = "Tratamiento",
                         selected_samples = NULL,
                         min_prevalence = 0,
                         title = NULL,
                         show_legend = TRUE,
                         group_colors = NULL) {
  
  table <- as.data.frame(table)
  metadata <- as.data.frame(metadata)
  
  
  # Determinar muestras a usar, si se da el vector selected_samples usa esas, si no 
  #usa todas las muestras de la tabla 
  all_samples <- if (is.null(selected_samples)) {
    metadata[,1]
  } else {
    selected_samples
  }
  
  # Sincronizar las muestras de la tabla y la metadata
  samples_in_both <- base::intersect(all_samples, base::colnames(table))
  if (length(samples_in_both) < length(all_samples)) {
    message("Se eliminaron ", length(all_samples) - length(samples_in_both), " muestras del metadata que no están en la tabla.")
  }
  if (length(samples_in_both) < (ncol(table) - 1)) {
    message("Se eliminaron ", (ncol(table) - 1) - length(samples_in_both), " muestras de la tabla que no están en el metadata.")
  }
  
  # Filtrar y ordenar
  metadata <- dplyr::filter(metadata, metadata[,1] %in% samples_in_both)
  metadata <- dplyr::arrange(metadata, base::match(metadata[,1], samples_in_both))
  table <- table[, c("taxonomy", samples_in_both), drop = FALSE]
  
  # Agrupar muestras por grupo basado en el parámetro merge_by
  metadata_split <- base::split(metadata[,1], metadata[[merge_by]])
  
  # Si min_prevalence > 0, aplicamos el filtro de prevalencia
  if (min_prevalence > 0) {
    prevalent_taxa <- base::lapply(metadata_split, function(samples) {
      sub_table <- table[, c("taxonomy", samples), drop = FALSE]
      mat <- base::as.matrix(sub_table[, -1, drop = FALSE]) > 0
      prevalence <- base::rowSums(mat) / base::ncol(mat)
      sub_table$taxonomy[prevalence >= min_prevalence]
    })
    
    keep_taxa <- base::unique(base::unlist(prevalent_taxa))
    table <- table[table$taxonomy %in% keep_taxa, , drop = FALSE]
  }
  
  # Identificar taxones presentes en cada grupo
  taxa_list <- base::lapply(metadata_split, function(samples) {
    sub_table <- table[, c("taxonomy", samples), drop = FALSE]
    mat <- base::as.matrix(sub_table[, -1, drop = FALSE]) > 0
    present_taxa <- sub_table$taxonomy[base::rowSums(mat) > 0]
    return(present_taxa)
  })
  
  # Si no se proporcionan colores, usar colores predeterminados
  if (is.null(group_colors)) {
    # Si hay más de 3 grupos, usamos una paleta con más colores
    n_groups <- length(taxa_list)
    if (n_groups >= 3) {
      group_colors <- RColorBrewer::brewer.pal(n_groups, "Set3")
    } else {
      # Para menos de 3 grupos, asignamos colores manualmente
      group_colors <- c("black", "gray19", "gray8")[1:n_groups]
    }
  }
  
  # Generar Venn
  venn_plot <- ggVennDiagram::ggVennDiagram(taxa_list,
                                            label_alpha = 0,
                                            set_color = group_colors) +
    
    ggplot2::scale_fill_gradient(low = "white", high = "grey") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text = ggplot2::element_text(size = 18,
                                        color = "black",
                                        face = "bold",
                                        family = "Arial") 
    )
  
  if (length(taxa_list) == 2) {
    venn_plot <- venn_plot + ggplot2::coord_flip()
  }
  
  # Mostrar u ocultar leyenda
  if (!show_legend) {
    venn_plot <- venn_plot + ggplot2::theme(legend.position = "none")
  }
  
  return(venn_plot)
}
