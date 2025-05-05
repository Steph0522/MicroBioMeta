#' Generate a Venn diagram of taxa shared between sample groups
#'
#' This function creates a Venn diagram using taxonomic abundance data and metadata,
#' showing the overlap of taxa among groups defined in the metadata. It uses the `microeco` package
#' to preprocess and merge samples by group, and `ggplot2` for visualization customization.
#'
#' @param table A data frame containing taxonomic abundance data with a column named `"taxonomy"` and subsequent columns as sample IDs.
#' @param metadata A data frame containing metadata with a column named `"SAMPLEID"` that matches the sample columns in `table`.
#' @param merge_by A character string specifying the metadata column by which to group and merge samples. Default is `"Tratamiento"`.
#' @param selected_samples Optional character vector specifying a subset of sample IDs to include in the analysis.
#' @param title Optional character string for the title of the plot.
#'
#' @return A `ggplot` object representing the Venn diagram.
#' @export
#'
#' @examples
#' 
#' 
venn_diagram <- function(table,
                         metadata,
                         merge_by = "Tratamiento",
                         selected_samples = NULL,
                         title = NULL) {
  
  # Asegurar que los datos son data frames
  table <- as.data.frame(table)
  metadata <- as.data.frame(metadata)
  
  # Verificar que la columna 'taxonomy' existe
  if (!"taxonomy" %in% colnames(table)) {
    stop("La tabla no contiene una columna llamada 'taxonomy'")
  }
  
  # Reordenar la tabla para colocar 'taxonomy' al principio
  table <- table[, c("taxonomy", setdiff(colnames(table), "taxonomy"))]
  
  # Eliminar entradas no informativas
  table <- dplyr::filter(table, taxonomy != "d__Bacteria;__;__;__;__;__")
  
  # Ordenar muestras en la tabla según metadata
  ordered_samples <- intersect(metadata$SAMPLEID, colnames(table)[-1])
  table <- table[, c("taxonomy", ordered_samples)]
  
  # Filtrar si se especifica un subconjunto
  if (!is.null(selected_samples)) {
    selected_samples <- intersect(selected_samples, colnames(table))
    metadata <- dplyr::filter(metadata, SAMPLEID %in% selected_samples)
    table <- table[, c("taxonomy", selected_samples), drop = FALSE]
  }
  
  # Verificar que merge_by esté presente en los metadatos
  if (!(merge_by %in% colnames(metadata))) {
    stop(paste("La columna", merge_by, "no existe en los metadatos."))
  }
  
  # Separar la columna taxonomy (si no se usará) y dejar solo las abundancias
  table <- table[, -1]
  
  # Establecer SAMPLEID como rownames
  metadata <- tibble::column_to_rownames(metadata, "SAMPLEID")
  
  # Crear objeto microeco
  dataset <- microeco::microtable$new(otu_table = table,
                                      sample_table = metadata,
                                      auto_tidy = TRUE)
  
  # Combinar muestras por grupo
  dataset_merged <- dataset$merge_samples(merge_by)
  t1 <- microeco::trans_venn$new(dataset_merged, ratio = NULL)
  
  # Verificar que hay datos para graficar
  if (nrow(t1$data_summary) == 0) {
    stop("No hay datos suficientes para generar el diagrama de Venn.")
  }
  
  # Crear gráfico
  venn_plot <- t1$plot_venn()
  
  # Personalizar con ggplot2
  venn_plot +
    ggplot2::scale_fill_gradient(low = "lightyellow", high = "red") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::ggtitle(title)
}
