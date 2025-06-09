#' Generate a Venn diagram of taxa shared between sample groups
#'
#' This function creates a Venn diagram using either the `microeco` or `ggVennDiagram` package
#' based on the user's preference. It shows the overlap of taxa among groups defined in the metadata.
#'
#' @param table A data frame containing taxonomic abundance data with a column named `"taxonomy"` and subsequent columns as sample IDs.
#' @param metadata A data frame containing metadata with a column named `"SAMPLEID"` that matches the sample columns in `table`.
#' @param merge_by A character string specifying the metadata column by which to group and merge samples. Default is `"Tratamiento"`.
#' @param selected_samples Optional character vector specifying a subset of sample IDs to include in the analysis.
#' @param title Optional character string for the title of the plot.
#' @param method Character string: `"microeco"` (default) or `"ggvenn"`, specifying the package to use for Venn diagram generation.
#'
#' @return A `ggplot` object or other plot depending on the method.
#' @export
venn_diagram <- function(table,
                         metadata,
                         merge_by = "Tratamiento",
                         selected_samples = NULL,
                         title = NULL,
                         method = "microeco") {
  
  library(dplyr)
  library(tibble)
  library(ggplot2)
  
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
  
  # Sincronizar muestras si selected_samples no está definido
  if (is.null(selected_samples)) {
    sample_ids_table <- colnames(table)[-1]
    sample_ids_metadata <- metadata$SAMPLEID
    common_samples <- intersect(sample_ids_table, sample_ids_metadata)
    
    removed_from_table <- setdiff(sample_ids_table, common_samples)
    removed_from_metadata <- setdiff(sample_ids_metadata, common_samples)
    
    # Mostrar mensaje si se eliminan muestras
    if (length(removed_from_table) > 0 || length(removed_from_metadata) > 0) {
      message("Se han eliminado muestras para sincronizar 'table' y 'metadata':")
      if (length(removed_from_table) > 0) {
        message(" - De la tabla: ", paste(removed_from_table, collapse = ", "))
      }
      if (length(removed_from_metadata) > 0) {
        message(" - De los metadatos: ", paste(removed_from_metadata, collapse = ", "))
      }
    }
    
    table <- table[, c("taxonomy", common_samples), drop = FALSE]
    metadata <- dplyr::filter(metadata, SAMPLEID %in% common_samples)
    
  } else {
    # Filtrar si se especifica un subconjunto
    selected_samples <- intersect(selected_samples, colnames(table))
    metadata <- dplyr::filter(metadata, SAMPLEID %in% selected_samples)
    table <- table[, c("taxonomy", selected_samples), drop = FALSE]
  }
  
  # Verificar columna de agrupación
  if (!(merge_by %in% colnames(metadata))) {
    stop(paste("La columna", merge_by, "no existe en los metadatos."))
  }
  
  if (method == "ggvenn") {
    if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
      stop("El paquete 'ggVennDiagram' no está instalado.")
    }
    
    metadata_split <- split(metadata$SAMPLEID, metadata[[merge_by]])
    taxa_list <- lapply(metadata_split, function(samples) {
      sub_table <- table[, c("taxonomy", samples), drop = FALSE]
      present_taxa <- sub_table$taxonomy[rowSums(as.matrix(sub_table[,-1, drop = FALSE]) > 0) > 0]
      unique(present_taxa)
    })
    
    # Eliminar grupos vacíos
    taxa_list <- taxa_list[sapply(taxa_list, length) > 0]
    
    if (length(taxa_list) < 2 || length(taxa_list) > 7) {
      stop("El número de grupos con datos debe estar entre 2 y 7 para usar 'ggvenn'.")
    }
    
    venn_plot <- ggVennDiagram::ggVennDiagram(taxa_list, label_alpha = 0)
    
    if (!inherits(venn_plot, "gg")) {
      warning("El resultado no es un objeto ggplot. Se devolverá sin personalizar.")
      return(venn_plot)
    }
    
    venn_plot +
      ggtitle(title) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.title = element_blank())
    
  } else if (method == "microeco") {
    if (!requireNamespace("microeco", quietly = TRUE)) {
      stop("El paquete 'microeco' no está instalado.")
    }
    
    abund_table <- table[, -1]
    metadata_df <- tibble::column_to_rownames(metadata, "SAMPLEID")
    
    dataset <- microeco::microtable$new(otu_table = abund_table,
                                        sample_table = metadata_df,
                                        auto_tidy = TRUE)
    
    dataset_merged <- dataset$merge_samples(merge_by)
    t1 <- microeco::trans_venn$new(dataset_merged, ratio = NULL)
    
    if (nrow(t1$data_summary) == 0) {
      stop("No hay datos suficientes para generar el diagrama de Venn.")
    }
    
    t1$plot_venn() +
      scale_fill_gradient(low = "lightyellow", high = "red") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.title = element_blank()) +
      ggtitle(title)
    
  } else {
    stop("El parámetro 'method' debe ser 'microeco' o 'ggvenn'.")
  }
}
