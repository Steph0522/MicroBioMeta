

venn_diagram_plot <- function(table, metadata, merge_by = NULL,
                              selected_samples = NULL, min_prevalence = 0,
                              title = NULL, method = "ggvenn",
                              group_colors = NULL) {
  # Conversión a data.frame
  table <- as.data.frame(table)
  metadata <- as.data.frame(metadata)

  # Verificación de columna de taxonomía
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if (length(tax_col) != 1) stop("There is no taxonomy column in the table.")
  
  # Verificación de merge_by
  if (!merge_by %in% colnames(metadata)) stop("Group or merge column is not in the metadata file.")
  
  # Alineación de muestras
  common_samples <- intersect(colnames(table), metadata$SAMPLEID)
  table <- table[, c(common_samples), drop = FALSE]
  metadata <- metadata[metadata$SAMPLEID %in% common_samples, , drop = FALSE]
  rownames(metadata) <- NULL
  
  # Selección de muestras específicas (si aplica)
  if (!is.null(selected_samples)) {
    common_samples <- intersect(selected_samples, colnames(table))
    table <- table[, c(common_samples), drop = FALSE]
    metadata <- metadata[metadata$SAMPLEID %in% common_samples, , drop = FALSE]
  }
  
  # Agrupamiento
  metadata_split <- split(metadata$SAMPLEID, metadata[[merge_by]])
  metadata_split <- metadata_split[sapply(metadata_split, length) > 0]
  num_groups <- length(metadata_split)
  
  # Colores
  use_manual <- !is.null(group_colors)
  if (!use_manual) {
    group_colors <- scales::hue_pal()(num_groups)
  } else {
    group_colors <- rep(group_colors, length.out = num_groups)
  }
  
  # Crear lista de OTUs por grupo
  lista <- lapply(metadata_split, function(samps) {
    subset <- table[, samps, drop = FALSE]
    subset_core <- if (min_prevalence > 0) {
      subset[rowMeans(subset > 0) >= min_prevalence, ]
    } else {
      subset[rowSums(subset) != 0, ]
    }
    rownames(subset_core)
  })
  names(lista) <- names(metadata_split)
  
  # ===============================
  # MÉTODO: ggVennDiagram
  # ===============================
  if (method == "ggvenndiagram") {
    if (!requireNamespace("ggVennDiagram", quietly = TRUE))
      stop("Install ggVennDiagram package.")
    
    venn_plot <- ggVennDiagram::ggVennDiagram(lista, label_alpha = 0, edge_size = 1) +
      ggplot2::scale_fill_gradient(low = "white", high = "#5A5A5A", na.value = NA) +
      ggplot2::theme_void(base_family = "serif", base_size = 10) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 10, family = "serif"),
        legend.text  = ggplot2::element_text(size = 10, family = "serif"),
        plot.title   = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5, family = "serif")
      )
    
  }
    # ===============================
    # MÉTODO: ggvenn
    # ===============================
 else if (method == "ggvenn") {
    if (!requireNamespace("ggvenn", quietly = TRUE))
      stop("Install 'ggvenn' package.")
    
    venn_plot <- ggvenn::ggvenn(
      lista,
      fill_color = group_colors,
      set_name_size = 10,      # Tamaño etiquetas de grupos
      text_size = 10    # Tamaño números en intersecciones
    ) +
      ggplot2::theme_void(base_family = "serif", base_size = 10) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(size = 10, family = "serif"),
        legend.text  = ggplot2::element_text(size = 10, family = "serif"),
        plot.title   = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5, family = "serif")
      )
  } else {
    stop("Method must be 'ggvenn' or 'ggvenndiagram'.")
  }
  
  venn_plot <- venn_plot + ggplot2::ggtitle(title)
  
  return(venn_plot)
}








