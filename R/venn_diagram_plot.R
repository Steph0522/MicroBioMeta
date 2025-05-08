#' Generate a Venn diagram of taxa shared between sample groups
#'
#' This function creates a Venn diagram using either the `microeco` or `ggVennDiagram` package
#' based on the user's preference. It includes additional customization options like
#' minimum prevalence filtering and custom group color scales (manual).
#'
#' @param table A data frame containing taxonomic abundance data with a column named `taxonomy` and subsequent columns as sample IDs.
#' @param metadata A data frame containing metadata with a column named `SAMPLEID` that matches the sample columns in `table`.
#' @param merge_by A character string specifying the metadata column by which to group and merge samples. Default is `Tratamiento`.
#' @param selected_samples Optional character vector specifying a subset of sample IDs to include in the analysis.
#' @param min_prevalence Optional numeric value (0-1) to filter taxa based on minimum prevalence across groups.
#' @param title Optional character string for the title of the plot.
#' @param method Character string: `microeco` (default) or `ggvenn`, specifying the package to use for Venn diagram generation.
#' @param group_colors Optional vector of colors for the groups. If NULL, a default `distiller` scale with `Set3` palette will be used.
#' @return A ggplot object or other plot depending on the method.
#' @export
venn_diagram_plot <- function(table, metadata, merge_by = "Tratamiento",
                                    selected_samples = NULL, min_prevalence = 0,
                                    title = NULL, method = "microeco",
                                    group_colors = NULL) {
  table <- as.data.frame(table)
  metadata <- as.data.frame(metadata)

  if (!"taxonomy" %in% colnames(table)) stop("La tabla debe contener 'taxonomy'.")
  if (!merge_by %in% colnames(metadata)) stop("La columna de agrupamiento no existe.")

  common_samples <- intersect(colnames(table)[-1], metadata$SAMPLEID)
  table <- table[, c("taxonomy", common_samples), drop = FALSE]
  metadata <- metadata[metadata$SAMPLEID %in% common_samples, , drop=FALSE]
  rownames(metadata) <- NULL

  if (!is.null(selected_samples)) {
    common_samples <- intersect(selected_samples, colnames(table))
    table <- table[, c("taxonomy", common_samples), drop = FALSE]
    metadata <- metadata[metadata$SAMPLEID %in% common_samples, , drop=FALSE]
  }

  metadata_split <- split(metadata$SAMPLEID, metadata[[merge_by]])
  num_groups <- length(metadata_split)

  if (min_prevalence > 0) {
    table <- table %>%
      dplyr::rowwise() %>%
      dplyr::mutate(prev = sum(c_across(-taxonomy) > 0) / length(metadata$SAMPLEID)) %>%
      dplyr::filter(prev >= min_prevalence) %>%
      dplyr::select(-prev)
  }

  # Default palette - Distiller Set3
  use_manual <- !is.null(group_colors)
  if (!use_manual) {
    # default distiller palette
    group_colors <- scales::hue_pal()(num_groups)
  } else {
    group_colors <- rep(group_colors, length.out = num_groups)
  }

  if (method == "ggvenn") {
    if (!requireNamespace("ggVennDiagram", quietly=TRUE)) stop("Instala ggVennDiagram.")
    taxa_list <- lapply(metadata_split, function(samps) unique(table$taxonomy[rowSums(table[,samps,drop=FALSE]>0)>0]))
    names(taxa_list) <- names(metadata_split)
    # Crear conjuntos con filas no nulas dinámicamente para cada grupo
    lista <- lapply(metadata_split, function(samps) {
      subset <- table[, samps, drop = FALSE]
      subset_core <- subset[rowSums(subset) != 0, ]
      rownames(subset_core)
    })
    names(lista) <- names(metadata_split)
    venn_plot <- ggVennDiagram::ggVennDiagram(lista, label_alpha = 0) 
    
    if (use_manual) {
      venn_plot <- venn_plot + ggplot2::scale_fill_manual(values = group_colors)
    } else {
      venn_plot <- venn_plot + ggplot2::scale_fill_distiller(palette = "Set3", direction = 1)
    }
    

  } else if (method == "microeco") {
    if (!requireNamespace("microeco", quietly=TRUE)) stop("Instala microeco.")
    abund <- table[, -1, drop=FALSE]
    samp_df <- tibble::column_to_rownames(metadata, "SAMPLEID")
    ds <- microeco::microtable$new(abund, samp_df, auto_tidy=TRUE)
    merged <- ds$merge_samples(merge_by)
    tv <- microeco::trans_venn$new(merged, ratio = "seqratio")
    if (nrow(tv$data_summary)==0) stop("No hay datos para Venn.")
    
    if (use_manual) {
      venn_plot <- tv$plot_venn(color_circle = group_colors,linesize = 3)
    } else {
      venn_plot <- tv$plot_venn(color_circle = RColorBrewer::brewer.pal(10, "Set3"),  linesize = 3)
    }
    
  } else stop("Método debe ser 'microeco' o 'ggvenn'.")

  venn_plot <- venn_plot + ggtitle(title) #+ theme_minimal(base_size=14)
   venn_plot <- venn_plot + theme(legend.position="none")
  return(venn_plot)
}
