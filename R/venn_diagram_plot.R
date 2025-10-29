#' Generate a Venn diagram of taxa shared between sample groups
#'
#' This function creates a Venn diagram using either the `ggVennDiagram` or `ggenn` package
#' based on the user's preference. It includes additional customization options like
#' minimum prevalence filtering and custom group color scales (manual).
#'
#' @param table A data frame containing taxonomic abundance data with a column named `taxonomy` and subsequent columns as sample IDs.
#' @param metadata A data frame containing metadata with a column named `SAMPLEID` that matches the sample columns in `table`.
#' @param merge_by A character string specifying the metadata column by which to group and merge samples. Default is `Tratamiento`.
#' @param selected_samples Optional character vector specifying a subset of sample IDs to include in the analysis.
#' @param min_prevalence Optional numeric value (0-1) to filter taxa based on minimum prevalence across groups.
#' @param title Optional character string for the title of the plot.
#' @param method Character string: `ggVennDiagram` (default) or `ggvenn`, specifying the package to use for Venn diagram generation.
#' @param group_colors Optional vector of colors for the groups. If NULL, a default `distiller` scale with `Set3` palette will be used.

#' @return A ggplot object or other plot depending on the method.
#' @export
venn_diagram_plot <- function(table, metadata, merge_by = NULL,
                              selected_samples = NULL, min_prevalence = 0,
                              title = NULL, method = "ggvenn",
                              group_colors = NULL) {
  table <- as.data.frame(table)
  metadata <- as.data.frame(metadata)
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  if (!merge_by %in% colnames(metadata)) stop("Group or merge column is not in the metadata file.")
  
  common_samples <- intersect(colnames(table), metadata[[1]])
  table <- table[, c(common_samples), drop = FALSE]
  metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
  rownames(metadata) <- NULL
  
  if (!is.null(selected_samples)) {
    common_samples <- intersect(selected_samples, colnames(table))
    table <- table[, c(common_samples), drop = FALSE]
    metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
  }
  
  metadata_split <- split(metadata[[1]], metadata[[merge_by]])
  metadata_split <- metadata_split[sapply(metadata_split, length) > 0]
  num_groups <- length(metadata_split)
  
  # Default palette - Distiller Set3
  use_manual <- !is.null(group_colors)
  if (!use_manual) {
    group_colors <- scales::hue_pal()(num_groups)
  } else {
    group_colors <- rep(group_colors, length.out = num_groups)
  }
  
  # Crear lista según el método
  lista <- lapply(metadata_split, function(samps) {
    subset <- table[, samps, drop = FALSE]
    subset_core <- if (min_prevalence > 0) {
      subset_core <- subset[rowMeans(subset > 0) >= min_prevalence, ]
      
    } else {
      subset_core <- subset[rowSums(subset) != 0, ]

    }
    rownames(subset_core)
  })
  names(lista) <- names(metadata_split)
  
  # Selección del método
  if (method == "ggvenndiagram") {
    if (!requireNamespace("ggVennDiagram", quietly = TRUE)) stop("Install ggVennDiagram package.")
    
    if (use_manual) {
      venn_plot <- ggVennDiagram::ggVennDiagram(lista, label_alpha = 0,
                                                set_color = group_colors,
                                                edge_size = 1) +
        ggplot2::scale_fill_gradient(low = "white", high = "#5A5A5A", na.value = NA)
    } else {
      venn_plot <- ggVennDiagram::ggVennDiagram(lista, label_alpha = 0,
                                                edge_size = 1) +
        ggplot2::scale_fill_gradientn(colours = c(
          "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032",
          "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
          "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
          "#882D17", "#8DB600", "#654522", "#E25822"
        ))
    }
    
  } else if (method == "ggvenn") {
    if (!requireNamespace("ggvenn", quietly = TRUE)) stop("Install 'ggvenn' package.")
    
    if (use_manual) {
      venn_plot <- ggvenn::ggvenn(lista, fill_color = group_colors)
    } else {
      venn_plot <- ggvenn::ggvenn(lista) + 
        ggplot2::scale_fill_manual(values = c(
          "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032",
          "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
          "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
          "#882D17", "#8DB600", "#654522", "#E25822"
        ))
    }
    
  } else {
    stop("Method must be 'ggvenn' or 'ggvenndiagram'.")
  }
  
  venn_plot <- venn_plot + 
               ggtitle(title) + 
               theme(legend.position = "right",
                     legend.text = ggplot2::element_text(size = 28, family = "Times New Roman"))
  return(venn_plot)
}

