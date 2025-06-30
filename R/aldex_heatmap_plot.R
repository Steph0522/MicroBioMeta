#' Heatmap ALDEx
#'
#' This function generates a heatmap to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table, Data frame with taxonomy, where, the columns are the samples and rows are ASV's or taxa.
#' @param conditions, Vector that defines the categories or classes to compare.
#' @param effect, Effect size (default: >= 0.8)
#' @param pvalue_BH, Value p-ajusted (optional)
#'
#' @return A plot with the deferentially abundant taxonomic groups between two categories or groups of samples.
#' @export
#'
#' @examples aldex_plot(table = table,
#'                      conditions = conditions,
#'                      effect = 0.8,
#'                      pvalue_BH = NULL)
#'

aldex_heatmap_plot <- function(table,
                               metadata,
                               col_cond,
                               effect_threshold = 0.8,
                               pvalue_BH = NULL,
                               cluster_rows = FALSE,
                               cluster_columns = FALSE,
                               heatmap_colors = circlize::colorRamp2(c(0, 0.5, 1), c("#cbdbcd", "#759c8a", "#00544d")),
                               effect_colors = circlize::colorRamp2(c(-1.5, 0, 1.5), c("lightsalmon4", "white", "lightseagreen")),
                               pvalue_colors = list('p-value' = c("<0.001" = '#C70039', "<0.01" = '#FF5733', "<0.05" = "#FFC300", ">0.05" = "#F9E79F")),
                               treatment_colors = c("Higher" = "#808000", "Lower" = "#1B5E20")) {
  
  #check complexheatmap package
  
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("El paquete 'ComplexHeatmap' no está instalado. Instalando desde Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("ComplexHeatmap")
  }
  library(ComplexHeatmap)
  
  # Verifica que la columna de condición exista en metadata
  if (!col_cond %in% colnames(metadata)) {
    stop(paste("Column", col_cond, "not found in metadata."))
  }
  
  # Extrae condiciones
  conditions <- metadata[[col_cond]]
  unique_conditions <- unique(conditions)
  if (length(unique_conditions) != 2) {
    stop("Exactly two conditions are required for the analysis.")
  }
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  names(table)[ncol(table)] <- "taxonomy"
  
  # Prepara tabla de conteos
  table <- table %>% rownames_to_column(var = "OTUID")
  table_counts <- table %>%
    dplyr::select(-taxonomy) %>%
    tibble::column_to_rownames("OTUID")
  
  # Validar correspondencia condiciones vs muestras
  if (length(conditions) != ncol(table_counts)) {
    stop("Number of conditions does not match number of samples.")
  }
  
  # Corre ALDEx2
  aldex_results <- ALDEx2::aldex(
    reads = table_counts,
    conditions = conditions,
    mc.samples = 128,
    effect = TRUE,
    test = "t",
    verbose = TRUE,
    denom = "all",
    include.sample.summary = FALSE
  )
  
  # Verifica que existan las columnas necesarias
  if (!all(c("effect", "wi.eBH") %in% colnames(aldex_results))) {
    stop("Columns 'effect' or 'wi.eBH' missing in ALDEx2 results.")
  }
  
  # Filtra resultados según thresholds
  aldex_filtered <- aldex_results
  
  if (effect_threshold > 0 && !is.null(pvalue_BH)) {
    aldex_filtered <- aldex_results %>%
      dplyr::filter(abs(effect) >= effect_threshold, wi.eBH <= pvalue_BH)
  } else if (effect_threshold > 0) {
    aldex_filtered <- aldex_results %>%
      dplyr::filter(abs(effect) >= effect_threshold)
  } else if (!is.null(pvalue_BH)) {
    aldex_filtered <- aldex_results %>%
      dplyr::filter(wi.eBH <= pvalue_BH)
  }
  
  # Prepare data for heatmap
  aldex_plot <- aldex_filtered %>%
    tibble::rownames_to_column("OTUID") %>%
    dplyr::left_join(table %>% dplyr::select(OTUID, taxonomy), by = "OTUID") %>%
    dplyr::mutate(
      seccion = dplyr::case_when(
        diff.btw < 0 ~ paste("Lower in", unique_conditions[2]),
        diff.btw > 0 ~ paste("Higher in", unique_conditions[1]),
        TRUE ~ "No Change"
      ),
      taxonomy = dplyr::case_when(
        grepl("g__[^;]*", taxonomy) &
          !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
        grepl("f__[^;]*", taxonomy) &
          !grepl("f__uncultured|f__$", taxonomy) ~ paste0(
            "other ",
            stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)
          ),
        grepl("o__[^;]*", taxonomy) &
          !grepl("o__uncultured|o__$", taxonomy) ~ paste0(
            "other ",
            stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)
          ),
        grepl("c__[^;]*", taxonomy) &
          !grepl("c__uncultured|c__$", taxonomy) ~ paste0(
            "other ",
            stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)
          ),
        grepl("p__[^;]*", taxonomy) &
          !grepl("p__uncultured|p__$", taxonomy) ~ paste0(
            "other ",
            stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)
          ),
        TRUE ~ "Unclassified"
      ),
      taxonomy = stringr::str_trim(taxonomy),
      taxonomy = make.unique(taxonomy),
      p.value = dplyr::case_when(
        wi.eBH <= 0.001 ~ "<0.001",
        wi.eBH <= 0.01 ~ "<0.01",
        wi.eBH <  0.05 ~ "<0.05",
        TRUE ~ ">0.05"
      )
    ) %>%
    dplyr::arrange(diff.btw)
  
  rab_cols <- paste0("rab.win.", unique_conditions)
  heat_data <- aldex_plot %>%
    dplyr::select(taxonomy, all_of(rab_cols)) %>%
    dplyr::rename_with(~ unique_conditions, all_of(rab_cols)) %>%
    tibble::column_to_rownames(var = "taxonomy") %>%
    as.matrix()
  
  treatment_colors_full <- setNames(
    treatment_colors,
    c(
      paste("Higher in", unique_conditions[1]),
      paste("Lower in", unique_conditions[2])
    )
  )
  
  barpl <- ComplexHeatmap::rowAnnotation(
    "difference \nbetween groups" = ComplexHeatmap::anno_barplot(
      aldex_plot$diff.btw,
      which = "row",
      gp = grid::gpar(fill = treatment_colors_full[aldex_plot$seccion]),
      width = unit(4, "cm")
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 8),
    annotation_name_rot = 0
  )
  
  annP <- ComplexHeatmap::rowAnnotation(
    "p-value" = aldex_plot$p.value,
    simple_anno_size = unit(0.45, "cm"),
    annotation_name_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    annotation_legend_param = list(
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8),
      direction = "vertical"
    ),
    col = pvalue_colors,
    show_legend = TRUE,
    gp = grid::gpar(col = "white"),
    show_annotation_name = TRUE
  )
  
  

  left_annotation <- ComplexHeatmap::rowAnnotation(
    "Effect size" = aldex_plot$effect,
    col = list("Effect size" = effect_colors),
    simple_anno_size = unit(0.45, "cm"),
    annotation_name_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    annotation_legend_param = list(
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8),
      direction = "vertical"
    ),
    show_legend = TRUE,
    gp = grid::gpar(col = "white"),
    show_annotation_name = TRUE
  )
  
  
  
  heatmap <- ComplexHeatmap::Heatmap(
    heat_data,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    width = grid::unit(ncol(heat_data) * 7, "mm"),
    height = grid::unit(nrow(heat_data) * 6, "mm"),
    column_names_rot = 90,
    rect_gp = grid::gpar(col = "white", lwd = 2),
    left_annotation = left_annotation,
    right_annotation = barpl,
    name = "Median clr value",
    heatmap_legend_param = list(
      direction = "vertical",
      labels_gp = grid::gpar(fontsize = 8),
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      legend_height = unit(2, "cm")
    ),
    column_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    col = heatmap_colors,
    row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
    show_heatmap_legend = TRUE
  )
  
  ComplexHeatmap::draw(heatmap,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right")
  return(invisible(heatmap))
}
