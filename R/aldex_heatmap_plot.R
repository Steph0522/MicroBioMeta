#' ALDEx2 differential abundance heatmap
#'
#' Runs ALDEx2 on a counts table and metadata and returns a ComplexHeatmap
#' showing differentially abundant taxa, their effect size, p-value, and
#' difference between groups.
#'
#' @param table Data frame with taxa as rows and samples as columns. Must
#'   contain exactly one taxonomy column (named "taxonomy", "Taxonomy",
#'   "taxon", "taxa", "Taxa", or "Taxon").
#' @param metadata Data frame with one row per sample. Must contain the column
#'   specified in \code{col_cond}.
#' @param col_cond Character. Name of the column in \code{metadata} that
#'   defines the two groups to compare. Exactly two unique values are required.
#' @param effect_threshold Numeric. Minimum absolute effect size to retain
#'   (default \code{0.8}).
#' @param pvalue_BH Numeric or NULL. Maximum BH-adjusted p-value to retain.
#'   If NULL (default) only \code{effect_threshold} is applied.
#' @param cluster_rows Logical. Cluster heatmap rows (default \code{FALSE}).
#' @param cluster_columns Logical. Cluster heatmap columns (default \code{FALSE}).
#' @param diverging_palette Character. Name of a built-in colorblind-friendly
#'   diverging palette for the main heatmap body. One of \code{"BuOr"}
#'   (blue-orange, default), \code{"BuVm"}, \code{"BuPk"}, \code{"GnPk"}.
#'   The color range is computed automatically from the data.
#' @param effect_colors Color function for effect size annotation.
#' @param pvalue_colors Named list of colors for p-value annotation.
#' @param col_higher Character. Color for bars where taxa are higher in the
#'   first condition (default \code{"#E69F00"}, Okabe-Ito orange).
#' @param col_lower Character. Color for bars where taxa are lower in the
#'   first condition (default \code{"#0072B2"}, Okabe-Ito blue).
#'
#' @return A \code{ComplexHeatmap} object (returned invisibly; drawn as a
#'   side effect).
#' @export
#'
#' @examples
#' \dontrun{
#' aldex_heatmap_plot(
#'   table            = feature_table,
#'   metadata         = sample_metadata,
#'   col_cond         = "Sample_type",
#'   effect_threshold = 2
#' )
#' }
#'

aldex_heatmap_plot <- function(table,
                               metadata,
                               col_cond,
                               effect_threshold  = 0.8,
                               pvalue_BH         = NULL,
                               cluster_rows      = FALSE,
                               cluster_columns   = FALSE,
                               diverging_palette = "BuOr",
                               effect_colors = circlize::colorRamp2(
                                 c(-1.5, 0, 1.5),
                                 c("#0072B2", "white", "#E69F00")
                               ),
                               pvalue_colors = list(
                                 "p-value" = c(
                                   "<0.001" = "#0072B2",
                                   "<0.01"  = "#56B4E9",
                                   "<0.05"  = "#E69F00",
                                   ">0.05"  = "grey85"
                                 )
                               ),
                               col_higher = "#E69F00",   # Okabe-Ito orange
                               col_lower  = "#0072B2") { # Okabe-Ito blue

  # Check ComplexHeatmap
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("Package 'ComplexHeatmap' not installed. Installing from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
  }

  # Verify condition column
  if (!col_cond %in% colnames(metadata))
    stop(paste("Column", col_cond, "not found in metadata."))

  conditions       <- metadata[[col_cond]]
  unique_conditions <- unique(conditions)
  if (length(unique_conditions) != 2)
    stop("Exactly two conditions are required for the analysis.")

  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table),
                  ignore.case = TRUE)
  if (length(tax_col) != 1) stop("There is no taxonomy column in the table")
  names(table)[tax_col] <- "taxonomy"

  # Prepare count table
  table        <- table %>% tibble::rownames_to_column(var = "OTUID")
  table_counts <- table %>%
    dplyr::select(-taxonomy) %>%
    tibble::column_to_rownames("OTUID")

  if (length(conditions) != ncol(table_counts))
    stop("Number of conditions does not match number of samples.")

  # Run ALDEx2
  aldex_results <- ALDEx2::aldex(
    reads                  = table_counts,
    conditions             = conditions,
    mc.samples             = 128,
    effect                 = TRUE,
    test                   = "t",
    verbose                = TRUE,
    denom                  = "all",
    include.sample.summary = FALSE
  )

  if (!all(c("effect", "wi.eBH") %in% colnames(aldex_results)))
    stop("Columns 'effect' or 'wi.eBH' missing in ALDEx2 results.")

  # Filter by thresholds
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

  # Prepare plot data
  aldex_plot <- aldex_filtered %>%
    tibble::rownames_to_column("OTUID") %>%
    dplyr::left_join(table %>% dplyr::select(OTUID, taxonomy), by = "OTUID") %>%
    dplyr::mutate(
      seccion = dplyr::case_when(
        diff.btw < 0 ~ paste("Lower in",  unique_conditions[2]),
        diff.btw > 0 ~ paste("Higher in", unique_conditions[1]),
        TRUE         ~ "No Change"
      ),
      taxonomy = dplyr::case_when(
        grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~
          sub(".*g__([^;]*).*", "\\1", taxonomy),
        grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~
          paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
        grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~
          paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
        grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~
          paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
        grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~
          paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
        TRUE ~ "Unclassified"
      ),
      taxonomy = stringr::str_trim(taxonomy),
      taxonomy = make.unique(taxonomy),
      p.value  = dplyr::case_when(
        wi.eBH <= 0.001 ~ "<0.001",
        wi.eBH <= 0.01  ~ "<0.01",
        wi.eBH <  0.05  ~ "<0.05",
        TRUE            ~ ">0.05"
      )
    ) %>%
    dplyr::arrange(diff.btw)

  rab_cols <- paste0("rab.win.", unique_conditions)
  heat_data <- aldex_plot %>%
    dplyr::select(taxonomy, dplyr::all_of(rab_cols)) %>%
    dplyr::rename_with(~ unique_conditions, dplyr::all_of(rab_cols)) %>%
    tibble::column_to_rownames(var = "taxonomy") %>%
    as.matrix()

  # Compute diverging heatmap colors from actual data range
  pal_colors <- .mbm_div_palettes[[diverging_palette]]
  if (is.null(pal_colors)) {
    warning("Unknown diverging_palette '", diverging_palette,
            "'. Using 'BuOr'. Valid options: ",
            paste(names(.mbm_div_palettes), collapse = ", "))
    pal_colors <- .mbm_div_palettes[["BuOr"]]
  }
  data_max <- max(abs(heat_data), na.rm = TRUE)
  heatmap_colors <- circlize::colorRamp2(
    c(-data_max, 0, data_max),
    pal_colors
  )

  # Barplot fill: direct assignment avoids fragile name-lookup
  bar_fills <- ifelse(aldex_plot$diff.btw > 0, col_higher, col_lower)

  # Shared gpar helpers (consistent font/color across all annotations)
  gp_title  <- grid::gpar(fontsize = 12, fontface = "bold",
                           fontfamily = "serif", col = "black")
  gp_labels <- grid::gpar(fontsize = 11, fontfamily = "serif", col = "black")
  gp_border <- grid::gpar(col = "black")

  # --- Left annotation: Effect size ---
  left_annotation <- ComplexHeatmap::rowAnnotation(
    "Effect size" = aldex_plot$effect,
    col           = list("Effect size" = effect_colors),
    simple_anno_size    = grid::unit(0.5, "cm"),
    annotation_name_gp  = gp_title,
    annotation_legend_param = list(
      title_gp  = gp_title,
      labels_gp = gp_labels,
      direction = "vertical"
    ),
    show_legend          = TRUE,
    gp                   = gp_border,
    show_annotation_name = TRUE
  )

  # --- Right annotation 1: p-value tiles ---
  annP <- ComplexHeatmap::rowAnnotation(
    "p-value"        = aldex_plot$p.value,
    simple_anno_size = grid::unit(0.5, "cm"),
    annotation_name_gp = gp_title,
    annotation_legend_param = list(
      title_gp  = gp_title,
      labels_gp = gp_labels,
      direction = "vertical"
    ),
    col          = pvalue_colors,
    show_legend  = TRUE,
    gp           = gp_border,
    show_annotation_name = TRUE
  )

  # --- Right annotation 2: difference barplot ---
  barpl <- ComplexHeatmap::rowAnnotation(
    "difference\nbetween groups" = ComplexHeatmap::anno_barplot(
      aldex_plot$diff.btw,
      which = "row",
      gp    = grid::gpar(fill = bar_fills, col = "black"),
      width = grid::unit(4, "cm")
    ),
    show_annotation_name = TRUE,
    annotation_name_gp   = gp_title,
    annotation_name_rot  = 0
  )

  # --- Main heatmap ---
  heatmap <- ComplexHeatmap::Heatmap(
    heat_data,
    cluster_rows     = cluster_rows,
    cluster_columns  = cluster_columns,
    width            = grid::unit(ncol(heat_data) * 7, "mm"),
    height           = grid::unit(nrow(heat_data) * 6, "mm"),
    column_names_rot = 90,
    rect_gp          = grid::gpar(col = "black", lwd = 1),
    left_annotation  = left_annotation,
    name             = "Median\nclr value",
    heatmap_legend_param = list(
      direction     = "vertical",
      labels_gp     = gp_labels,
      title_gp      = gp_title,
      legend_height = grid::unit(2.5, "cm")
    ),
    column_names_gp  = gp_title,
    col              = heatmap_colors,
    row_names_gp     = grid::gpar(fontsize = 11, fontface = "italic",
                                  fontfamily = "serif", col = "black"),
    show_heatmap_legend = TRUE
  )

  # Draw: main heatmap + p-value annotation + barplot side by side
  ht_list <- heatmap + annP + barpl

  ComplexHeatmap::draw(
    ht_list,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right",
    merge_legend           = FALSE
  )
  return(invisible(ht_list))
}
