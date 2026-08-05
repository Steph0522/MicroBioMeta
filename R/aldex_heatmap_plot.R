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
#' @param heatmap_colors Controls the color scale of the main heatmap body (CLR values).
#'   Three options: \code{NULL} (default, same as \code{"viridis"}) uses a sequential
#'   colorblind-friendly viridis scale (the same family used in
#'   \code{abundance_heatmap_plot}), with range computed automatically from the data;
#'   a preset name string — \code{"viridis"} or one of the diverging presets
#'   \code{"BuOr"} (blue-orange), \code{"BuVm"} (blue-vermillion), \code{"BuPk"}
#'   (blue-pink), \code{"GnPk"} (green-pink); or a \code{circlize::colorRamp2}
#'   function for full manual control.
#' @param effect_colors Color function for the effect size annotation strip
#'   (default: green-white-pink colorblind-friendly scale, distinct from the
#'   heatmap body and p-value defaults).
#' @param pvalue_colors Named list of colors for the p-value annotation strip
#'   (default: black/vermillion/yellow/grey categorical scale, distinct from
#'   the heatmap body and effect size defaults).
#' @param group_colors Optional character vector of colors for the difference
#'   barplot annotation, one color per condition in the order they appear in
#'   \code{metadata[[col_cond]]}. If \code{NULL} (default) a sky blue/yellow
#'   colorblind-friendly palette is used (distinct from the other annotations'
#'   defaults), cycling through the rest of the Okabe-Ito palette as needed
#'   for more than two groups.
#' @param save_table Logical. If \code{TRUE}, saves the underlying ALDEx2
#'   results table (filtered taxa, effect size, diff.btw, p-value category)
#'   to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"aldex_pval_effect.txt"}.
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
                               heatmap_colors    = NULL,
                               effect_colors = circlize::colorRamp2(
                                 c(-1.5, 0, 1.5),
                                 c("#009E73", "white", "#CC79A7")
                               ),
                               pvalue_colors = list(
                                 "p-value" = c(
                                   "<0.001" = "#000000",
                                   "<0.01"  = "#D55E00",
                                   "<0.05"  = "#F0E442",
                                   ">0.05"  = "grey85"
                                 )
                               ),
                               group_colors = NULL,
                               save_table = FALSE,
                               table_filename = "aldex_pval_effect.txt") {

  # Check ComplexHeatmap
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("Package 'ComplexHeatmap' not installed. Installing from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("ComplexHeatmap")
  }

  # Verify condition column
  if (!col_cond %in% colnames(metadata))
    stop(paste("Column", col_cond, "not found in metadata."))

  conditions       <- as.character(metadata[[col_cond]])
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

  if (nrow(aldex_filtered) == 0)
    stop(paste0(
      "No taxa passed the filtering thresholds ",
      "(effect_threshold = ", effect_threshold,
      if (!is.null(pvalue_BH)) paste0(", pvalue_BH = ", pvalue_BH) else "",
      ").\nTry lowering effect_threshold and/or raising pvalue_BH."
    ))

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

  if (save_table) {
    utils::write.table(aldex_plot, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  rab_cols <- paste0("rab.win.", unique_conditions)
  heat_data <- aldex_plot %>%
    dplyr::select(taxonomy, dplyr::all_of(rab_cols)) %>%
    dplyr::rename_with(~ unique_conditions, dplyr::all_of(rab_cols)) %>%
    tibble::column_to_rownames(var = "taxonomy") %>%
    as.matrix()

  # Resolve heatmap body colors
  data_max <- max(abs(heat_data), na.rm = TRUE)
  heatmap_colors_fn <- if (is.null(heatmap_colors) ||
                          (is.character(heatmap_colors) && length(heatmap_colors) == 1 &&
                           tolower(heatmap_colors) == "viridis")) {
    # Default: sequential viridis scale (same family as abundance_heatmap_plot)
    circlize::colorRamp2(
      seq(-data_max, data_max, length.out = 13),
      viridis::viridis(13, option = "C", direction = -1)
    )
  } else if (is.character(heatmap_colors) && length(heatmap_colors) == 1) {
    # Diverging preset name string
    pal <- .mbm_div_palettes[[heatmap_colors]]
    if (is.null(pal)) {
      warning("Unknown heatmap_colors preset '", heatmap_colors,
              "'. Using 'viridis'. Valid options: 'viridis', ",
              paste(names(.mbm_div_palettes), collapse = ", "))
      pal <- NULL
    }
    if (is.null(pal)) {
      circlize::colorRamp2(
        seq(-data_max, data_max, length.out = 13),
        viridis::viridis(13, option = "C", direction = -1)
      )
    } else {
      circlize::colorRamp2(c(-data_max, 0, data_max), pal)
    }
  } else {
    heatmap_colors   # assume already a colorRamp2 function
  }

  # Resolve group colors for barplot: default starts at sky blue/yellow
  # (distinct from the blue/orange heatmap and green/pink effect size
  # defaults), then cycles through the rest of Okabe-Ito for N conditions
  if (is.null(group_colors)) {
    group_default <- .mbm_colors[c(2, 4, 3, 7, 6, 8, 5, 1)]
    group_colors  <- rep_len(group_default, length(unique_conditions))
  } else {
    group_colors <- rep_len(group_colors, length(unique_conditions))
  }
  # Map: diff.btw > 0 → condition 1 color; diff.btw < 0 → condition 2 color
  bar_fills <- ifelse(aldex_plot$diff.btw > 0, group_colors[1], group_colors[2])

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
    col              = heatmap_colors_fn,
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
