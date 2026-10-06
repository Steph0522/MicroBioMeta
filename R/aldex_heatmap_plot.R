#' ALDEx2 differential abundance heatmap
#'
#' Runs ALDEx2 on a table and returns a ComplexHeatmap
#'
#' @param table A data frame with taxa in rows and samples in columns. 
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. 
#' Must include a `SAMPLEID` column matching sample names in `table`.
#' @param group_col Character. Name of the column in \code{metadata} that
#'   defines the two groups to compare. Exactly two unique values are required.
#' @param effect_threshold Numeric. Minimum absolute effect size to retain.
#'   Default \code{0} (no effect-size filter), so by default the heatmap shows
#'   the same taxa that \code{aldex_volcano_plot} colors as significant.
#' @param pval_threshold Numeric or NULL. Maximum p-value to retain (adjusted
#'   or not, depending on \code{p_adjust_method}). Default \code{0.05}. If
#'   NULL only \code{effect_threshold} is applied.
#' @param p_adjust_method \code{"BH"} (default) or \code{"none"}: whether
#'   \code{pval_threshold} and the p-value annotation use ALDEx2's
#'   Benjamini-Hochberg adjusted p-values (\code{wi.eBH}) or the raw ones
#'   (\code{wi.ep}). ALDEx2 only computes the BH correction.
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
#' @param effect_colors Three colors for the effect size annotation strip:
#'   negative, zero and positive effect (mapped to -1.5, 0 and 1.5). Default
#'   \code{c("#009E73", "white", "#CC79A7")} (green-white-pink,
#'   colorblind-friendly, distinct from the heatmap body and p-value
#'   defaults). A \code{circlize::colorRamp2()} function is also accepted.
#' @param pvalue_colors Named vector of colors for the p-value annotation
#'   strip, with names \code{"<0.001"}, \code{"<0.01"}, \code{"<0.05"} and
#'   \code{">0.05"}. Default black/vermillion/yellow/grey (distinct from the
#'   heatmap body and effect size defaults). A list with one such vector named
#'   \code{"p-value"} is also accepted.
#' @param group_colors Optional character vector of colors for the difference
#'   barplot annotation, either named after the conditions (e.g.
#'   \code{c(Rhizosphere = "#56B4E9", Roots = "#009E73")}) or one color per
#'   condition in the order they appear in
#'   \code{metadata[[group_col]]}. If \code{NULL} (default) an orange/blue
#'   colorblind-friendly palette is used (matching the col_sup/col_inf
#'   convention used elsewhere, e.g. \code{aldex_volcano_plot}), cycling
#'   through the rest of the Okabe-Ito palette as needed
#'   for more than two groups.
#' @param save_table Logical. If \code{TRUE}, saves the underlying ALDEx2
#'   results table (filtered taxa, effect size, diff.btw, p-value category)
#'   to disk. Default \code{FALSE}. \code{effect} and \code{diff.btw} are
#'   saved as ALDEx2 returns them (alphabetically second group minus the
#'   first); the \code{seccion} column says in which group each taxon is higher.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"aldex_pval_effect.txt"}.
#'
#' @param draw Logical. If \code{TRUE} (default), the heatmap is drawn on the
#'   current device. Use \code{FALSE} to only build the returned grob without
#'   drawing it.
#' @return Invisibly, a \code{gTree} (grid grob) with the heatmap. Printing it
#'   (e.g. typing its name) draws the heatmap; it can also be combined with
#'   other plots (e.g. \code{cowplot::plot_grid()},
#'   \code{patchwork::wrap_elements()}).
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' # group_col must have exactly two groups; Location has two
#' # (Rhizosphere and Roots) in the bundled example data. Here taxa are
#' # filtered by effect size alone (pval_threshold = NULL), since this small
#' # (46-sample) dataset has few taxa with significant p-values.
#' aldex_heatmap_plot(
#'   table            = table,
#'   metadata         = metadata,
#'   group_col        = "Location",
#'   effect_threshold = 0.5,
#'   pval_threshold   = NULL
#' )
#'
#' # All the colors have colorblind-friendly defaults, but each can be set by
#' # hand: group_colors (difference bars, named after the groups),
#' # effect_colors (negative, zero and positive effect size), pvalue_colors
#' # (one color per p-value class) and heatmap_colors (median clr values)
#' \donttest{
#' aldex_heatmap_plot(
#'   table            = table,
#'   metadata         = metadata,
#'   group_col        = "Location",
#'   effect_threshold = 0.5,
#'   pval_threshold   = NULL,
#'   group_colors     = c(Rhizosphere = "#56B4E9", Roots = "#009E73"),
#'   effect_colors    = c("#0072B2", "white", "#E69F00"),
#'   pvalue_colors    = c("<0.001" = "black", "<0.01" = "grey30",
#'                        "<0.05" = "grey60", ">0.05" = "grey90"),
#'   heatmap_colors   = "BuOr"
#' )
#' }
#'

aldex_heatmap_plot <- function(table,
                               metadata,
                               group_col,
                               effect_threshold  = 0,
                               pval_threshold    = 0.05,
                               p_adjust_method   = "BH",
                               cluster_rows      = FALSE,
                               cluster_columns   = FALSE,
                               heatmap_colors    = NULL,
                               effect_colors = c("#009E73", "white", "#CC79A7"),
                               pvalue_colors = c(
                                 "<0.001" = "#000000",
                                 "<0.01"  = "#D55E00",
                                 "<0.05"  = "#F0E442",
                                 ">0.05"  = "grey85"
                               ),
                               group_colors = NULL,
                               save_table = FALSE,
                               table_filename = "aldex_pval_effect.txt",
                               draw = TRUE) {

  # Check ComplexHeatmap
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop(
      "Package 'ComplexHeatmap' is required but not installed.\n",
      "Install it with: BiocManager::install(\"ComplexHeatmap\")",
      call. = FALSE
    )
  }

  # verify condition column
  if (!group_col %in% colnames(metadata))
    stop("Column ", group_col, " not found in metadata.")

  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table),
                  ignore.case = TRUE)
  if (length(tax_col) != 1) stop("There is no taxonomy column in the table")
  names(table)[tax_col] <- "taxonomy"

  # prepare count table
  table        <- table %>% tibble::rownames_to_column(var = "OTUID")
  table_counts <- table %>%
    dplyr::select(-taxonomy) %>%
    tibble::column_to_rownames("OTUID")

  
  metadata     <- .mbm_align_metadata(colnames(table_counts), metadata)
  table_counts <- table_counts[, metadata[[1]], drop = FALSE]

  conditions       <- as.character(metadata[[group_col]])
  unique_conditions <- unique(conditions)
  if (length(unique_conditions) != 2)
    stop("Exactly two conditions are required for the analysis.")

  if (length(conditions) != ncol(table_counts))
    stop("Number of conditions does not match number of samples.")

  # run ALDEx2
  aldex_results <- ALDEx2::aldex(
    reads                  = table_counts,
    conditions             = conditions,
    mc.samples             = 128,
    effect                 = TRUE,
    test                   = "t",
    verbose                = FALSE,
    denom                  = "all",
    include.sample.summary = FALSE
  )

  if (!all(c("effect", "wi.eBH") %in% colnames(aldex_results)))
    stop("Columns 'effect' or 'wi.eBH' missing in ALDEx2 results.")

  
  flip_sign <- unique_conditions[1] == sort(unique_conditions)[1]
  if (flip_sign) {
    aldex_results$diff.btw <- -aldex_results$diff.btw
    aldex_results$effect   <- -aldex_results$effect
  }

  p_adjust_method <- match.arg(p_adjust_method, c("BH", "none"))
  aldex_results$.p <- if (p_adjust_method == "BH") aldex_results$wi.eBH else aldex_results$wi.ep
  aldex_filtered <- aldex_results
  if (effect_threshold > 0 && !is.null(pval_threshold)) {
    aldex_filtered <- aldex_results %>%
      dplyr::filter(abs(effect) >= effect_threshold, .p <= pval_threshold)
  } else if (effect_threshold > 0) {
    aldex_filtered <- aldex_results %>%
      dplyr::filter(abs(effect) >= effect_threshold)
  } else if (!is.null(pval_threshold)) {
    aldex_filtered <- aldex_results %>%
      dplyr::filter(.p <= pval_threshold)
  }

  if (nrow(aldex_filtered) == 0)
    stop(
      "No taxa passed the filtering thresholds ",
      "(effect_threshold = ", effect_threshold,
      if (!is.null(pval_threshold)) ", pval_threshold = " else "",
      if (!is.null(pval_threshold)) pval_threshold else "",
      ").\nTry lowering effect_threshold and/or raising pval_threshold."
    )

  # Prepare plot data
  aldex_plot <- aldex_filtered %>%
    tibble::rownames_to_column("OTUID") %>%
    dplyr::left_join(table %>% dplyr::select(OTUID, taxonomy), by = "OTUID") %>%
    dplyr::mutate(
      seccion = dplyr::case_when(
        diff.btw < 0 ~ paste("Higher in", unique_conditions[2]),
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
      taxonomy = .mbm_composite_genus(stringr::str_trim(taxonomy)),
      taxonomy = make.unique(taxonomy),
      p.value  = dplyr::case_when(
        .p <= 0.001 ~ "<0.001",
        .p <= 0.01  ~ "<0.01",
        .p <  0.05  ~ "<0.05",
        TRUE            ~ ">0.05"
      )
    ) %>%
    dplyr::arrange(diff.btw)

  if (save_table) {
    aldex_saved <- aldex_plot
    aldex_saved$.p <- NULL
    if (flip_sign) {
      aldex_saved$diff.btw <- -aldex_saved$diff.btw
      aldex_saved$effect   <- -aldex_saved$effect
    }
    utils::write.table(aldex_saved, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message("Table saved as: ", table_filename)
  }

  rab_cols <- paste0("rab.win.", unique_conditions)
  heat_data <- aldex_plot %>%
    dplyr::select(taxonomy, dplyr::all_of(rab_cols)) %>%
    dplyr::rename_with(~ unique_conditions, dplyr::all_of(rab_cols)) %>%
    tibble::column_to_rownames(var = "taxonomy") %>%
    as.matrix()

  data_max <- max(abs(heat_data), na.rm = TRUE)
  heatmap_colors_fn <- if (is.null(heatmap_colors) ||
                          (is.character(heatmap_colors) && length(heatmap_colors) == 1 &&
                           tolower(heatmap_colors) == "viridis")) {
    circlize::colorRamp2(
      seq(-data_max, data_max, length.out = 13),
      viridis::viridis(13, option = "C", direction = -1)
    )
  } else if (is.character(heatmap_colors) && length(heatmap_colors) == 1) {
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
    heatmap_colors   
  }


  if (is.null(group_colors)) {
    group_default <- .mbm_colors[c(1, 5, 3, 7, 6, 8, 2, 4)]
    group_colors  <- rep_len(group_default, length(unique_conditions))
  } else if (!is.null(names(group_colors))) {
    # Named colors are matched to the conditions by name
    group_colors <- .mbm_match_colors(group_colors, as.character(unique_conditions),
                                      "group_colors")
  } else {
    group_colors <- rep_len(group_colors, length(unique_conditions))
  }
  bar_fills <- ifelse(aldex_plot$diff.btw > 0, group_colors[1], group_colors[2])

  gp_title  <- grid::gpar(fontsize = 12, fontface = "bold",
                           fontfamily = "serif", col = "black")
  gp_legend_title <- grid::gpar(fontsize = 14, fontface = "bold",
                                fontfamily = "serif", col = "black")
  gp_labels <- grid::gpar(fontsize = 12, fontfamily = "serif", col = "black")
 
  gp_border <- grid::gpar(col = "white", lwd = 1.5)

  
  if (is.character(effect_colors)) {
    if (length(effect_colors) != 3)
      stop("`effect_colors` must be 3 colors (negative, zero, positive effect).")
    effect_colors <- circlize::colorRamp2(c(-1.5, 0, 1.5), effect_colors)
  }
  if (!is.list(pvalue_colors)) pvalue_colors <- list("p-value" = pvalue_colors)

  left_annotation <- ComplexHeatmap::rowAnnotation(
    "Effect size" = aldex_plot$effect,
    col           = list("Effect size" = effect_colors),
    simple_anno_size    = grid::unit(0.5, "cm"),
    annotation_name_gp  = gp_title,
    annotation_legend_param = list(
      title_gp  = gp_legend_title,
      labels_gp = gp_labels,
      direction = "vertical"
    ),
    show_legend          = TRUE,
    gp                   = gp_border,
    show_annotation_name = TRUE
  )

  annP <- ComplexHeatmap::rowAnnotation(
    "p-value"        = aldex_plot$p.value,
    simple_anno_size = grid::unit(0.5, "cm"),
    annotation_name_gp = gp_title,
    annotation_legend_param = list(
      title_gp  = gp_legend_title,
      labels_gp = gp_labels,
      direction = "vertical"
    ),
    col          = pvalue_colors,
    show_legend  = TRUE,
    gp           = gp_border,
    show_annotation_name = TRUE
  )

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

  heatmap <- ComplexHeatmap::Heatmap(
    heat_data,
    cluster_rows     = cluster_rows,
    cluster_columns  = cluster_columns,
    width            = grid::unit(ncol(heat_data) * 7, "mm"),
    height           = grid::unit(nrow(heat_data) * 8, "mm"),
    column_names_rot = 90,
    rect_gp          = grid::gpar(col = "white", lwd = 1.5),
    left_annotation  = left_annotation,
    name             = "Median\nclr value",
    heatmap_legend_param = list(
      direction     = "vertical",
      labels_gp     = gp_labels,
      title_gp      = gp_legend_title,
      legend_height = grid::unit(2.5, "cm")
    ),
    column_names_gp  = gp_title,
    col              = heatmap_colors_fn,
    
    show_row_names   = FALSE,
    show_heatmap_legend = TRUE
  )


  taxon_names <- sub("[.][0-9]+$", "", rownames(heat_data))
  taxon_face  <- ifelse(grepl("^other |^Unclassified", taxon_names), "plain", "italic")

  
  taxon_labels <- ComplexHeatmap::rowAnnotation(
    taxon = ComplexHeatmap::anno_text(
      taxon_names,
      gp   = grid::gpar(fontsize = 12, fontfamily = "serif", col = "black",
                        fontface = taxon_face),
      just = "left"
    ),
    show_annotation_name = FALSE
  )

 
  ht_list <- heatmap + annP + barpl + taxon_labels

  draw_heatmap <- function() ComplexHeatmap::draw(
    ht_list,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right",
    merge_legend           = FALSE,
    
    padding = grid::unit(c(2, 2, 2, 12), "mm"),
    background = "transparent"
  )

  
  if (draw) {
    grid::grid.newpage()
    draw_heatmap()
    heatmap_output <- grid::grid.grab(wrap.grobs = TRUE)
  } else {
    heatmap_output <- grid::grid.grabExpr(draw_heatmap())
  }
  class(heatmap_output) <- c("mbm_heatmap", class(heatmap_output))
  return(invisible(heatmap_output))
}
