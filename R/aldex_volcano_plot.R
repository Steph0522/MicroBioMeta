#' This function generates either an effect size plot or a volcano plot based on ALDEx2 results.
#'
#' @param table Data frame with count data; columns represent samples, rows represent features.
#' @param metadata Data frame containing metadata for the samples.
#' @param group_col Name of the column in `metadata` that contains the experimental conditions.
#' @param type Type of plot to generate: `"volcano"` (default, volcano plot)
#'   or `"effect"` (effect-size plot). Case-insensitive.
#' @param col_inf Color for points lower than threshold. Default `'#56B4E9'` (Okabe-Ito blue, matching the package's 2-group default).
#' @param col_sup Color for points higher than threshold. Default `'#E69F00'` (Okabe-Ito orange).
#' @param threshold_lower Lower threshold for effect size/difference (x-axis).
#' @param threshold_upper Upper threshold for effect size/difference (x-axis).
#' @param cond Name of the reference condition for the plot labels: positive
#'   values on the x-axis mean higher in \code{cond}. Default \code{NULL}, the
#'   condition of the first sample in \code{table}.
#' @param pval_threshold p-value cutoff for significance (default = 0.05), drawn
#'   as the dashed horizontal line.
#' @param p_adjust_method \code{"BH"} (default) or \code{"none"}. With
#'   \code{"BH"}, the y-axis and the significance cutoff use ALDEx2's
#'   Benjamini-Hochberg adjusted p-values (\code{wi.eBH}), since thousands of
#'   taxa are tested at once; \code{"none"} uses the raw p-values
#'   (\code{wi.ep}). ALDEx2 only computes the BH correction, so no other
#'   method is available here.
#' @param show_labels Logical. Whether to display "Higher/Lower in cond" labels (for "effect" plot only, default is TRUE).
#' @param taxa Data frame with taxonomic information (required for "volcano" plot only).
#' @param x_axis_title,y_axis_title Titles for the x- and y-axis. Default
#'   \code{NULL}: \code{"Effect size"} (\code{"effect"}) or log2 fold change
#'   (\code{"volcano"}) on x, and the (adjusted) p-value on y.
#' @param label_size Numeric. Font size of the taxon labels drawn on
#'   significant points in the "volcano" plot. Default \code{3.5}.
#' @param filter_uncultured Logical. If \code{TRUE}, taxa whose name contains
#'   "uncultured"/"unculture" are still plotted as points but not labelled,
#'   keeping the volcano plot's text annotations to named taxa. Default
#'   \code{FALSE}.
#' @param save_table Logical. If \code{TRUE}, saves the ALDEx2 result table to disk. Default \code{FALSE}.
#'   \code{effect} and \code{diff.btw} are saved as ALDEx2 returns them
#'   (alphabetically second group minus the first); in the plot they are
#'   shown so that positive means higher in \code{cond}.
#' @param table_filename Character. File path/name for the saved table. Default \code{"aldex_pval_effect.txt"}.
#'
#' @param ... Old names of renamed arguments (\code{col_cond}, \code{cutoff.pval}, \code{adjusted_p}), still accepted
#'   with a warning. Any other extra argument is an error.
#' @return A `ggplot` object with the selected plot.
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
#' # (Rhizosphere and Roots) in the bundled example data
#' aldex_volcano_plot(
#'   table           = table,
#'   metadata        = metadata,
#'   group_col        = "Location",
#'   type            = "effect",
#'   col_inf         = "#56B4E9",
#'   col_sup         = "#E69F00",
#'   threshold_lower = -0.5,
#'   threshold_upper = 0.5,
#'   cond            = "Rhizosphere",
#'   show_labels     = TRUE
#' )

aldex_volcano_plot <- function(table,
                               metadata,
                               group_col,
                               type = "volcano",
                               col_inf = "#56B4E9",
                               col_sup = "#E69F00",
                               threshold_lower = -1.5,
                               threshold_upper = 1.5,
                               cond = NULL,
                               pval_threshold = 0.05,
                               p_adjust_method = "BH",
                               show_labels = TRUE,
                               taxa = NULL,
                               x_axis_title = NULL,
                               y_axis_title = NULL,
                               label_size = 3.5,
                               filter_uncultured = FALSE,
                               save_table = FALSE,
                               table_filename = "aldex_pval_effect.txt",
                               ...) {
  # Old argument names still work, with a warning (see .mbm_renamed_args)
  renamed <- .mbm_renamed_args(list(...), c(col_cond = "group_col", cutoff.pval = "pval_threshold", adjusted_p = "p_adjust_method"), "aldex_volcano_plot")
  for (nm in names(renamed)) assign(nm, renamed[[nm]])
  # the old adjusted_p was TRUE/FALSE
  if (is.logical(p_adjust_method)) p_adjust_method <- if (isTRUE(p_adjust_method)) "BH" else "none"

  # Verify that type has a valid value (case-insensitive)
  type <- tolower(type)
  if (!type %in% c("effect", "volcano")) {
    stop("`type` must be either \"effect\" or \"volcano\".", call. = FALSE)
  }
  
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("Package 'ALDEx2' needed for this function to work. Please install it.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.")
  }
  
  # Verify that the condition column exists
  if (!group_col %in% colnames(metadata)) {
    stop(
      "The column ", group_col,
      " does not exist in the object 'metadata'. Check the name is written correctly."
    )
  }
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  
  #Automatically identify the taxonomy column (last column)
  if (is.null(taxa)) {
    last_col <- ncol(table)
   taxa_colname <- colnames(table)[last_col]
   taxa <- data.frame(
    Feature.ID = rownames(table),
   Taxon = table[[taxa_colname]],
   stringsAsFactors = FALSE
  )
  #  Remove the last column of table for the analysis
    table <- table[, -last_col, drop = FALSE]
  }

  # Align samples between table and metadata (first column = sample ID,
  # regardless of its original name)
  metadata <- .mbm_align_metadata(colnames(table), metadata)
  table    <- table[, metadata[[1]], drop = FALSE]

  conditions <- as.character(metadata[[group_col]])
  groups <- unique(conditions)
  if (is.null(cond)) cond <- groups[1]
  if (!cond %in% groups)
    stop("`cond` must be one of: ", paste(groups, collapse = ", "), call. = FALSE)
  other_cond <- setdiff(groups, cond)[1]
  
  aldex_clr <- ALDEx2::aldex(table, conditions, mc.samples = 128, denom = "all")
  
  
  
  
  
  # Process taxonomic data for both plot types
  processed_data <- aldex_clr %>%
    tibble::rownames_to_column(var = "Feature.ID") %>%
    dplyr::left_join(taxa, by = "Feature.ID") %>%
    dplyr::mutate(
      taxa = dplyr::case_when(
        stringr::str_detect(Taxon, "g__") ~ stringr::str_extract(Taxon, "(?<=g__)[^_;]+"),
        stringr::str_detect(Taxon, "f__") ~ stringr::str_extract(Taxon, "(?<=f__)[^_;]+"),
        stringr::str_detect(Taxon, "c__") ~ stringr::str_extract(Taxon, "(?<=c__)[^_;]+"),
        stringr::str_detect(Taxon, "o__") ~ stringr::str_extract(Taxon, "(?<=o__)[^_;]+"),
        TRUE ~ Feature.ID 
      ),
      taxa = .mbm_composite_genus(taxa)
    )
  
  # Guardar tabla si se solicita
  if (save_table) {
    utils::write.table(
      processed_data,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    ) 
    message("Table saved as: ", table_filename)
  }

  # The table above is saved exactly as ALDEx2 returns it. ALDEx2 computes
  # diff.btw/effect as the alphabetically second group minus the first, so
  # for the plot flip the sign when needed: a positive value then always
  # means higher in `cond`, which is what the "Higher/Lower in" labels assume.
  if (cond == sort(groups)[1]) {
    processed_data$diff.btw <- -processed_data$diff.btw
    processed_data$effect   <- -processed_data$effect
  }
  
  
  # p-value used for the y-axis and the significance cutoff
  p_adjust_method <- match.arg(p_adjust_method, c("BH", "none"))
  adjusted <- p_adjust_method == "BH"
  processed_data$.p <- if (adjusted) processed_data$wi.eBH else processed_data$wi.ep
  y_lab <- if (adjusted) {
    expression("-Log"[10]~"adjusted p-value (BH)")
  } else {
    expression("-Log"[10]~"p-value")
  }

  if (type == "effect") {
    # Preparar datos para effect plot
    plot_data <- processed_data %>%
      dplyr::mutate(
        grupo = dplyr::case_when(
          effect <= threshold_lower ~ paste("Lower in", cond),
          effect >= threshold_upper ~ paste("Higher in", cond),
          TRUE ~ "Not significant"
        ),
        log_pvalue = -log10(.p + min(.p[.p > 0])/10)
      )
    
    # Obtener los top taxones para etiquetar
    top_taxa <- plot_data %>%
      dplyr::filter(grupo != "Not significant") %>%
      dplyr::group_by(grupo) %>%
      dplyr::slice_max(order_by = abs(effect), n = 3)
    
    lim_x <- max(abs(plot_data$effect), na.rm = TRUE)
    
    # Create a color vector with dynamic names
    color_values <- c(col_sup, col_inf, "gray")
    names(color_values) <- c(paste("Higher in", cond), paste("Lower in", cond), "Not significant")
    
    # Create the base plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect, y = log_pvalue, color = grupo)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::geom_vline(
        xintercept = c(threshold_lower, threshold_upper),
        linetype = 2,
        color = "black"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(pval_threshold),
        linetype = 2,
        color = "black"
      ) +
      ggplot2::labs(
        x = if (is.null(x_axis_title)) "Effect size" else x_axis_title,
        y = if (is.null(y_axis_title)) y_lab else y_axis_title,
        color = NULL
      ) +
      .mbm_theme(legend_position = "none") +
      ggplot2::scale_x_continuous(limits = c(-lim_x, lim_x)) +
      # Generous top headroom: the corner condition labels sit in this
      # padding, strictly above the highest data point, so they don't
      # compete with that point's own taxon-name label for space.
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30)))

    # Add labels for significant taxa
    if (nrow(top_taxa) > 0) {
      plot_taxa <- if (filter_uncultured) {
        top_taxa[!grepl("uncultured|unculture", top_taxa$taxa, ignore.case = TRUE), ]
      } else top_taxa
      if (nrow(plot_taxa) > 0) {
        p <- p +
          # repelled so nearby taxa don't overlap or get cut at the edges
          ggrepel::geom_text_repel(
            data = plot_taxa,
            ggplot2::aes(label = taxa),
            color = "black",
            family = "serif",
            size = label_size,
            fontface = "italic",
            box.padding = 0.5,
            min.segment.length = 0.3,
            segment.color = "grey50",
            max.overlaps = Inf,
            seed = 1
          )
      }
    }

    # Condition labels: anchored to the panel's own top-left/top-right
    # corners (-Inf/Inf on x, Inf on y) with exactly hjust = 0/1 (the only
    # values that keep text fully inside an Inf-anchored edge without
    # clipping - an inset like 0.05/0.95 still overhangs past the boundary
    # on the outward side). vjust > 1 plus the generous top expansion above
    # push them into the padding strictly above the highest data point, so
    # they don't compete with that point's own taxon-name label for space.
    if (show_labels) {
      p <- p +
        ggplot2::annotate(
          "text",
          x = -Inf,
          y = Inf,
          label = paste0("Lower in\n", cond),
          color = col_inf,
          family = "serif",
          fontface = "bold",
          size = 5,
          hjust = 0,
          vjust = 1.3
        ) +
        ggplot2::annotate(
          "text",
          x = Inf,
          y = Inf,
          label = paste0("Higher in\n", cond),
          color = col_sup,
          family = "serif",
          fontface = "bold",
          size = 5,
          hjust = 1,
          vjust = 1.3
        )
    }
    
  } else if (type == "volcano") {
    # Preparar datos para volcano plot
    plot_data <- processed_data %>%
      dplyr::mutate(
        log_pvalue = -log10(.p + min(.p[.p > 0])/10),
        significant = .p <= pval_threshold,
        direction = ifelse(diff.btw < 0, 
                           paste("Lower in", cond),
                           paste("Higher in", cond))
      )
    
    # Obtener los top taxones para etiquetar
    top_taxa <- plot_data %>%
      dplyr::filter(significant) %>%
      dplyr::group_by(direction) %>%
      dplyr::slice_max(order_by = abs(diff.btw), n = 3)
    
    # Create a color vector with dynamic names
    color_values <- c(col_sup, col_inf)
    names(color_values) <- c(paste("Higher in", cond), paste("Lower in", cond))
    
    # Create the base plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = diff.btw, y = log_pvalue)) +
      ggplot2::geom_point(
        data = dplyr::filter(plot_data, !significant),
        color = "gray",
        size = 3,
        alpha = 0.7
      ) +
      ggplot2::geom_point(
        data = dplyr::filter(plot_data, significant),
        ggplot2::aes(color = direction),
        size = 3
      ) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::geom_vline(
        xintercept = c(threshold_lower, threshold_upper),
        color = 'black',
        linetype = 'dashed'
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(pval_threshold),
        color = 'black',
        linetype = 'dashed'
      ) +
      ggplot2::labs(
        x = if (is.null(x_axis_title)) expression("Log"[2]~"Fold Change") else x_axis_title,
        y = if (is.null(y_axis_title)) y_lab else y_axis_title,
        color = NULL
      ) +
      .mbm_theme(legend_position = "none") +
      # Generous top headroom: the corner condition labels sit in this
      # padding, strictly above the highest data point, so they don't
      # compete with that point's own taxon-name label for space.
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.1)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30)))

    # Add taxa labels
    if (nrow(top_taxa) > 0) {
      plot_taxa <- if (filter_uncultured) {
        top_taxa[!grepl("uncultured|unculture", top_taxa$taxa, ignore.case = TRUE), ]
      } else top_taxa
      if (nrow(plot_taxa) > 0) {
        p <- p +
          # repelled so nearby taxa don't overlap or get cut at the edges
          ggrepel::geom_text_repel(
            data = plot_taxa,
            ggplot2::aes(label = taxa),
            color = "black",
            family = "serif",
            size = label_size,
            fontface = "italic",
            box.padding = 0.5,
            min.segment.length = 0.3,
            segment.color = "grey50",
            max.overlaps = Inf,
            seed = 1
          )
      }
    }

    # Condition labels: anchored to the panel's own top-left/top-right
    # corners (see the "effect" branch above for why: exact hjust = 0/1, and
    # the generous top expansion above keeps them clear of the highest
    # point's own taxon-name label).
    p <- p +
      ggplot2::annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = paste0("Lower in\n", cond),
        color = col_inf,
        family = "serif",
        fontface = "bold",
        size = 5,
        hjust = 0,
        vjust = 1.3
      ) +
      ggplot2::annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0("Higher in\n", cond),
        color = col_sup,
        family = "serif",
        fontface = "bold",
        size = 5,
        hjust = 1,
        vjust = 1.3
      )
  }
  
  return(p)
}
