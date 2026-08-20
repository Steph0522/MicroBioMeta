#' Alpha diversity plot
#'
#' This function generates a boxplot or barplot to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows.
#'              The first column must contain the OTUID, ASV, or species name.
#' @param metadata A data frame with metadata. The first column must match sample names in `table`.
#' @param type Type of plot: either "boxplot" or "barplot". Default is "boxplot".
#' @param stat Optional. A string indicating the test used for comparing means (e.g., "wilcox.test").
#' @param x_col Column in `metadata` to be used on the x-axis.
#' @param fill_col Column in `metadata` to define fill color.
#' @param facet_by Optional. A metadata column to facet (e.g., Treatment, Site).
#' @param facet_by2 Optional. A metadata column to double facet (e.g., Treatment, Site).
#' @param facet_orientation Whether `facet_by` appears in columns ("horizontal", default) or rows ("vertical").
#' @param palette Color palette to use: "colorb", "grey", "viridis", or "brewer". Default: "colorb".
#' @param group_colors A vector of custom colors. Overrides `palette` if provided.
#' @param n_cols Number of columns in facet wrap (optional).
#' @param n_rows Number of rows in facet wrap (optional).
#' @param strip_color Background color of facet strips. Default: "grey".
#' @param show_legend Logical. Show legend? Default: TRUE.
#' @param legend_title Title for the legend.
#' @param legend_position Position of the legend: "bottom", "top", "right", or "left". Default is "bottom".
#' @param title Title for the entire plot.
#' @param x_axis_title Title for the x-axis.
#' @param y_axis_title Title for the y-axis.
#' @param free_y Logical. Whether y-axis scales are free across facets. Default \code{FALSE}.
#' @param rarefy_depth Integer. If provided, rarefies samples to this depth before diversity estimation. Default \code{NULL}.
#' @param panel_label_case Character. Case of the auto-generated panel tags
#'   (A, B, C... added per panel when \code{stat} is used). One of
#'   \code{"upper"} (default, "A", "B", "C") or \code{"lower"} ("a", "b", "c").
#'   Ignored if \code{panel_labels} is supplied.
#' @param panel_labels Optional character vector of custom panel tags, one per
#'   panel, used as-is (e.g. \code{c("(a)", "(b)", "(c)")} or
#'   \code{c("a.", "b.", "c.")}) — for journal styles that
#'   \code{panel_label_case} alone can't produce. Overrides
#'   \code{panel_label_case} when provided.
#' @param panel_label_bold Logical. If \code{TRUE} (default), panel tags are
#'   bold. Set to \code{FALSE} for journals that require plain (non-bold)
#'   panel tags.
#' @param x_label_angle Numeric. Rotation (in degrees) of the x-axis tick
#'   labels. Default \code{0} (horizontal); use e.g. \code{45} or \code{90} when
#'   group names are long enough to overlap.
#' @param strip_text_bold Logical. If \code{TRUE}, facet strip labels are bold.
#'   Default \code{FALSE} (plain).
#' @param aspect_ratio Numeric. Aspect ratio (height/width) of each panel.
#'   Default \code{NULL}, which lets the panels fill the available space.
#' @param save_table Logical. If \code{TRUE}, saves the diversity table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table. Default \code{"diversity.txt"}.
#'
#' @return A ggplot object showing alpha diversity with Hill numbers.
#' @export
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SAMPLEID"
#'
#' # facet_by + stat triggers the auto A/B/C panel tags (see panel_label_case
#' # and panel_labels to customize their case/format)
#' alpha_diversity_plot(
#'   table            = table,
#'   metadata         = metadata,
#'   type             = "barplot",
#'   x_col            = "Location",
#'   fill_col         = "Location",
#'   facet_by         = "Treatment",
#'   free_y           = TRUE,
#'   legend_position  = "top",
#'   stat             = "t.test",
#'   panel_label_case = "upper"
#' )
#' }

alpha_diversity_plot <- function(
    table,
    metadata,
    type = "boxplot",
    stat = NULL,
    x_col,
    fill_col,
    facet_by = NULL,
    facet_by2 = NULL,  # New parameter
    facet_orientation = "horizontal",
    palette = "colorb",
    group_colors = NULL,
    n_cols = NULL,
    n_rows = NULL,
    strip_color = "grey",
    show_legend = TRUE,
    title = NULL,
    legend_title = NULL,
    legend_position = "bottom",
    x_axis_title = NULL,
    y_axis_title = "Diversity measure",
    free_y = FALSE,
    rarefy_depth=NULL,
    panel_label_case = "upper",
    panel_labels = NULL,
    panel_label_bold = TRUE,
    x_label_angle = 0,
    strip_text_bold = FALSE,
    aspect_ratio = NULL,
    save_table = FALSE,
    table_filename = "diversity.txt") {

  # Treat metadata's first column as the sample ID regardless of its original name
  colnames(metadata)[1] <- "SAMPLEID"

  #this part may not be necessary; not sure if it's a good idea to remove ASVs
  #without a taxonomic ID, since that doesn't mean they shouldn't be considered
  #table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
  # Remove uninformative taxonomy strings
  #table <- table %>%
   # dplyr::filter(taxonomy != "d__Bacteria;__;__;__;__;__")
  
  # Keep only samples present in both table and metadata. Final row order
  # doesn't matter here: per-sample diversity indices are computed
  # independently, and results are joined back to metadata by SAMPLEID below.
  sample_order <- metadata[[1]]
  common_samples <- intersect(colnames(table), sample_order)
  if (length(common_samples) == 0) stop("No matching sample names between table and metadata.")
  table <- table[, common_samples, drop = FALSE]
  table <- data.frame(t(table))
  
  if (!is.null(rarefy_depth) && rarefy_depth > 0) {
    table <- vegan::rrarefy(table, sample = rarefy_depth)
  }
  
  # Estimate diversity indices
  est_richness <- vegan::estimateR(table)  # samples as rows
  chao1 <- est_richness["S.chao1", ]
  shannon <- vegan::diversity(table, index = "shannon")
  simpson <- vegan::diversity(table, index = "simpson")

  # Build unified data frame
  results <- data.frame(
    SAMPLEID = rownames(table),
    Chao1 = as.numeric(chao1),
    Shannon = shannon,
    Simpson = simpson
  )
  
  # Unir con metadata por SAMPLEID
  results <- dplyr::left_join(results, metadata, by = "SAMPLEID")
  
  # Guardar tabla si se solicita
  if (save_table) {
    utils::write.table(
      results,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message(paste("Table saved as:", table_filename))
  }
  
  results[[fill_col]] <- factor(results[[fill_col]], levels = unique(results[[fill_col]]))
  results[[x_col]] <- factor(results[[x_col]], levels = unique(results[[x_col]]))
  
  results_largo <- tidyr::pivot_longer(
    results,
    cols = c("Chao1", "Shannon", "Simpson"),
    names_to = "Index",
    values_to = "value"
  )
  
  colorb_default <- if (nlevels(results[[fill_col]]) == 2) .mbm_colors_2group else .mbm_colors
  fill_scale <- if (!is.null(group_colors)) {
    ggplot2::scale_fill_manual(values = group_colors)
  } else {
    switch(palette,
           "colorb" = ggplot2::scale_fill_manual(values = colorb_default),
           "grey" = ggplot2::scale_fill_grey(start = 0.9, end = 0.3),
           "viridis" = ggplot2::scale_fill_viridis_d(option = "plasma"),
           "brewer" = ggplot2::scale_fill_brewer(palette = "Set2"),
           ggplot2::scale_fill_grey(start = 0.9, end = 0.3)
    )
  }
  
  facet_scales <- if (free_y) "free_y" else "fixed"  # Esto ahora se usa correctamente

  # When panel tags are needed (stat set) and there's no nested double
  # facet, build the grid as separate cowplot-composed subplots instead of
  # a single faceted ggplot, so the A/B/C tags land truly outside each
  # panel - the same mechanism already used by alpha_hill_corrplot and
  # beta_partition_plot - instead of trying to carve out space inside one
  # shared facet gtable.
  use_grid_compose <- !is.null(stat) && is.null(facet_by2) &&
    identical(facet_orientation, "horizontal")

  facet_config <- if (!is.null(facet_by) && !is.null(facet_by2)) {
    formula_facet <- if (facet_orientation == "horizontal") {
      stats::as.formula(paste(facet_by2, "~",  "Index +", facet_by))
    } else {
      stats::as.formula(paste("Index +", facet_by, "~", facet_by2))
    }
    ggh4x::facet_nested(
      formula_facet,
      nest_line = ggplot2::element_line(colour = "black"),
      scales = facet_scales
    )
  } else if (!is.null(facet_by)) {
    formula_facet <- if (facet_orientation == "horizontal") {
      stats::as.formula(paste(facet_by, "~ Index"))
    } else {
      stats::as.formula(paste("Index ~", facet_by))
    }
    
    if (free_y) {
      ggh4x::facet_grid2(
        formula_facet,
        scales = "free_y",
        independent = "y"
      )
    } else {
      ggh4x::facet_grid2(
        formula_facet,
        scales = "fixed"
      )
    }
  } else {
    ggplot2::facet_wrap(
      ~Index,
      ncol = if (facet_orientation == "horizontal") {
        length(unique(results_largo$Index))
      } else {
        NULL
      },
      nrow = if (facet_orientation == "vertical") {
        length(unique(results_largo$Index))
      } else {
        NULL
      },
      scales = if (free_y) "free_y" else "fixed"
    )
  }
  
  
  capa_geom <- if (type == "boxplot") {
    ggplot2::geom_boxplot(
      width = 0.5,
      position = ggplot2::position_dodge(width = 0.9)
    )
  } else if (type == "barplot") {
    ggplot2::geom_bar(
      stat = "summary",
      fun = "mean",
      position = ggplot2::position_dodge(width = 0.9),
      width = 0.5,
      color = "black"
    )
  } else {
    stop("Invalid plot type. Use 'boxplot' or 'barplot'.")
  }
  
  capa_error <- if (type == "barplot") {
    ggplot2::stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.9)
    )
  } else {
    NULL
  }
  
  aspect_ratio_theme <- if (!is.null(aspect_ratio)) {
    ggplot2::theme(aspect.ratio = aspect_ratio)
  } else if (facet_orientation == "horizontal") {
    ggplot2::theme(aspect.ratio = 1.5)
  } else {
    ggplot2::theme(aspect.ratio = 0.7)
  }

  if (use_grid_compose) {
    # Build each panel as its own small ggplot and combine with cowplot, so
    # the A/B/C tags land in cowplot's own outside-the-panel margin - the
    # same mechanism alpha_hill_corrplot/beta_partition_plot already use -
    # instead of carving space out of one shared facet gtable.
    index_levels <- unique(results_largo$Index)
    has_facet_by <- !is.null(facet_by)
    col_levels <- index_levels
    # Alphabetical order (matching ggplot2's default factor-level order,
    # i.e. what facet_grid2 used before) - not unique()'s first-appearance
    # order, which follows the row order of the input data instead.
    row_levels <- if (has_facet_by) sort(unique(as.character(results_largo[[facet_by]]))) else "__all__"
    n_rows_grid <- length(row_levels)
    n_cols_grid <- length(col_levels)
    # Shared y-axis range across all panels when scales aren't free (the
    # default), so boxplots stay visually comparable across the grid.
    y_range <- if (!free_y) range(results_largo$value, na.rm = TRUE) else NULL

    panel_letters <- if (!is.null(panel_labels)) {
      panel_labels
    } else if (identical(panel_label_case, "lower")) {
      letters[seq_len(n_rows_grid * n_cols_grid)]
    } else {
      LETTERS[seq_len(n_rows_grid * n_cols_grid)]
    }

    build_cell <- function(row_idx, col_idx) {
      index_val <- col_levels[col_idx]
      fb_val    <- row_levels[row_idx]
      cell_data <- results_largo[results_largo$Index == index_val, ]
      if (has_facet_by) cell_data <- cell_data[cell_data[[facet_by]] == fb_val, ]

      p_cell <- ggplot2::ggplot(
        cell_data,
        ggplot2::aes(x = .data[[x_col]], y = value, fill = .data[[fill_col]])
      ) +
        capa_geom + capa_error + fill_scale +
        { if (!is.null(y_range)) ggplot2::coord_cartesian(ylim = y_range) } +
        ggplot2::labs(x = NULL, y = NULL, fill = legend_title) +
        .mbm_theme(
          legend_position = legend_position,
          extra = ggplot2::theme(
            panel.grid   = ggplot2::element_blank(),
            axis.text.x  = .mbm_x_text(x_label_angle)
          )
        )
      # Only applied when the user explicitly asks for a ratio: a fixed one
      # shrinks each cowplot cell's panel to fit its slot, leaving dead space
      # around it. Left NULL each panel stretches to fill its cell.
      if (!is.null(aspect_ratio)) {
        p_cell <- p_cell + ggplot2::theme(aspect.ratio = aspect_ratio)
      }
      # Note: no per-panel aspect.ratio here (unlike the single-facet path)
      # - forcing a fixed ratio on each independently-sized cowplot cell
      # shrinks the panel to fit within its cell, leaving dead space around
      # it instead of filling the cell the way facet_grid2 used to.

      y_val <- max(cell_data$value, na.rm = TRUE) * 0.98
      x_lvls <- levels(factor(cell_data[[x_col]]))
      x_numeric <- match(cell_data[[x_col]], x_lvls)
      x_center <- mean(range(x_numeric, na.rm = TRUE))
      p_cell <- p_cell + ggpubr::stat_compare_means(
        data = cell_data, method = stat,
        mapping = ggplot2::aes(
          label = paste0("p = ", scales::label_pvalue(accuracy = 0.001)(ggplot2::after_stat(p)))
        ),
        size = 3.5, family = "serif", hide.ns = TRUE,
        label.y = y_val, label.x = x_center
      )

      show_top_strip   <- row_idx == 1
      show_right_strip <- has_facet_by && col_idx == n_cols_grid

      if (show_top_strip && show_right_strip) {
        p_cell <- p_cell +
          ggplot2::facet_grid(
            rows     = ggplot2::vars(.data[[facet_by]]),
            cols     = ggplot2::vars(Index)
          ) +
          ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = strip_color, color = "black"),
            strip.text.x = .mbm_strip_text(strip_text_bold, size = 12),
            strip.text.y = .mbm_strip_text(strip_text_bold, size = 11, angle = -90)
          )
      } else if (show_top_strip) {
        p_cell <- p_cell +
          ggplot2::facet_wrap(~Index) +
          ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = strip_color, color = "black"),
            strip.text = .mbm_strip_text(strip_text_bold, size = 12)
          )
      } else if (show_right_strip) {
        p_cell <- p_cell +
          ggplot2::facet_grid(rows = ggplot2::vars(.data[[facet_by]])) +
          ggplot2::theme(
            strip.background.y = ggplot2::element_rect(fill = strip_color, color = "black"),
            strip.text.y = .mbm_strip_text(strip_text_bold, size = 11, angle = -90)
          )
      }

      # Shared axes: y-axis text only on the first column, x-axis text only
      # on the last row, matching the previous facet_grid2 look. Made
      # invisible (colour = NA) rather than element_blank(), which collapses
      # its allotted space to zero and would make rows/columns without
      # visible axis text shorter/narrower than the ones that have it.
      if (col_idx != 1) {
        p_cell <- p_cell + ggplot2::theme(
          axis.text.y  = ggplot2::element_text(colour = NA),
          axis.ticks.y = ggplot2::element_line(colour = NA)
        )
      }
      if (row_idx != n_rows_grid) {
        p_cell <- p_cell + ggplot2::theme(
          axis.text.x  = .mbm_x_text(x_label_angle, colour = NA),
          axis.ticks.x = ggplot2::element_line(colour = NA)
        )
      }
      p_cell
    }

    grid_ncol <- n_cols_grid
    plots <- vector("list", n_rows_grid * grid_ncol)
    labels_full <- character(length(plots))
    letter_i <- 1L
    pos <- 1L
    for (row_idx in seq_len(n_rows_grid)) {
      for (col_idx in seq_len(n_cols_grid)) {
        plots[[pos]] <- build_cell(row_idx, col_idx)
        labels_full[pos] <- panel_letters[letter_i]
        letter_i <- letter_i + 1L
        pos <- pos + 1L
      }
    }

    has_legend <- show_legend
    if (has_legend) {
      leg <- cowplot::get_legend(build_cell(1, 1) + ggplot2::theme(legend.position = legend_position))
      plots <- lapply(plots, function(pl) pl + ggplot2::theme(legend.position = "none"))
    }

    panel_grid <- cowplot::plot_grid(
      plotlist         = plots,
      ncol             = grid_ncol,
      labels           = labels_full,
      label_fontfamily = "serif",
      label_fontface   = if (panel_label_bold) "bold" else "plain",
      label_size       = 14,
      label_x          = 0,
      label_y          = 1,
      hjust            = -0.2,
      vjust            = 1.3
    )

    p <- if (has_legend) {
      switch(legend_position,
        "top"    = cowplot::plot_grid(leg, panel_grid, ncol = 1, rel_heights = c(0.1, 1)),
        "left"   = cowplot::plot_grid(leg, panel_grid, nrow = 1, rel_widths = c(0.2, 1)),
        "right"  = cowplot::plot_grid(panel_grid, leg, nrow = 1, rel_widths = c(1, 0.2)),
        cowplot::plot_grid(panel_grid, leg, ncol = 1, rel_heights = c(1, 0.1))
      )
    } else {
      panel_grid
    }

    if (!is.null(title)) {
      title_grob <- cowplot::ggdraw() +
        cowplot::draw_label(title, fontface = "bold", fontfamily = "serif", size = 14)
      p <- cowplot::plot_grid(title_grob, p, ncol = 1, rel_heights = c(0.08, 1))
    }

    return(p)
  }

  # Calcular label.y ajustado por q
  p <- ggplot2::ggplot(
    results_largo,
    ggplot2::aes(x = .data[[x_col]], y = value, fill = .data[[fill_col]])
  ) +
    capa_geom +
    capa_error +
    fill_scale +
    facet_config +
    ggplot2::labs(
      title = title,
      x = x_axis_title,
      y = y_axis_title,
      fill = legend_title
    ) +
    .mbm_theme(
      legend_position = if (show_legend) legend_position else "none",
      extra = ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        panel.spacing    = grid::unit(1, "lines"),
        axis.text.x      = .mbm_x_text(x_label_angle),
        strip.text       = .mbm_strip_text(strip_text_bold),
        strip.background = ggplot2::element_rect(fill = strip_color, color = "black")
      )
    ) + aspect_ratio_theme
  
  if (!is.null(stat)) {
    split_vars <- if (!is.null(facet_by)) c("Index", facet_by) else "Index"
    p_vals_layers <- results_largo %>%
      dplyr::group_split(dplyr::across(dplyr::all_of(split_vars))) %>%
      purrr::map(~ {
        y_val <- max(.x$value, na.rm = TRUE) * 0.98
        
        # Get the unique x-axis levels and compute the central value
        x_levels <- levels(factor(.x[[x_col]]))
        x_numeric <- match(.x[[x_col]], x_levels)
        x_center <- mean(range(x_numeric, na.rm = TRUE))
        
        ggpubr::stat_compare_means(
          data = .x,
          method = stat,
          mapping = ggplot2::aes(
            label = paste0("p = ", scales::label_pvalue(accuracy = 0.001)(ggplot2::after_stat(p)))
          ),
          size = 3.5,
          family= "serif",
          hide.ns = TRUE,
          label.y = y_val,
          label.x = x_center
        )
      })
    
    p <- p + p_vals_layers
  }

  if (!is.null(stat) || !is.null(panel_labels)) {
    # Add A, B, C... style labels to the panels. Rather than reusing an
    # existing gtable row (which may not exist, e.g. the 2nd/3rd row of a
    # facet_by grid has no strip above it, only panel.spacing), a brand new,
    # dedicated, guaranteed-empty row is inserted directly above every
    # panel, so the tag never overlaps the strip text or the plotted data.
    panel_letters <- if (!is.null(panel_labels)) {
      panel_labels
    } else if (identical(panel_label_case, "lower")) letters else LETTERS
    gb <- ggplot2::ggplot_build(p)
    lay <- gb$layout$layout
    panel_names <- sprintf("panel-%d-%d", lay$ROW, lay$COL)

    g <- ggplot2::ggplotGrob(p)
    tag_height <- grid::unit(10, "mm")

    # Panels sharing the same facet row (e.g. Chao1/Shannon/Simpson side by
    # side) share the same gtable row index, so insert exactly one tag-row
    # per unique row - not one per panel, which would stack up redundant
    # blank rows. Process bottom-to-top so inserting above a lower row never
    # shifts the row index of rows still to be processed further up.
    orig_row_of_panel <- vapply(panel_names, function(nm) g$layout$t[g$layout$name == nm][1], numeric(1))
    unique_rows <- sort(unique(orig_row_of_panel), decreasing = TRUE)
    row_reps <- panel_names[match(unique_rows, orig_row_of_panel)]
    for (nm in row_reps) {
      t_now <- g$layout$t[g$layout$name == nm][1]
      g <- gtable::gtable_add_rows(g, tag_height, pos = t_now - 1)
    }

    for (i in seq_len(nrow(lay))) {
      panel_cell <- g$layout[g$layout$name == panel_names[i], ]
      if (nrow(panel_cell) == 1) {
        tag_row <- panel_cell$t - 1
        g <- gtable::gtable_add_grob(
          g,
          grid::textGrob(
            panel_letters[lay$PANEL[i]],
            x = grid::unit(2, "mm"),
            y = grid::unit(0.5, "npc"),
            hjust = 0,
            vjust = 0.5,
            gp = grid::gpar(
              fontface = if (panel_label_bold) "bold" else "plain",
              fontfamily = "serif",
              fontsize = 13
            )
          ),
          t = tag_row, l = max(1, panel_cell$l - 1), b = tag_row, r = panel_cell$r,
          z = Inf,
          name = paste0("panel-tag-", i)
        )
      }
    }
    p <- cowplot::ggdraw() + cowplot::draw_grob(g)
  }

  return(p)
}
