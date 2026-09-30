#' Alpha diversity plot
#'
#' This function generates a boxplot or barplot to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows.
#'              The first column must contain the OTUID, ASV, or species name.
#' @param metadata A data frame with metadata. The first column must match sample names in `table`.
#' @param type Type of plot: `"boxplot"` (default) or `"barplot"`. Case-insensitive.
#' @param stat Optional. Statistical test to compare the groups within each
#'   panel. \code{"wilcox.test"} or \code{"t.test"} compare every pair of
#'   groups, each with its own bracket and p-value (\code{ggpubr::stat_pwc()});
#'   \code{"kruskal.test"} or \code{"anova"} give one global p-value per panel.
#' @param p_adjust_method Multiple-comparison correction for the pairwise
#'   tests (\code{stat = "wilcox.test"} or \code{"t.test"}), applied within
#'   each panel; any method of \code{stats::p.adjust()}. Default
#'   \code{"holm"}; \code{"none"} shows the raw p-values.
#'   Panel tags (A, B, C...) are added regardless of whether \code{stat} is set;
#'   \code{stat} only adds the p-value annotations on top of them.
#' @param x_col Column in `metadata` to be used on the x-axis.
#' @param fill_col Column in `metadata` to define fill color.
#' @param facet_by Optional. A metadata column to facet (e.g., Treatment, Site).
#' @param facet_by2 Optional. A metadata column to double facet (e.g., Treatment, Site).
#' @param facet_orientation Whether `facet_by` appears in columns
#'   (`"horizontal"`, default) or rows (`"vertical"`). Case-insensitive.
#' @param palette Color palette to use: `"colorb"` (default, colorblind-friendly
#'   Okabe-Ito), `"grey"`, `"viridis"`, or `"brewer"`. Case-insensitive.
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
#' @param panel_label_case Character. Case of the auto-generated panel tags
#'   (A, B, C... added to every panel by default). One of
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
#' @param table_filename Character. File path/name for the saved table. Default \code{"hill.txt"}.
#'
#' @return A ggplot object showing alpha diversity with Hill numbers.
#' @export
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' # Panel tags (A/B/C) are always added (see panel_label_case and
#' # panel_labels to customize their case/format); stat additionally
#' # overlays p-value annotations on each panel

#' alpha_hill_plot(
#'   table            = table,
#'   metadata         = metadata,
#'   type             = "boxplot",
#'   x_col            = "Location",
#'   fill_col         = "Location",
#'   facet_by         = "Treatment",
#'   legend_position  = "top",
#'   stat             = "kruskal.test",
#'   panel_label_case = "upper"
#' )
alpha_hill_plot <- function(
    table,
    metadata,
    type = "boxplot",
    stat = NULL,
    p_adjust_method = "holm",
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
    y_axis_title = "Effective number of features",
    free_y = FALSE,
    panel_label_case = "upper",
    panel_labels = NULL,
    panel_label_bold = TRUE,
    x_label_angle = 0,
    strip_text_bold = FALSE,
    aspect_ratio = NULL,
    save_table = FALSE,
    table_filename = "hill.txt") {

  # Accept the fixed-choice arguments case-insensitively.
  type              <- tolower(type)
  facet_orientation <- tolower(facet_orientation)
  palette           <- tolower(palette)
  panel_label_case  <- tolower(panel_label_case)

  sample_order <- metadata[[1]]
  common_samples <- intersect(colnames(table), sample_order)
  if (length(common_samples) == 0) stop("No matching sample names between table and metadata.")

  # Keep and align only the samples present in both table and metadata,
  # in the same order, so the positional cbind below (results <-
  # data.frame(q0, q1, q2, metadata)) lines up correctly.
  table <- table[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata[[1]]), , drop = FALSE]
  table <- data.frame(t(table))
  
  results <- data.frame(
    q0 = hillR::hill_taxa(comm = table, q = 0),
    q1 = hillR::hill_taxa(comm = table, q = 1),
    q2 = hillR::hill_taxa(comm = table, q = 2),
    metadata
  )
  
  # Guardar tabla si se solicita
  if (save_table) {
    utils::write.table(
      results,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message("Table saved as: ", table_filename)
  }
  
  
  #results[[fill_col]] <- factor(results[[fill_col]], levels = unique(results[[fill_col]]))
  #results[[x_col]] <- factor(results[[x_col]], levels = unique(results[[x_col]]))
  
  results_largo <- tidyr::pivot_longer(
    results,
    cols = tidyselect::starts_with("q"),
    names_to = "q",
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
  q_labeller <- .mbm_q_labeller(strip_text_bold)

  # Panel tags (A/B/C) are always added. When there's no nested double facet,
  # build the grid as separate cowplot-composed subplots instead of a single
  # faceted ggplot, so the tags land truly outside each panel - the same
  # mechanism already used by alpha_hill_corr_plot and beta_partition_ord_plot -
  # instead of trying to carve out space inside one shared facet gtable.
  # Also used (with q stacked in rows instead of columns) for the vertical
  # orientation when there's no facet_by - the plain facet_wrap() version of
  # that case could never get closer than a small panel.spacing gap between
  # q's stacked strip+panel units, unlike this grid's seamless cowplot look.
  use_grid_compose <- is.null(facet_by2) &&
    (identical(facet_orientation, "horizontal") || is.null(facet_by))

  facet_config <- if (!is.null(facet_by) && !is.null(facet_by2)) {
    formula_facet <- if (facet_orientation == "horizontal") {
      stats::as.formula(paste(facet_by2, "~",  "q +", facet_by))
    } else {
      stats::as.formula(paste("q +", facet_by, "~", facet_by2))
    }
    ggh4x::facet_nested(
      formula_facet,
      nest_line = ggplot2::element_line(colour = "black"),
      scales = facet_scales,
      labeller = ggplot2::labeller(q = q_labeller),
      axes = "x"
    )
  } else if (!is.null(facet_by)) {
    formula_facet <- if (facet_orientation == "horizontal") {
      stats::as.formula(paste(facet_by, "~ q"))
    } else {
      stats::as.formula(paste("q ~", facet_by))
    }

    if (free_y) {
      ggh4x::facet_grid2(
        formula_facet,
        scales = "free_y",
        independent = "y",
        labeller = ggplot2::labeller(q = q_labeller),
        axes = "x"
      )
    } else {
      ggh4x::facet_grid2(
        formula_facet,
        scales = "fixed",
        labeller = ggplot2::labeller(q = q_labeller),
        axes = "x"
      )
    }
  } else {
    ggplot2::facet_wrap(
      ~q,
      ncol = if (facet_orientation == "horizontal") {
        length(unique(results_largo$q))
      } else {
        NULL
      },
      nrow = if (facet_orientation == "vertical") {
        length(unique(results_largo$q))
      } else {
        NULL
      },
      scales = if (free_y) "free_y" else "fixed",
      labeller = q_labeller
      # Default axes ("margins"): x-axis text only on the bottom-most panel
      # of each column - already true for the horizontal layout (every panel
      # sits at the bottom of its own column) and, for vertical, matches the
      # use_grid_compose look (facet_by rows only labelled on the last row)
      # instead of repeating the x categories on every stacked panel, which
      # was inflating the visual gap between q0/q1/q2.
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
    # same mechanism alpha_hill_corr_plot/beta_partition_ord_plot already use -
    # instead of carving space out of one shared facet gtable.
    q_levels <- intersect(c("q0", "q1", "q2"), unique(results_largo$q))
    has_facet_by <- !is.null(facet_by)
    # q goes in rows only for the no-facet_by + vertical case; otherwise it
    # stays in columns (facet_by, when present, always takes the rows).
    q_in_rows <- !has_facet_by && identical(facet_orientation, "vertical")
    col_levels <- if (q_in_rows) "__all__" else q_levels
    # Alphabetical order (matching ggplot2's default factor-level order,
    # i.e. what facet_grid2 used before) - not unique()'s first-appearance
    # order, which follows the row order of the input data instead.
    row_levels <- if (has_facet_by) {
      sort(unique(as.character(results_largo[[facet_by]])))
    } else if (q_in_rows) {
      q_levels
    } else {
      "__all__"
    }
    n_rows_grid <- length(row_levels)
    n_cols_grid <- length(col_levels)
    # Shared y-axis range across all panels when scales aren't free (the
    # default), so boxplots stay visually comparable across the grid.
    y_range <- if (!free_y) range(results_largo$value, na.rm = TRUE) else NULL
    # long y titles go on two lines when there are several (short) rows
    y_title_cell <- .mbm_row_axis_title(y_axis_title, n_rows_grid)

    panel_letters <- if (!is.null(panel_labels)) {
      panel_labels
    } else if (identical(panel_label_case, "lower")) {
      letters[seq_len(n_rows_grid * n_cols_grid)]
    } else {
      LETTERS[seq_len(n_rows_grid * n_cols_grid)]
    }

    build_cell <- function(row_idx, col_idx) {
      q_val  <- if (q_in_rows) row_levels[row_idx] else col_levels[col_idx]
      fb_val <- if (has_facet_by) row_levels[row_idx] else NULL
      cell_data <- results_largo[results_largo$q == q_val, ]
      if (has_facet_by) cell_data <- cell_data[cell_data[[facet_by]] == fb_val, ]

      p_cell <- ggplot2::ggplot(
        cell_data,
        ggplot2::aes(x = .data[[x_col]], y = value, fill = .data[[fill_col]])
      ) +
        capa_geom + capa_error + fill_scale +
        { if (!is.null(y_range)) ggplot2::coord_cartesian(ylim = y_range) } +
        # q in rows (a2-style, no facet_by): every panel is one of only 3
        # total, so the title repeats on all of them. q in columns (a1-style,
        # facet_by present): only the first (leftmost) column shows it -
        # same column whose axis.text.y ticks are the only ones left visible
        # below, since every facet_by row down that column reads the same
        # units/legend as the row above it.
        ggplot2::labs(
          # x title only on the bottom row (patchwork then merges the
          # identical ones into one), as the y title is only on column 1
          x = if (row_idx == n_rows_grid) x_axis_title else NULL,
          y = if (q_in_rows || col_idx == 1) y_title_cell else NULL,
          fill = legend_title
        ) +
        .mbm_theme(
          legend_position = legend_position,
          extra = ggplot2::theme(
            panel.grid   = ggplot2::element_blank(),
            axis.text.x  = .mbm_x_text(x_label_angle),
            # Smaller than the package default (12) so the topmost y-axis
            # tick number doesn't collide with the panel tag letter, which
            # cowplot draws in the panel's top-left corner overlapping the
            # axis area rather than in a reserved margin.
            axis.text.y  = ggplot2::element_text(size = 9, color = "black")
          )
        )
      # A fixed aspect.ratio shrinks each cowplot cell's panel to fit inside
      # its slot, leaving dead space around it, so it is only applied when the
      # user explicitly asks for one. Left NULL (the default) each panel
      # stretches to fill its cell the way facet_grid2 used to.
      if (!is.null(aspect_ratio)) {
        p_cell <- p_cell + ggplot2::theme(aspect.ratio = aspect_ratio)
      }

      if (!is.null(stat)) {
        # label.y is placed above the data max, not at 0.98x it (which sits
        # right at the top whisker/outlier and collides with it); extra top
        # expansion gives it room so it isn't clipped by the panel border.
        data_range <- range(cell_data$value, na.rm = TRUE)
        y_val <- data_range[2] + diff(data_range) * 0.12
        x_lvls <- levels(factor(cell_data[[x_col]]))
        x_numeric <- match(cell_data[[x_col]], x_lvls)
        x_center <- mean(range(x_numeric, na.rm = TRUE))
        p_cell <- p_cell +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
          .mbm_stat_layer(stat, p_adjust_method, data = cell_data,
                         label.x = x_center, label.y = y_val)
      }

      # Every row needs its own q strip when q is stacked in rows (each row
      # is a different q, unlike the facet_by case where only row 1 sits
      # under the shared q header row).
      show_top_strip   <- if (q_in_rows) TRUE else row_idx == 1
      show_right_strip <- has_facet_by && col_idx == n_cols_grid

      if (show_top_strip && show_right_strip) {
        p_cell <- p_cell +
          ggplot2::facet_grid(
            rows     = ggplot2::vars(.data[[facet_by]]),
            cols     = ggplot2::vars(q),
            labeller = ggplot2::labeller(q = q_labeller)
          ) +
          ggplot2::theme(
            strip.background = ggplot2::element_rect(fill = strip_color, color = "black"),
            strip.text.x = .mbm_strip_text(strip_text_bold, size = 12),
            strip.text.y = .mbm_strip_text(strip_text_bold, size = 11, angle = -90)
          )
      } else if (show_top_strip) {
        p_cell <- p_cell +
          ggplot2::facet_wrap(~q, labeller = q_labeller) +
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

    plots <- list()
    for (row_idx in seq_len(n_rows_grid)) {
      for (col_idx in seq_len(n_cols_grid)) {
        plots[[length(plots) + 1L]] <- build_cell(row_idx, col_idx)
      }
    }

    # Joined with patchwork (see .mbm_patchwork_grid()), so the result stays
    # modifiable: `p & theme(...)`, `p[[2]] + labs(...)`, `+ plot_annotation()`.
    # Axis titles stay on each panel; to merge them into one, add
    # `+ patchwork::plot_layout(axis_titles = "collect")` to the result.
    p <- .mbm_patchwork_grid(plots, ncol = n_cols_grid, tags = panel_letters,
                             title = title, show_legend = show_legend,
                             legend_position = legend_position,
                             tag_bold = panel_label_bold)

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
        # The strip is always fused to its own panel below (ggplot never
        # inserts a gap there); this is the only gap ggplot has, and it sits
        # between a panel and the NEXT facet's strip above it. Too small
        # (e.g. 0.15) and that gap reads as "attached on both sides" instead
        # of only to the panel below.
        panel.spacing    = grid::unit(0.5, "lines"),
        axis.text.x      = .mbm_x_text(x_label_angle),
        strip.text       = .mbm_strip_text(strip_text_bold),
        strip.background = ggplot2::element_rect(fill = strip_color, color = "black")
      )
    ) + aspect_ratio_theme
  
  if (!is.null(stat)) {
    split_vars <- if (!is.null(facet_by)) c("q", facet_by) else "q"
    p_vals_layers <- results_largo %>%
      dplyr::group_split(dplyr::across(dplyr::all_of(split_vars))) %>%
      purrr::map(~ {
        # Above the data max, not at 0.98x it (which sits right at the top
        # whisker/outlier and collides with it).
        data_range <- range(.x$value, na.rm = TRUE)
        y_val <- data_range[2] + diff(data_range) * 0.12

        # Get the unique x-axis levels and compute the central value
        x_levels <- levels(factor(.x[[x_col]]))
        x_numeric <- match(.x[[x_col]], x_levels)
        x_center <- mean(range(x_numeric, na.rm = TRUE))

        .mbm_stat_layer(stat, p_adjust_method, data = .x,
                       label.x = x_center, label.y = y_val)
      })

    p <- p +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
      p_vals_layers
  }

  # A, B, C... tags on each panel. The plot stays a normal (faceted) ggplot,
  # modifiable with + theme()/labs(), instead of a gtable edited into an image.
  panel_letters <- if (!is.null(panel_labels)) {
    panel_labels
  } else if (identical(panel_label_case, "lower")) {
    letters
  } else {
    LETTERS
  }
  p <- .mbm_facet_tags(p, panel_letters, bold = panel_label_bold)

  return(p)
}
