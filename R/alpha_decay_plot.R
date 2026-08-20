#' Alpha diversity along a continuous gradient
#'
#' Computes Hill numbers (q = 0, 1, 2) from an ASV/OTU table and plots them
#' against a continuous metadata variable (e.g. distance to urban center,
#' elevation, pH). Each Hill order is shown in its own facet panel with an
#' ordinary least-squares regression line and an annotation reporting the
#' Spearman rank correlation coefficient (rho, or Pearson r), its p-value,
#' and optionally the linear regression R² and slope.
#' An optional grouping variable adds per-group coloring and separate
#' regression lines.
#'
#' @param table A data frame with taxa as rows and samples as columns.
#'   Must contain a column named \code{taxonomy} (any position).
#' @param metadata A data frame whose \strong{first column} contains sample
#'   identifiers matching the column names of \code{table}.
#' @param cont_var Character. Name of the continuous variable column in
#'   \code{metadata} to place on the x-axis (e.g. \code{"dist_km"}).
#' @param group_col Character or \code{NULL}. Optional column in
#'   \code{metadata} used to color points and fit separate regression lines
#'   per group (e.g. \code{"estado2"}).
#' @param method Character. Correlation method for the statistic annotation.
#'   One of \code{"spearman"} (default) or \code{"pearson"}.
#' @param show_lm_stats Logical. If \code{TRUE} (default), adds R² and slope
#'   from the linear model to the annotation label.
#' @param facet_orientation Character. \code{"horizontal"} (default) places
#'   q-panels in a single row; \code{"vertical"} stacks them in one column.
#' @param palette Character. Built-in palette name: \code{"colorb"}
#'   (default), \code{"grey"}, \code{"viridis"}, or \code{"brewer"}.
#' @param group_colors A named or unnamed character vector of colors.
#'   Overrides \code{palette} when provided.
#' @param x_axis_title Character. X-axis label. Defaults to the value of
#'   \code{cont_var}.
#' @param y_axis_title Character. Y-axis label.
#'   Default: \code{"Effective number of features"}.
#' @param title Character or \code{NULL}. Overall plot title.
#' @param show_legend Logical. Show the color legend? Default \code{TRUE}.
#' @param legend_position Character. Legend position: \code{"bottom"}
#'   (default), \code{"top"}, \code{"right"}, or \code{"left"}.
#' @param free_y Logical. Use free y-axis scales across facets?
#'   Default \code{TRUE}.
#' @param point_size Numeric. Size of scatter points. Default \code{2}.
#' @param line_width Numeric. Width of regression lines. Default \code{0.9}.
#' @param point_alpha Numeric (0–1). Transparency of points. Default \code{0.8}.
#' @param annotation_size Numeric. Font size for the stats annotation.
#'   Default \code{3.5}.
#' @param panel_label_case Character. Case of the auto-generated A/B/C panel
#'   tags. One of \code{"upper"} (default, "A", "B", "C") or \code{"lower"}
#'   ("a", "b", "c"). Ignored if \code{panel_labels} is supplied.
#' @param panel_labels Optional character vector of 3 custom panel tags (one
#'   per q0/q1/q2 panel), used as-is (e.g. \code{c("(a)", "(b)", "(c)")} or
#'   \code{c("a.", "b.", "c.")}) — for journal styles that
#'   \code{panel_label_case} alone can't produce. Overrides
#'   \code{panel_label_case} when provided.
#' @param panel_label_bold Logical. If \code{TRUE} (default), panel tags are
#'   bold. Set to \code{FALSE} for journals that require plain (non-bold)
#'   panel tags.
#'
#' @return A \code{ggplot} object (a \code{cowplot} composite of the three
#'   q0/q1/q2 panels, always tagged A/B/C).
#' @export
#'
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' # alpha_decay_plot expects a continuous environmental gradient (e.g.
#' # distance, elevation, pH). This bundled example dataset doesn't include
#' # one, so this creates an illustrative synthetic "distance to a reference
#' # point" (km), spread across the 'Loc' site codes with a little
#' # per-sample jitter - substitute your own real gradient variable.
#' set.seed(1)
#' loc_dist <- data.frame(Loc = 1:7, dist_km = seq(0, 12, length.out = 7))
#' metadata$dist_km <- loc_dist$dist_km[match(metadata$Loc, loc_dist$Loc)] +
#'   stats::rnorm(nrow(metadata), sd = 0.3)
#'
#' # All samples, no grouping. Always tagged A/B/C (see panel_label_case and
#' # panel_labels to customize)
#' alpha_decay_plot(
#'   table        = table,
#'   metadata     = metadata,
#'   cont_var     = "dist_km",
#'   x_axis_title = "Distance (km)"
#' )
#'
#' # Separate regression lines by location
#' alpha_decay_plot(
#'   table        = table,
#'   metadata     = metadata,
#'   cont_var     = "dist_km",
#'   x_axis_title = "Distance (km)",
#'   group_col    = "Location"
#' )
#' }
alpha_decay_plot <- function(
    table,
    metadata,
    cont_var,
    group_col         = NULL,
    method            = "spearman",
    show_lm_stats     = TRUE,
    facet_orientation = "horizontal",
    palette           = "colorb",
    group_colors      = NULL,
    x_axis_title      = NULL,
    y_axis_title      = "Effective number of features",
    title             = NULL,
    show_legend       = TRUE,
    legend_position   = "bottom",
    free_y            = TRUE,
    point_size        = 2,
    line_width        = 0.9,
    point_alpha       = 0.8,
    annotation_size   = 3.5,
    panel_label_case  = "upper",
    panel_labels      = NULL,
    panel_label_bold  = TRUE
) {

  # ---- 0. Validate inputs ----
  method <- match.arg(method, c("spearman", "pearson"))

  if (!cont_var %in% colnames(metadata))
    stop("`cont_var` '", cont_var, "' not found in metadata.")
  if (!is.null(group_col) && !group_col %in% colnames(metadata))
    stop("`group_col` '", group_col, "' not found in metadata.")

  # ---- 1. Align table and metadata ----
  sample_col <- colnames(metadata)[1]
  tax_idx    <- grep("^taxonomy$", colnames(table), ignore.case = TRUE)
  if (length(tax_idx) != 1)
    stop("Table must contain exactly one column named 'taxonomy'.")

  common <- intersect(colnames(table)[-tax_idx],
                      as.character(metadata[[sample_col]]))
  if (length(common) == 0)
    stop("No matching sample names between table and metadata.")

  counts_t <- data.frame(t(table[, common, drop = FALSE]))

  # ---- 2. Compute Hill numbers ----
  hills <- data.frame(
    .sample = rownames(counts_t),
    q0      = hillR::hill_taxa(comm = counts_t, q = 0),
    q1      = hillR::hill_taxa(comm = counts_t, q = 1),
    q2      = hillR::hill_taxa(comm = counts_t, q = 2),
    stringsAsFactors = FALSE
  )

  # ---- 3. Merge with metadata ----
  meta_sub              <- metadata
  meta_sub[[sample_col]] <- as.character(meta_sub[[sample_col]])
  meta_sub[[cont_var]]   <- suppressWarnings(as.numeric(meta_sub[[cont_var]]))

  hills_meta <- dplyr::left_join(hills, meta_sub,
                                 by = stats::setNames(sample_col, ".sample"))

  # ---- 4. Pivot to long ----
  hills_long <- tidyr::pivot_longer(
    hills_meta,
    cols      = c("q0", "q1", "q2"),
    names_to  = "q",
    values_to = "hill"
  )
  hills_long <- hills_long[!is.na(hills_long[[cont_var]]) &
                              !is.na(hills_long$hill), ]

  # ---- 5. Compute per-group stats ----
  group_vars <- if (!is.null(group_col)) c("q", group_col) else "q"
  rho_sym    <- if (method == "spearman") "rho" else "r"

  # Helper that receives one sub-data frame and returns a 1-row stats data frame
  .compute_stats <- function(df) {
    x  <- df[[cont_var]]
    y  <- df$hill
    ok <- complete.cases(x, y)
    x  <- x[ok]; y <- y[ok]

    ct <- tryCatch(
      stats::cor.test(x, y, method = method),
      error = function(e) NULL
    )
    m <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)

    rho   <- if (is.null(ct)) NA_real_ else unname(ct$estimate)
    pval  <- if (is.null(ct)) NA_real_ else ct$p.value
    r2    <- if (is.null(m))  NA_real_ else summary(m)$r.squared
    slope <- if (is.null(m))  NA_real_ else unname(stats::coef(m)[2])

    label <- if (show_lm_stats) {
      sprintf("%s = %.3f\np = %s\nR² = %.3f\nslope = %.3f",
              rho_sym, rho, .mbm_format_pval(pval), r2, slope)
    } else {
      sprintf("%s = %.3f\np = %s", rho_sym, rho, .mbm_format_pval(pval))
    }

    data.frame(rho = rho, pval = pval, r2 = r2, slope = slope,
               label = label, stringsAsFactors = FALSE)
  }

  stats_df <- hills_long %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::group_modify(~ .compute_stats(.x)) %>%
    dplyr::ungroup()

  # ---- 5b. Position annotations: right if alone; alternate left/right
  #          (group 1 left, group 2 right, ...) when there is more than one ----
  line_gap <- if (show_lm_stats) 5.5 else 2.7

  if (!is.null(group_col)) {
    group_levels <- levels(factor(hills_long[[group_col]]))
    n_groups     <- length(group_levels)

    stats_df <- stats_df %>%
      dplyr::mutate(
        .grp_idx = match(.data[[group_col]], group_levels),
        side     = if (n_groups == 1) "right" else
          ifelse(.grp_idx %% 2 == 1, "left", "right")
      ) %>%
      dplyr::arrange(q, side, .grp_idx) %>%
      dplyr::group_by(q, side) %>%
      dplyr::mutate(slot = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        x_pos    = ifelse(side == "left", -Inf, Inf),
        hjust    = ifelse(side == "left", -0.05, 1.05),
        vjust    = 1.1 + (slot - 1) * line_gap,
        y_pos    = Inf,
        .grp_idx = NULL,
        side     = NULL,
        slot     = NULL
      )
  } else {
    stats_df <- stats_df %>%
      dplyr::mutate(x_pos = Inf, hjust = 1.05, vjust = 1.1, y_pos = Inf)
  }

  # ---- 6. Color scale ----
  colorb_default <- if (!is.null(group_col) && length(unique(hills_long[[group_col]])) == 2) {
    .mbm_colors_2group
  } else {
    .mbm_colors
  }
  color_scale <- if (!is.null(group_colors)) {
    ggplot2::scale_color_manual(values = group_colors)
  } else {
    switch(palette,
      "colorb"  = ggplot2::scale_color_manual(values = colorb_default),
      "grey"    = ggplot2::scale_color_grey(start = 0.7, end = 0.2),
      "viridis" = ggplot2::scale_color_viridis_d(option = "plasma"),
      "brewer"  = ggplot2::scale_color_brewer(palette = "Set2"),
      ggplot2::scale_color_manual(values = colorb_default)
    )
  }

  # ---- 7. Build plot ----
  aes_pts <- if (!is.null(group_col)) {
    ggplot2::aes(
      x     = .data[[cont_var]],
      y     = hill,
      color = .data[[group_col]]
    )
  } else {
    ggplot2::aes(x = .data[[cont_var]], y = hill)
  }

  # Annotation mapping (color per group if requested)
  aes_ann <- if (!is.null(group_col)) {
    ggplot2::aes(x = x_pos, y = y_pos, label = label,
                 hjust = hjust, vjust = vjust,
                 color = .data[[group_col]])
  } else {
    ggplot2::aes(x = x_pos, y = y_pos, label = label,
                 hjust = hjust, vjust = vjust)
  }

  facet_ncol <- if (facet_orientation == "horizontal") 3L else 1L
  q_labeller <- ggplot2::as_labeller(.mbm_q_labels, default = ggplot2::label_parsed)

  # When ungrouped, match alpha_hill_corrplot's fixed reg.line/CI colors;
  # when grouped, let each group keep its own palette color.
  smooth_layer <- if (is.null(group_col)) {
    ggplot2::geom_smooth(
      method    = "lm",
      formula   = y ~ x,
      se        = TRUE,
      linewidth = line_width,
      color     = "#D55E00",
      fill      = "#56B4E9"
    )
  } else {
    ggplot2::geom_smooth(
      method    = "lm",
      formula   = y ~ x,
      se        = TRUE,
      linewidth = line_width
    )
  }

  # Shared y-axis range across the three q panels when scales aren't free
  y_range <- if (!free_y) range(hills_long$hill, na.rm = TRUE) else NULL

  # ---- 7. Build one panel per q level, then combine with cowplot so the
  # A/B/C tags land outside each panel (matching alpha_hill_corrplot) ----
  q_levels <- intersect(c("q0", "q1", "q2"), unique(hills_long$q))

  build_q_panel <- function(q_level) {
    ggplot2::ggplot(hills_long[hills_long$q == q_level, ], aes_pts) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      smooth_layer +
      ggplot2::geom_text(
        data        = stats_df[stats_df$q == q_level, ],
        mapping     = aes_ann,
        size        = annotation_size,
        family      = "serif",
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      ggplot2::facet_wrap(~q, ncol = 1, labeller = q_labeller) +
      color_scale +
      { if (!is.null(y_range)) ggplot2::coord_cartesian(ylim = y_range) } +
      ggplot2::labs(
        x     = if (!is.null(x_axis_title)) x_axis_title else cont_var,
        y     = y_axis_title,
        color = group_col
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        strip.text       = ggplot2::element_text(
          size = 12, color = "black", family = "serif", face = "bold"
        ),
        strip.background = ggplot2::element_rect(fill = "grey"),
        axis.title.x     = ggplot2::element_text(
          size = 14, color = "black", family = "serif"
        ),
        axis.title.y     = ggplot2::element_text(
          size = 14, color = "black", family = "serif"
        ),
        axis.text.x      = ggplot2::element_text(
          size = 12, colour = "black", family = "serif"
        ),
        axis.text.y      = ggplot2::element_text(
          size = 12, color = "black", family = "serif"
        ),
        legend.title     = ggplot2::element_text(
          size = 14, color = "black", family = "serif", face = "bold"
        ),
        legend.text      = ggplot2::element_text(
          size = 12, color = "black", family = "serif"
        ),
        legend.position  = legend_position
      )
  }

  panel_letters <- if (!is.null(panel_labels)) {
    panel_labels
  } else if (identical(panel_label_case, "lower")) {
    letters[seq_along(q_levels)]
  } else {
    LETTERS[seq_along(q_levels)]
  }

  plots <- lapply(q_levels, build_q_panel)

  has_legend <- !is.null(group_col) && show_legend
  if (has_legend) {
    leg <- cowplot::get_legend(plots[[1]])
    plots <- lapply(plots, function(pl) pl + ggplot2::theme(legend.position = "none"))
  }

  panel_grid <- cowplot::plot_grid(
    plotlist         = plots,
    ncol              = facet_ncol,
    labels            = panel_letters,
    label_fontfamily  = "serif",
    label_fontface    = if (panel_label_bold) "bold" else "plain",
    label_size        = 14,
    label_x           = 0,
    label_y           = 1,
    hjust             = -0.2,
    vjust             = 1.3
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

  p
}
