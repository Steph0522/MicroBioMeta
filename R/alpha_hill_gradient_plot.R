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
#' @param fill_palette Character. Built-in palette name: \code{"colorb"}
#'   (default), \code{"grey"}, \code{"viridis"}, or \code{"brewer"}.
#' @param custom_palette A named or unnamed character vector of colors.
#'   Overrides \code{fill_palette} when provided.
#' @param x_title Character. X-axis label. Defaults to the value of
#'   \code{cont_var}.
#' @param y_title Character. Y-axis label.
#'   Default: \code{"Effective number of features"}.
#' @param figure_title Character or \code{NULL}. Overall plot title.
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
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' # All samples, no grouping
#' alpha_hill_gradient_plot(
#'   table    = table_taxa2,
#'   metadata = metas2,
#'   cont_var = "dist_km"
#' )
#'
#' # Separate regression lines by state
#' alpha_hill_gradient_plot(
#'   table     = table_taxa2,
#'   metadata  = metas2,
#'   cont_var  = "dist_km",
#'   group_col = "estado2"
#' )
alpha_hill_gradient_plot <- function(
    table,
    metadata,
    cont_var,
    group_col         = NULL,
    method            = "spearman",
    show_lm_stats     = TRUE,
    facet_orientation = "horizontal",
    fill_palette      = "colorb",
    custom_palette    = NULL,
    x_title           = NULL,
    y_title           = "Effective number of features",
    figure_title      = NULL,
    show_legend       = TRUE,
    legend_position   = "bottom",
    free_y            = TRUE,
    point_size        = 2,
    line_width        = 0.9,
    point_alpha       = 0.8,
    annotation_size   = 3.5
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
  rho_sym    <- paste0(tools::toTitleCase(method), " r")

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
      sprintf("%s = %.3f\np = %.4f\nR² = %.3f\nslope = %.3f",
              rho_sym, rho, pval, r2, slope)
    } else {
      sprintf("%s = %.3f\np = %.4f", rho_sym, rho, pval)
    }

    data.frame(rho = rho, pval = pval, r2 = r2, slope = slope,
               label = label, x_pos = -Inf, y_pos = Inf,
               stringsAsFactors = FALSE)
  }

  stats_df <- hills_long %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::group_modify(~ .compute_stats(.x)) %>%
    dplyr::ungroup()

  # ---- 6. Color scale ----
  default_colors <- c(
    "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032",
    "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
    "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26"
  )

  color_scale <- if (!is.null(custom_palette)) {
    ggplot2::scale_color_manual(values = custom_palette)
  } else {
    switch(fill_palette,
      "colorb"  = ggplot2::scale_color_manual(values = default_colors),
      "grey"    = ggplot2::scale_color_grey(start = 0.7, end = 0.2),
      "viridis" = ggplot2::scale_color_viridis_d(option = "plasma"),
      "brewer"  = ggplot2::scale_color_brewer(palette = "Set2"),
      ggplot2::scale_color_manual(values = default_colors)
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
                 color = .data[[group_col]])
  } else {
    ggplot2::aes(x = x_pos, y = y_pos, label = label)
  }

  facet_ncol <- if (facet_orientation == "horizontal") 3L else 1L

  p <- ggplot2::ggplot(hills_long, aes_pts) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::geom_smooth(
      method    = "lm",
      formula   = y ~ x,
      se        = TRUE,
      linewidth = line_width
    ) +
    ggplot2::geom_text(
      data        = stats_df,
      mapping     = aes_ann,
      hjust       = -0.05,
      vjust       = 1.1,
      size        = annotation_size,
      family      = "serif",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ~q,
      ncol   = facet_ncol,
      scales = if (free_y) "free_y" else "fixed"
    ) +
    color_scale +
    ggplot2::labs(
      title = figure_title,
      x     = if (!is.null(x_title)) x_title else cont_var,
      y     = y_title,
      color = group_col
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.spacing    = grid::unit(1, "lines"),
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
      legend.position  = if (show_legend) legend_position else "none"
    )

  p
}
