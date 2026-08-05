#' Distance-decay of community similarity
#'
#' Computes pairwise community dissimilarity (Jaccard, Horn/Morisita-Horn,
#' Bray-Curtis, or other \code{vegan::vegdist} methods) and pairwise
#' geographic distances (Haversine formula, km) from sample coordinates
#' stored in \code{metadata}. It runs a Mantel test to evaluate the
#' relationship between community similarity (1 − dissimilarity) and
#' geographic distance, fits a linear regression, and returns a scatter plot
#' annotated with the Mantel statistic, p-value, and regression slope.
#'
#' @param table A data frame with taxa as rows and samples as columns.
#'   Must contain a column named \code{taxonomy} (any position).
#' @param metadata A data frame whose \strong{first column} contains sample
#'   identifiers matching the column names of \code{table}. Must also contain
#'   latitude and longitude columns (see \code{lat_col} and \code{lon_col}).
#' @param lat_col Character. Name of the latitude column in \code{metadata}
#'   (decimal degrees).
#' @param lon_col Character. Name of the longitude column in \code{metadata}
#'   (decimal degrees).
#' @param group_col Character or \code{NULL}. Optional categorical column in
#'   \code{metadata} (e.g. \code{"estado2"}). When supplied, sample pairs from
#'   different groups are dropped, and a separate Mantel test, regression
#'   line, and annotation are computed \strong{within each group} (matching
#'   what you'd get running \code{beta_decay_plot} once per group), all drawn
#'   on the same plot colored by group. When \code{NULL} (default), a single
#'   global Mantel test is run on all samples, as before.
#' @param palette Only used when \code{group_col} is supplied. Either a
#'   palette name (\code{"colorb"} default, \code{"grey"}, \code{"viridis"},
#'   \code{"brewer"}) or a vector of fixed colors, one per group level.
#' @param distance Character. Dissimilarity metric passed to
#'   \code{vegan::vegdist}. Common options: \code{"jaccard"} (default),
#'   \code{"horn"} (Morisita-Horn / Hill q = 1 analogue), \code{"bray"}
#'   (Bray-Curtis). Any method accepted by \code{vegdist} is valid.
#' @param method Character. Correlation method for the Mantel test.
#'   \code{"spearman"} (default) or \code{"pearson"}.
#' @param permutations Integer. Number of permutations for the Mantel test.
#'   Default \code{999}.
#' @param show_lm_stats Logical. If \code{TRUE} (default), adds R² to the
#'   annotation label in addition to the Mantel r, p-value, and slope.
#' @param point_color Character. Color of scatter points. Default
#'   \code{"black"}, matching \code{alpha_hill_corrplot}/\code{alpha_decay_plot}.
#' @param line_color Character. Color of the regression line. Default
#'   \code{"#D55E00"}, matching \code{alpha_hill_corrplot}/\code{alpha_decay_plot}'s
#'   ungrouped color scheme. The confidence-interval ribbon uses that same
#'   scheme's fill, \code{"#56B4E9"}.
#' @param point_size Numeric. Size of scatter points. Default \code{1}.
#' @param point_alpha Numeric (0–1). Transparency of scatter points.
#'   Default \code{0.5}.
#' @param annotation_size Numeric. Font size for the stats annotation.
#'   Default \code{3.5}.
#' @param x_title Character. X-axis label.
#'   Default \code{"Spatial distance (km)"}.
#' @param y_title Character or \code{NULL}. Y-axis label. If \code{NULL}
#'   (default), it is built automatically from the \code{distance} method,
#'   e.g. \code{"Jaccard similarity (1 − dissimilarity)"}.
#' @param figure_title Character or \code{NULL}. Plot title. Default
#'   \code{NULL} (no title).
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' # Jaccard + Spearman Mantel (default)
#' beta_decay_plot(
#'   table    = table_taxa2r,
#'   metadata = metas2r,
#'   lat_col  = "lat",
#'   lon_col  = "lon"
#' )
#'
#' # Horn dissimilarity + Spearman Mantel
#' beta_decay_plot(
#'   table    = table_taxa2r,
#'   metadata = metas2r,
#'   lat_col  = "lat",
#'   lon_col  = "lon",
#'   distance = "horn",
#'   method   = "spearman"
#' )
#'
#' # Separate Mantel test per geographical area
#' beta_decay_plot(
#'   table     = table_taxa2r,
#'   metadata  = metas2r,
#'   lat_col   = "lat",
#'   lon_col   = "lon",
#'   group_col = "estado2"
#' )
beta_decay_plot <- function(
    table,
    metadata,
    lat_col,
    lon_col,
    group_col       = NULL,
    palette         = "colorb",
    distance        = "jaccard",
    method          = "spearman",
    permutations    = 999,
    show_lm_stats   = TRUE,
    point_color     = "black",
    line_color      = "#D55E00",
    point_size      = 1,
    point_alpha     = 0.5,
    annotation_size = 3.5,
    x_title         = "Spatial distance (km)",
    y_title         = NULL,
    figure_title    = NULL
) {

  # ---- 0. Validate inputs ----
  method <- match.arg(method, c("spearman", "pearson"))

  if (!lat_col %in% colnames(metadata))
    stop("`lat_col` '", lat_col, "' not found in metadata.")
  if (!lon_col %in% colnames(metadata))
    stop("`lon_col` '", lon_col, "' not found in metadata.")
  if (!is.null(group_col) && !group_col %in% colnames(metadata))
    stop("`group_col` '", group_col, "' not found in metadata.")

  requireNamespace("vegan",     quietly = TRUE)
  requireNamespace("geosphere", quietly = TRUE)

  # ---- 1. Extract and align abundance table ----
  sample_col <- colnames(metadata)[1]
  tax_idx    <- grep("^taxonomy$", colnames(table), ignore.case = TRUE)
  if (length(tax_idx) != 1)
    stop("Table must contain exactly one column named 'taxonomy'.")

  common <- intersect(
    colnames(table)[-tax_idx],
    as.character(metadata[[sample_col]])
  )
  if (length(common) == 0)
    stop("No matching sample names between table and metadata.")

  # samples × taxa  (vegan expects rows = samples)
  otu_t        <- t(as.matrix(table[, common, drop = FALSE]))
  mode(otu_t)  <- "numeric"

  # align metadata rows to the same order as otu_t
  meta_sub <- metadata[as.character(metadata[[sample_col]]) %in% common, ,
                       drop = FALSE]
  meta_sub <- meta_sub[match(rownames(otu_t),
                             as.character(meta_sub[[sample_col]])), ,
                       drop = FALSE]

  # groups, aligned to the same sample order as otu_t/dis_mat/geo_mat
  groups <- if (!is.null(group_col)) as.character(meta_sub[[group_col]]) else NULL

  # ---- 2. Pairwise dissimilarity matrix ----
  dis_mat <- vegan::vegdist(otu_t, method = distance)

  # ---- 3. Geographic distance matrix (Haversine, km) ----
  coords  <- cbind(
    as.numeric(meta_sub[[lon_col]]),
    as.numeric(meta_sub[[lat_col]])
  )
  geo_mat <- geosphere::distm(coords) / 1000   # metres → km
  rownames(geo_mat) <- rownames(otu_t)
  colnames(geo_mat) <- rownames(otu_t)

  # ---- 4. Build long-format pairwise data frame ----
  snames       <- rownames(otu_t)
  dis_full     <- as.matrix(dis_mat)

  idx          <- which(lower.tri(dis_full), arr.ind = TRUE)
  pairs        <- data.frame(
    s1         = snames[idx[, 1]],
    s2         = snames[idx[, 2]],
    dissim     = dis_full[idx],
    geo_km     = geo_mat[idx],
    stringsAsFactors = FALSE
  )
  pairs$similarity <- 1 - pairs$dissim

  # With group_col, keep only within-group pairs (cross-group pairs aren't
  # part of any single group's Mantel test / regression).
  if (!is.null(group_col)) {
    g1 <- groups[idx[, 1]]
    g2 <- groups[idx[, 2]]
    pairs$group <- ifelse(g1 == g2, g1, NA_character_)
    pairs <- pairs[!is.na(pairs$group), ]
  }

  # remove pairs with NA coordinates
  pairs <- pairs[stats::complete.cases(pairs$similarity, pairs$geo_km), ]

  # ---- 5. Mantel test + regression (global, or one per group) ----
  method_sym <- paste0(tools::toTitleCase(method), " r")

  .fmt_label <- function(mantel_res, slope, r2) {
    if (show_lm_stats) {
      sprintf(
        "Mantel %s = %.3f\np = %.4f\nslope = %.4f\nR² = %.3f",
        method_sym, mantel_res$statistic, mantel_res$signif, slope, r2
      )
    } else {
      sprintf(
        "Mantel %s = %.3f\np = %.4f\nslope = %.4f",
        method_sym, mantel_res$statistic, mantel_res$signif, slope
      )
    }
  }

  # samp_idx: indices (into dis_full/geo_mat) of the samples belonging to
  # this group (or all samples, when ungrouped).
  .group_stats <- function(samp_idx, sub_pairs) {
    dis_sub <- stats::as.dist(dis_full[samp_idx, samp_idx, drop = FALSE])
    geo_sub <- stats::as.dist(geo_mat[samp_idx, samp_idx, drop = FALSE])
    mantel_res <- vegan::mantel(dis_sub, geo_sub, method = method,
                                permutations = permutations)
    lm_fit <- stats::lm(similarity ~ geo_km, data = sub_pairs)
    list(
      label = .fmt_label(mantel_res, unname(stats::coef(lm_fit)[2]),
                         summary(lm_fit)$r.squared)
    )
  }

  if (is.null(group_col)) {
    res    <- .group_stats(seq_along(snames), pairs)
    ann_df <- data.frame(
      geo_km = Inf, similarity = Inf, label = res$label,
      hjust = 1.05, vjust = 1.1,
      stringsAsFactors = FALSE
    )
  } else {
    group_levels <- levels(factor(groups))
    n_groups     <- length(group_levels)
    line_gap     <- if (show_lm_stats) 5.5 else 4.2

    ann_rows <- lapply(seq_along(group_levels), function(i) {
      g         <- group_levels[i]
      samp_idx  <- which(groups == g)
      sub_pairs <- pairs[pairs$group == g, ]
      if (length(samp_idx) < 3 || nrow(sub_pairs) < 3) return(NULL)

      res  <- .group_stats(samp_idx, sub_pairs)
      side <- if (n_groups == 1) "right" else if (i %% 2 == 1) "left" else "right"
      data.frame(
        group  = g,
        geo_km = if (side == "left") -Inf else Inf,
        similarity = Inf,
        label  = res$label,
        hjust  = if (side == "left") -0.05 else 1.05,
        side   = side,
        stringsAsFactors = FALSE
      )
    })
    ann_df <- do.call(rbind, ann_rows)
    if (is.null(ann_df) || nrow(ann_df) == 0) {
      stop("No group in `group_col` has at least 3 samples with valid ",
           "coordinates; cannot run a per-group Mantel test.")
    }
    ann_df <- ann_df %>%
      dplyr::group_by(side) %>%
      dplyr::mutate(slot = dplyr::row_number(),
                    vjust = 1.1 + (slot - 1) * line_gap) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }

  # ---- 6. Y-axis label ----
  dist_labels <- c(
    jaccard    = "Jaccard",
    horn       = "Horn",
    bray       = "Bray-Curtis",
    morisita   = "Morisita",
    kulczynski = "Kulczynski",
    raup       = "Raup-Crick",
    binomial   = "Binomial",
    cao        = "Cao",
    chao       = "Chao"
  )
  dist_name <- if (distance %in% names(dist_labels)) dist_labels[distance] else distance
  y_lab     <- if (!is.null(y_title)) y_title else
    paste0(dist_name, " similarity (1 − dissimilarity)")

  # ---- 7. Build ggplot ----
  aes_pts <- if (is.null(group_col)) {
    ggplot2::aes(x = geo_km, y = similarity)
  } else {
    ggplot2::aes(x = geo_km, y = similarity, color = group)
  }
  aes_ann <- if (is.null(group_col)) {
    ggplot2::aes(x = geo_km, y = similarity, label = label,
                 hjust = hjust, vjust = vjust)
  } else {
    ggplot2::aes(x = geo_km, y = similarity, label = label,
                 hjust = hjust, vjust = vjust, color = group)
  }

  point_layer <- if (is.null(group_col)) {
    ggplot2::geom_point(shape = 16, size = point_size, alpha = point_alpha,
                        color = point_color)
  } else {
    ggplot2::geom_point(shape = 16, size = point_size, alpha = point_alpha)
  }

  smooth_layer <- if (is.null(group_col)) {
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         color = line_color, fill = "#56B4E9", linewidth = 0.9)
  } else {
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         linewidth = 0.9)
  }

  color_scale <- NULL
  if (!is.null(group_col)) {
    is_named_palette <- is.character(palette) && length(palette) == 1
    color_scale <- if (!is_named_palette) {
      ggplot2::scale_color_manual(name = group_col, values = palette)
    } else {
      switch(palette,
        "grey"    = ggplot2::scale_color_grey(name = group_col, start = 0.7, end = 0.2),
        "viridis" = ggplot2::scale_color_viridis_d(name = group_col),
        "brewer"  = ggplot2::scale_color_brewer(name = group_col, palette = "Set2"),
        ggplot2::scale_color_manual(name = group_col, values = .mbm_colors)
      )
    }
  }

  p <- ggplot2::ggplot(pairs, aes_pts) +
    point_layer +
    smooth_layer +
    ggplot2::geom_text(
      data        = ann_df,
      mapping     = aes_ann,
      size        = annotation_size,
      family      = "serif",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    color_scale +
    ggplot2::labs(
      title = figure_title,
      x     = x_title,
      y     = y_lab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
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
      plot.title       = ggplot2::element_text(
        size = 16, color = "black", family = "serif", face = "bold"
      ),
      legend.title     = ggplot2::element_text(
        size = 14, color = "black", family = "serif", face = "bold"
      ),
      legend.text      = ggplot2::element_text(
        size = 12, color = "black", family = "serif"
      )
    )

  p
}
