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
#'   \code{"#566573"}.
#' @param line_color Character. Color of the regression line. Default
#'   \code{"#2166AC"}.
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
beta_decay_plot <- function(
    table,
    metadata,
    lat_col,
    lon_col,
    distance        = "jaccard",
    method          = "spearman",
    permutations    = 999,
    show_lm_stats   = TRUE,
    point_color     = "#566573",
    line_color      = "#2166AC",
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

  # ---- 4. Mantel test ----
  geo_dist_obj <- stats::as.dist(geo_mat)
  mantel_res   <- vegan::mantel(dis_mat, geo_dist_obj,
                                method      = method,
                                permutations = permutations)

  # ---- 5. Build long-format pairwise data frame ----
  snames       <- rownames(otu_t)
  n            <- length(snames)
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

  # remove pairs with NA coordinates
  pairs <- pairs[stats::complete.cases(pairs$similarity, pairs$geo_km), ]

  # ---- 6. Linear regression ----
  lm_fit <- stats::lm(similarity ~ geo_km, data = pairs)
  slope  <- unname(stats::coef(lm_fit)[2])
  r2     <- summary(lm_fit)$r.squared

  # ---- 7. Annotation label ----
  method_sym <- paste0(tools::toTitleCase(method), " r")
  ann_label  <- if (show_lm_stats) {
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

  ann_df <- data.frame(
    geo_km     = -Inf,
    similarity = Inf,
    label      = ann_label,
    stringsAsFactors = FALSE
  )

  # ---- 8. Y-axis label ----
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

  # ---- 9. Build ggplot ----
  p <- ggplot2::ggplot(pairs, ggplot2::aes(x = geo_km, y = similarity)) +
    ggplot2::geom_point(
      shape = 16, size = point_size, alpha = point_alpha, color = point_color
    ) +
    ggplot2::geom_smooth(
      method    = "lm",
      formula   = y ~ x,
      se        = TRUE,
      color     = line_color,
      linewidth = 0.9
    ) +
    ggplot2::geom_text(
      data        = ann_df,
      mapping     = ggplot2::aes(x = geo_km, y = similarity, label = label),
      hjust       = -0.05,
      vjust       = 1.1,
      size        = annotation_size,
      family      = "serif",
      inherit.aes = FALSE
    ) +
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
      )
    )

  p
}
