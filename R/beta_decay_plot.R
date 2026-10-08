#' Distance-decay of community similarity plot
#'
#' Computes the pairwise community dissimilarity (Jaccard, Morisita-Horn,
#' Bray-Curtis or another vegdist method) and the pairwise geographic distance
#' from the sample coordinates in metadata, and tests their relationship with a
#' Mantel test.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. Its first column
#'   must hold the sample IDs (the column names of `table`). It must also have
#'   the latitude and longitude columns.
#' @param lat_col Character. Name of the latitude column in \code{metadata}
#'   (decimal degrees).
#' @param lon_col Character. Name of the longitude column in \code{metadata}
#'   (decimal degrees).
#' @param group_col Character. Name of the column in \code{metadata} that
#'   defines the groups. If given, only pairs of samples from the same group are
#'   kept, with one color per group. Optional; \code{NULL} (default) for no
#'   groups.
#' @param palette Palette name (\code{"colorb"} (default), \code{"grey"},
#'   \code{"viridis"} or \code{"brewer"}) or a vector of colors, one per group.
#'   Only used when \code{group_col} is given.
#' @param distance Character. Dissimilarity metric passed to
#'   \code{vegan::vegdist}. Common options: \code{"jaccard"} (default),
#'   \code{"horn"} (Morisita-Horn / Hill q = 1 analogue), \code{"bray"}
#'   (Bray-Curtis). Any method accepted by \code{vegdist} is valid.
#'   Case-insensitive.
#' @param method Character. Correlation method for the Mantel test:
#'   \code{"spearman"} (default) or \code{"pearson"}. Case-insensitive.
#' @param permutations Integer. Number of permutations for the Mantel test.
#'   Default \code{999}.
#' @param show_lm_stats Logical. If \code{TRUE} (default), adds \eqn{R^2} to the
#'   annotation label in addition to the Mantel r, p-value, and slope.
#' @param point_color Character. Color of scatter points. Default
#'   \code{"black"}.
#' @param line_color Character. Color of the regression line. Default
#'   \code{"#D55E00"}.
#' @param point_size Numeric. Size of the points. Default \code{1}.
#' @param point_alpha Numeric (0-1). Transparency of the points. Default
#'   \code{0.5}.
#' @param annotation_size Numeric. Font size for the stats annotation.
#'   Default \code{3.5}.
#' @param x_axis_title Character. Title of the x-axis. Default \code{"Spatial
#'   distance (km)"}.
#' @param y_axis_title Character. Title of the y-axis. If \code{NULL}
#'   (default), it is built automatically.
#' @param title Character. Plot title. \code{NULL} (default) shows no title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' loc_coords <- data.frame(
#'     Loc = 1:7,
#'     lat = 19.0 + seq(0, 0.6, length.out = 7),
#'     lon = -99.0 + seq(0, 0.6, length.out = 7)
#' )
#' metadata$lat <- loc_coords$lat[match(metadata$Loc, loc_coords$Loc)]
#' metadata$lon <- loc_coords$lon[match(metadata$Loc, loc_coords$Loc)]
#'
#' # Jaccard + Spearman Mantel (default)
#' beta_decay_plot(
#'     table    = table,
#'     metadata = metadata,
#'     lat_col  = "lat",
#'     lon_col  = "lon"
#' )
#'
#' # Horn dissimilarity + Spearman Mantel
#' beta_decay_plot(
#'     table    = table,
#'     metadata = metadata,
#'     lat_col  = "lat",
#'     lon_col  = "lon",
#'     distance = "horn",
#'     method   = "spearman"
#' )
#'
#' # Separate Mantel test per location
#' beta_decay_plot(
#'     table     = table,
#'     metadata  = metadata,
#'     lat_col   = "lat",
#'     lon_col   = "lon",
#'     group_col = "Location"
#' )
beta_decay_plot <- function(
  table,
  metadata,
  lat_col,
  lon_col,
  group_col = NULL,
  palette = "colorb",
  distance = "jaccard",
  method = "spearman",
  permutations = 999,
  show_lm_stats = TRUE,
  point_color = "black",
  line_color = "#D55E00",
  point_size = 1,
  point_alpha = 0.5,
  annotation_size = 3.5,
  x_axis_title = "Spatial distance (km)",
  y_axis_title = NULL,
  title = NULL
) {
    method <- match.arg(tolower(method), c("spearman", "pearson"))
    distance <- tolower(distance)

    if (!lat_col %in% colnames(metadata)) {
        stop("`lat_col` '", lat_col, "' not found in metadata.")
    }
    if (!lon_col %in% colnames(metadata)) {
        stop("`lon_col` '", lon_col, "' not found in metadata.")
    }
    if (!is.null(group_col) && !group_col %in% colnames(metadata)) {
        stop("`group_col` '", group_col, "' not found in metadata.")
    }

    requireNamespace("vegan", quietly = TRUE)
    requireNamespace("geosphere", quietly = TRUE)

    sample_col <- colnames(metadata)[1]
    tax_idx <- grep("^taxonomy$", colnames(table), ignore.case = TRUE)
    if (length(tax_idx) != 1) {
        stop("Table must contain exactly one column named 'taxonomy'.")
    }

    common <- intersect(
        colnames(table)[-tax_idx],
        as.character(metadata[[sample_col]])
    )
    if (length(common) == 0) {
        stop("No matching sample names between table and metadata.")
    }

    otu_t <- t(as.matrix(table[, common, drop = FALSE]))
    mode(otu_t) <- "numeric"

    meta_sub <- metadata[as.character(metadata[[sample_col]]) %in% common, ,
        drop = FALSE
    ]
    meta_sub <- meta_sub[
        match(
            rownames(otu_t),
            as.character(meta_sub[[sample_col]])
        ), ,
        drop = FALSE
    ]

    groups <- if (!is.null(group_col)) as.character(meta_sub[[group_col]]) else NULL

    dis_mat <- vegan::vegdist(otu_t, method = distance)

    coords <- cbind(
        as.numeric(meta_sub[[lon_col]]),
        as.numeric(meta_sub[[lat_col]])
    )
    geo_mat <- geosphere::distm(coords) / 1000
    rownames(geo_mat) <- rownames(otu_t)
    colnames(geo_mat) <- rownames(otu_t)

    snames <- rownames(otu_t)
    dis_full <- as.matrix(dis_mat)

    idx <- which(lower.tri(dis_full), arr.ind = TRUE)
    pairs <- data.frame(
        s1 = snames[idx[, 1]],
        s2 = snames[idx[, 2]],
        dissim = dis_full[idx],
        geo_km = geo_mat[idx],
        stringsAsFactors = FALSE
    )
    pairs$similarity <- 1 - pairs$dissim

    if (!is.null(group_col)) {
        g1 <- groups[idx[, 1]]
        g2 <- groups[idx[, 2]]
        pairs$group <- ifelse(g1 == g2, g1, NA_character_)
        pairs <- pairs[!is.na(pairs$group), ]
    }

    pairs <- pairs[stats::complete.cases(pairs$similarity, pairs$geo_km), ]

    method_sym <- paste0(tools::toTitleCase(method), " r")

    .fmt_label <- function(mantel_res, slope, r2) {
        if (show_lm_stats) {
            sprintf(
                "Mantel %s = %.3f\n%s\nslope = %.4f\nR\u00b2 = %.3f",
                method_sym, mantel_res$statistic, .mbm_p_label(mantel_res$signif), slope, r2
            )
        } else {
            sprintf(
                "Mantel %s = %.3f\n%s\nslope = %.4f",
                method_sym, mantel_res$statistic, .mbm_p_label(mantel_res$signif), slope
            )
        }
    }

    .group_stats <- function(samp_idx, sub_pairs) {
        dis_sub <- stats::as.dist(dis_full[samp_idx, samp_idx, drop = FALSE])
        geo_sub <- stats::as.dist(geo_mat[samp_idx, samp_idx, drop = FALSE])
        mantel_res <- vegan::mantel(dis_sub, geo_sub,
            method = method,
            permutations = permutations
        )
        lm_fit <- stats::lm(similarity ~ geo_km, data = sub_pairs)
        list(
            label = .fmt_label(
                mantel_res, unname(stats::coef(lm_fit)[2]),
                summary(lm_fit)$r.squared
            )
        )
    }

    if (is.null(group_col)) {
        res <- .group_stats(seq_along(snames), pairs)
        ann_df <- data.frame(
            geo_km = Inf, similarity = Inf, label = res$label,
            hjust = 1.05, vjust = 1.1,
            stringsAsFactors = FALSE
        )
    } else {
        group_levels <- levels(factor(groups))
        n_groups <- length(group_levels)
        line_gap <- if (show_lm_stats) 5.5 else 4.2

        ann_rows <- lapply(seq_along(group_levels), function(i) {
            g <- group_levels[i]
            samp_idx <- which(groups == g)
            sub_pairs <- pairs[pairs$group == g, ]
            if (length(samp_idx) < 3 || nrow(sub_pairs) < 3) {
                return(NULL)
            }

            res <- .group_stats(samp_idx, sub_pairs)
            side <- if (n_groups == 1) "right" else if (i %% 2 == 1) "left" else "right"
            data.frame(
                group = g,
                geo_km = if (side == "left") -Inf else Inf,
                similarity = Inf,
                label = res$label,
                hjust = if (side == "left") -0.05 else 1.05,
                side = side,
                stringsAsFactors = FALSE
            )
        })
        ann_df <- do.call(rbind, ann_rows)
        if (is.null(ann_df) || nrow(ann_df) == 0) {
            stop(
                "No group in `group_col` has at least 3 samples with valid ",
                "coordinates; cannot run a per-group Mantel test."
            )
        }
        ann_df <- ann_df %>%
            dplyr::group_by(side) %>%
            dplyr::mutate(
                slot = dplyr::row_number(),
                vjust = 1.1 + (slot - 1) * line_gap
            ) %>%
            dplyr::ungroup() %>%
            as.data.frame()
    }

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
    y_lab <- if (!is.null(y_axis_title)) {
        y_axis_title
    } else {
        paste0(dist_name, " similarity (1 - dissimilarity)")
    }

    aes_pts <- if (is.null(group_col)) {
        ggplot2::aes(x = geo_km, y = similarity)
    } else {
        ggplot2::aes(x = geo_km, y = similarity, color = group)
    }
    aes_ann <- if (is.null(group_col)) {
        ggplot2::aes(
            x = geo_km, y = similarity, label = label,
            hjust = hjust, vjust = vjust
        )
    } else {
        ggplot2::aes(
            x = geo_km, y = similarity, label = label,
            hjust = hjust, vjust = vjust, color = group
        )
    }

    point_layer <- if (is.null(group_col)) {
        ggplot2::geom_point(
            shape = 16, size = point_size, alpha = point_alpha,
            color = point_color
        )
    } else {
        ggplot2::geom_point(shape = 16, size = point_size, alpha = point_alpha)
    }

    smooth_layer <- if (is.null(group_col)) {
        ggplot2::geom_smooth(
            method = "lm", formula = y ~ x, se = TRUE,
            color = line_color, fill = "#56B4E9", linewidth = 0.9
        )
    } else {
        ggplot2::geom_smooth(
            method = "lm", formula = y ~ x, se = TRUE,
            linewidth = 0.9
        )
    }

    color_scale <- NULL
    if (!is.null(group_col)) {
        is_named_palette <- is.character(palette) && length(palette) == 1
        colorb_default <- if (length(unique(groups)) == 2) .mbm_colors_2group else .mbm_colors
        color_scale <- if (!is_named_palette) {
            ggplot2::scale_color_manual(name = group_col, values = palette)
        } else {
            switch(palette,
                "grey" = ggplot2::scale_color_grey(name = group_col, start = 0.7, end = 0.2),
                "viridis" = ggplot2::scale_color_viridis_d(name = group_col),
                "brewer" = ggplot2::scale_color_brewer(name = group_col, palette = "Set2"),
                ggplot2::scale_color_manual(name = group_col, values = colorb_default)
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
            title = title,
            x     = x_axis_title,
            y     = y_lab
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_text(
                size = 14, color = "black", family = "serif"
            ),
            axis.title.y = ggplot2::element_text(
                size = 14, color = "black", family = "serif"
            ),
            axis.text.x = ggplot2::element_text(
                size = 12, colour = "black", family = "serif"
            ),
            axis.text.y = ggplot2::element_text(
                size = 12, color = "black", family = "serif"
            ),
            plot.title = ggplot2::element_text(
                hjust = 0.5, size = 14, color = "black", family = "serif", face = "bold"
            ),
            legend.title = ggplot2::element_text(
                size = 14, color = "black", family = "serif", face = "bold"
            ),
            legend.text = ggplot2::element_text(
                size = 12, color = "black", family = "serif"
            )
        )

    p
}
