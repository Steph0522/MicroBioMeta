#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details. Re-exported
#' here so every function's \code{@examples} (most of which chain steps with
#' \code{\%>\%}) work with just \code{library(MicroBioMeta)}, without also
#' requiring \code{library(dplyr)}.
#'
#' @param lhs A value.
#' @param rhs A function call using the magrittr semantics.
#' @return \code{rhs(lhs)}, i.e. \code{lhs} piped into \code{rhs}.
#' @importFrom dplyr %>%
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @usage lhs \%>\% rhs
#' @examples
#' c(1, 2, 3) %>% sum()
NULL

# Column/aesthetic names used inside dplyr and ggplot2 non-standard evaluation
# (data-masked expressions) that R's static checker can't see as bound. Declared
# here so `R CMD check` doesn't flag them as undefined global variables.
utils::globalVariables(c(
  ".grp_idx", "Index", "Variable", "geo_km", "group", "hill", "hjust",
  "row_id", "side", "similarity", "slot", "vjust", "x_pos", "y_pos"
))

# --- Taxonomy string parser -------------------------------------------------
# Splits a QIIME2-style Feature.ID/Taxon table into one column per rank
# (Kingdom..Species). Adapted from qiime2R::parse_taxonomy() (MIT License,
# Copyright (c) 2018 Jordan Bisanz, https://github.com/jbisanz/qiime2R) so
# ancombc_plot() doesn't need qiime2R (a GitHub-only package) just to split a
# string already sitting in the `taxonomy` column of MicroBioMeta's own
# tables - no artifact reading involved.
.mbm_parse_taxonomy <- function(taxonomy, tax_sep = "; |;", trim_extra = TRUE) {
  if (sum(colnames(taxonomy) %in% c("Feature.ID", "Taxon")) != 2) {
    stop("Table does not match expected format, i.e. does not have columns Feature.ID and Taxon.")
  }

  taxonomy <- taxonomy[, c("Feature.ID", "Taxon")]
  if (trim_extra) {
    taxonomy$Taxon <- gsub("[kdpcofgs]__", "", taxonomy$Taxon) # GreenGenes/SILVA/Kraken2-style prefixes
    taxonomy$Taxon <- gsub("D_\\d__", "", taxonomy$Taxon)      # SILVA D_0__ style prefixes
  }
  taxonomy <- suppressWarnings(
    tidyr::separate(taxonomy, Taxon,
      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = tax_sep, fill = "right", extra = "merge")
  )
  taxonomy <- apply(taxonomy, 2, function(x) ifelse(x == "", NA_character_, x))
  taxonomy <- as.data.frame(taxonomy)
  rownames(taxonomy) <- taxonomy$Feature.ID
  taxonomy$Feature.ID <- NULL
  taxonomy
}

# --- Ordination "spider" plot data ------------------------------------------
# Builds the site scores, group centroids and per-sample-to-centroid segments
# needed for a ggplot2 spider plot from a vegan::betadisper() object - the
# specific subset of what ggordiplots::gg_ordiplot(spiders = TRUE, ellipse =
# FALSE, hull = FALSE) computes that beta_partition_ord_plot() actually uses.
# Written independently against vegan's public scores()/eigenvals() API
# (rather than adapted from ggordiplots' GPL-licensed source, which isn't
# compatible with this package's Artistic-2.0 license) so the two are
# expected to differ in implementation while matching in numeric output.
.mbm_betadisper_spider_df <- function(betadisper_obj, groups) {
  groups <- factor(groups)

  n_axes <- length(vegan::eigenvals(betadisper_obj))
  if (n_axes < 2) {
    stop(
      "This dissimilarity partition only has ", n_axes, " usable ordination axis, ",
      "so a 2D spider plot can't be drawn for it. This happens when one component ",
      "of the beta-diversity partition (often nestedness) is close to zero across ",
      "the whole dataset - i.e. it's a property of this data, not a fixable bug.",
      call. = FALSE
    )
  }

  df_ord <- vegan::scores(betadisper_obj, display = "sites", choices = c(1, 2))
  df_ord <- data.frame(x = df_ord[, 1], y = df_ord[, 2], Group = groups)

  df_mean.ord <- stats::aggregate(df_ord[, c("x", "y")], by = list(Group = df_ord$Group), mean)

  df_spiders <- df_ord
  df_spiders$cntr.x <- df_mean.ord$x[match(df_spiders$Group, df_mean.ord$Group)]
  df_spiders$cntr.y <- df_mean.ord$y[match(df_spiders$Group, df_mean.ord$Group)]

  # % variance explained, following the same convention used elsewhere in the
  # package for PCoA axes (positive eigenvalues only in the denominator);
  # falls back to a plain "PCoAn" label if eigenvalues aren't usable (e.g. all
  # non-positive), the same situation ord_labels() falls back to "DIMn" for.
  eig <- vegan::eigenvals(betadisper_obj)
  pos_sum <- sum(eig[eig > 0])
  axis_label <- function(i) {
    if (is.na(pos_sum) || pos_sum <= 0 || is.na(eig[i]) || eig[i] <= 0) {
      paste0("PCoA", i)
    } else {
      sprintf("PCoA%d (%.1f%%)", i, 100 * eig[i] / pos_sum)
    }
  }

  list(
    df_ord      = df_ord,
    df_mean.ord = df_mean.ord,
    df_spiders  = df_spiders,
    xlab        = axis_label(1),
    ylab        = axis_label(2)
  )
}

# --- Shared p-value formatter ---------------------------------------------
# Formats a p-value to 3 decimal places, rendering values below `accuracy` as
# "<0.001". Matches the formatting used for group-comparison p-values drawn by
# ggpubr::stat_compare_means() (which apply scales::label_pvalue() inline in
# their aes(label = ...) mapping), so every p-value in the package reads the
# same way.
.mbm_format_pval <- function(p, accuracy = 0.001) {
  scales::label_pvalue(accuracy = accuracy)(p)
}

# Same, as a full label: "p<0.001" below `accuracy`, "p=0.012" otherwise
# (never "p = <0.001").
.mbm_p_label <- function(p, accuracy = 0.001) {
  scales::pvalue(p, accuracy = accuracy, add_p = TRUE)
}

# --- Group-comparison layer for boxplots --------------------------------------
# Two-sample tests ("wilcox.test", "t.test") are run for every pair of boxes
# within each panel/facet, drawn with a bracket each and adjusted for
# multiple comparisons within the panel (p_adjust_method, "holm" by default;
# "none" to skip it). Global tests ("kruskal.test", "anova") give one p-value
# per panel. Before, stat_compare_means() with a two-sample test and more
# than two boxes drew every pairwise p-value on top of each other at the
# same spot, which read as a single (and misleading) p-value.
#   data, label.x, label.y: optional, for callers that add one layer per
#   panel with its own data subset (alpha_diversity_plot/alpha_hill_plot);
#   label.x/label.y only place the global-test label.
.mbm_stat_layer <- function(stat, p_adjust_method = "holm",
                            data = NULL, label.x = NULL, label.y = NULL) {
  if (stat %in% c("wilcox.test", "t.test")) {
    # Labels use the package's own p-value format ("p<0.001", "p=0.012"),
    # through stat_pwc()'s glue template, instead of rstatix's scientific
    # notation.
    p_col <- if (p_adjust_method == "none") "p" else "p.adj"
    ggpubr::stat_pwc(
      data            = data,
      method          = stat,
      p.adjust.method = p_adjust_method,
      p.adjust.by     = "panel",
      label           = paste0("{scales::pvalue(", p_col, ", 0.001, add_p = TRUE)}"),
      hide.ns         = p_col,
      size            = 0.4,
      label.size      = 3,
      family          = "serif",
      tip.length      = 0.01,
      step.increase   = 0.14,
      vjust           = -0.2
    )
  } else {
    ggpubr::stat_compare_means(
      method  = stat,
      mapping = ggplot2::aes(
        # kept short: ggpubr deparses this mapping and fails on multi-line code
        label = scales::pvalue(after_stat(p), 0.001, add_p = TRUE)
      ),
      size = 3.5, family = "serif", hide.ns = TRUE,
      data = data, label.x = label.x, label.y = label.y
    )
  }
}

# --- Multi-panel figures (patchwork) -------------------------------------------
# Joins a list of ggplots with patchwork instead of cowplot, so the result is
# still a modifiable plot: `p & theme(...)` changes every panel,
# `p[[2]] + labs(...)` a single one, `p + plot_annotation(...)` the figure.
# Identical legends are collected into one, and the A/B/C tags are patchwork
# tags (drawn outside each panel, like cowplot's labels were).
#   tags: character vector of panel tags (one per plot), or NULL for none.
.mbm_patchwork_grid <- function(plots, ncol, tags = NULL, title = NULL,
                                show_legend = TRUE, legend_position = "bottom",
                                tag_bold = TRUE, nrow = NULL, widths = NULL,
                                heights = NULL) {
  p <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow,
                             widths = widths, heights = heights) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title      = title,
      tag_levels = if (!is.null(tags)) list(tags[seq_along(plots)]) else NULL,
      theme      = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", family = "serif",
                                           size = 14, hjust = 0.5)
      )
    )
  p & ggplot2::theme(
    legend.position = if (show_legend) legend_position else "none",
    plot.tag = ggplot2::element_text(family = "serif", size = 14,
                                     face = if (tag_bold) "bold" else "plain")
  )
}

# Warns when a named `group_colors` matches none of the groups actually
# plotted (e.g. c("Roots" = ...) for boxes labelled "Rhizosphere_vs_Roots"):
# ggplot would otherwise silently draw every box grey.
.mbm_check_color_names <- function(group_colors, groups, arg = "group_colors") {
  if (!is.null(names(group_colors)) && !any(names(group_colors) %in% groups)) {
    warning("None of the names of `", arg, "` (", paste(names(group_colors), collapse = ", "),
            ") match the groups in the plot (", paste(groups, collapse = ", "),
            "), so they are drawn grey. Name the colors after these groups.",
            call. = FALSE)
  }
  invisible(NULL)
}

# Y-axis title for each row of a multi-panel grid: with more than one row the
# rows are short, so a long title (e.g. "Effective number of features") is
# wrapped onto two lines instead of running into the row above/below.
.mbm_row_axis_title <- function(title, n_rows, width = 18) {
  if (n_rows > 1 && is.character(title) && length(title) == 1 && nchar(title) > width) {
    paste(strwrap(title, width = width), collapse = "\n")
  } else {
    title
  }
}

# Adds A/B/C tags to the panels of a single faceted ggplot, as text in the
# top-left corner of each panel (in reading order: rows, then columns). The
# plot stays a normal ggplot, unlike the gtable-edited image used before.
.mbm_facet_tags <- function(p, tags, bold = TRUE) {
  layout <- ggplot2::ggplot_build(p)$layout$layout
  layout <- layout[order(layout$ROW, layout$COL), , drop = FALSE]
  facet_vars <- setdiff(names(layout), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y",
                                         "AXIS_X", "AXIS_Y", "COORD"))
  tag_df <- layout[, facet_vars, drop = FALSE]
  tag_df$.tag <- tags[seq_len(nrow(tag_df))]
  p + ggplot2::geom_text(
    data = tag_df,
    mapping = ggplot2::aes(x = -Inf, y = Inf, label = .tag),
    inherit.aes = FALSE, hjust = -0.4, vjust = 1.4,
    family = "serif", size = 4.5, fontface = if (bold) "bold" else "plain"
  )
}

# --- Renamed arguments ----------------------------------------------------------
# Several arguments were renamed so the same thing has the same name across
# the package (e.g. col_cond -> group_col). Functions take `...` and pass it
# here: an old name still works, with a warning saying what to use instead,
# and any other unknown argument is an error (so a typo isn't silently
# swallowed by `...`).
#   dots:    list(...) of the calling function.
#   renames: named character vector, c(old_name = "new_name").
#   fn:      the calling function's name, for the messages.
# Returns a named list of the values to use, keyed by the new names.
.mbm_renamed_args <- function(dots, renames, fn) {
  if (length(dots) == 0) return(list())
  nms <- names(dots)
  if (is.null(nms) || any(nms == "")) {
    stop("Unused unnamed argument(s) in ", fn, "().", call. = FALSE)
  }
  unknown <- setdiff(nms, names(renames))
  if (length(unknown) > 0) {
    stop("Unused argument(s) in ", fn, "(): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  out <- list()
  for (old in nms) {
    new <- renames[[old]]
    warning("In ", fn, "(), `", old, "` was renamed to `", new,
            "`; please use `", new, "` instead.", call. = FALSE)
    out[[new]] <- dots[[old]]
  }
  out
}

# --- Metadata / sample alignment ---------------------------------------------
# Matches metadata rows to the samples of a table by ID, never by position,
# since vegan::adonis2(), betadisper(), envfit(), ALDEx2 and randomForest all
# pair samples with metadata rows by position and would silently give wrong
# results on misaligned inputs.
#   sample_ids: sample IDs in the order the table/distance holds them.
#   metadata:   data frame whose first column holds the sample IDs.
# Stops if metadata has duplicated IDs or shares no sample with the table.
# Samples of the table that are not in metadata are left out, saying which
# ones with a message (subsetting metadata is a common way to keep only some
# groups); metadata rows for samples not in the table are dropped too. If
# the rows had to be reordered, says so with a message.
# Returns metadata in table order, for the samples present in both, after
# checking the IDs are identical row by row. Callers then keep only those
# samples of the table: table[, aligned[[1]]] (or the rows, or the dist).
.mbm_align_metadata <- function(sample_ids, metadata) {
  sample_ids <- trimws(as.character(sample_ids))
  meta_ids   <- trimws(as.character(metadata[[1]]))

  if (anyDuplicated(meta_ids)) {
    stop("Duplicated sample IDs in the first column of metadata: ",
         paste(head(unique(meta_ids[duplicated(meta_ids)])), collapse = ", "),
         call. = FALSE)
  }
  keep <- sample_ids[sample_ids %in% meta_ids]
  if (length(keep) == 0) {
    stop("None of the samples of the table are in the first column of metadata, ",
         "which must hold the sample IDs.\n",
         "Table samples (first 6): ", paste(head(sample_ids), collapse = ", "), "\n",
         "Metadata IDs (first 6): ", paste(head(meta_ids), collapse = ", "),
         call. = FALSE)
  }
  missing_ids <- setdiff(sample_ids, keep)
  if (length(missing_ids) > 0) {
    message(length(missing_ids), " sample(s) of the table are not in metadata ",
            "and were left out: ", paste(head(missing_ids), collapse = ", "),
            if (length(missing_ids) > 6) ", ..." else "")
  }

  if (!identical(meta_ids[meta_ids %in% keep], keep)) {
    message("metadata rows were reordered to match the sample order of the table.")
  }
  aligned <- metadata[match(keep, meta_ids), , drop = FALSE]
  aligned[[1]] <- keep

  if (!identical(as.character(aligned[[1]]), keep)) {
    stop("Internal error: metadata could not be aligned with the table.", call. = FALSE)
  }
  aligned
}

# --- Shared theme -----------------------------------------------------------
# Internal helper: unified ggplot2 theme for all MicroBioMeta plots.
# legend_position: passed through from each function's parameter.
# extra: additional theme() overrides specific to each plot type.
.mbm_theme <- function(legend_position = "bottom", extra = NULL) {
  base <- ggplot2::theme_bw(base_family = "serif") +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold",
                                               size = 14, color = "black"),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 12,
                                               color = "black"),
      axis.text.x      = ggplot2::element_text(size = 12, color = "black"),
      axis.text.y      = ggplot2::element_text(size = 12, color = "black"),
      axis.title.x     = ggplot2::element_text(size = 14, color = "black"),
      axis.title.y     = ggplot2::element_text(size = 14, color = "black"),
      legend.position  = legend_position,
      legend.title     = ggplot2::element_text(size = 14, face = "bold",
                                               color = "black"),
      legend.text      = ggplot2::element_text(size = 12, color = "black"),
      strip.text       = ggplot2::element_text(size = 12, face = "bold",
                                               color = "black"),
      strip.background = ggplot2::element_rect(fill = "white",
                                               color = "black"),
      panel.grid.minor = ggplot2::element_blank()
    )
  if (!is.null(extra)) base <- base + extra
  base
}

# --- Shared x-axis label element --------------------------------------------
# Internal helper: builds the axis.text.x element for a given rotation, taking
# care of the hjust/vjust that each angle needs (centered when horizontal,
# right/top-anchored once rotated). Every bar/boxplot function routes its
# `x_label_angle` parameter through here so rotation behaves identically
# package-wide. `colour = NA` is used to render the labels invisible while
# still reserving their space (shared-axis panels in the multi-panel grids).
.mbm_x_text <- function(angle = 0, colour = "black", size = 12) {
  ggplot2::element_text(
    angle = angle,
    hjust = if (angle == 0) 0.5 else 1,
    vjust = if (angle == 0) 0.5 else 1,
    size  = size,
    color = colour
  )
}

# --- Shared strip-text element ----------------------------------------------
# Internal helper: facet strip labels. Bold is opt-in (`strip_text_bold`) and
# defaults to FALSE across the package, so strips read as plain labels unless
# the user asks for emphasis.
.mbm_strip_text <- function(bold = FALSE, size = 12, angle = NULL, colour = "black") {
  ggplot2::element_text(
    size   = size,
    face   = if (bold) "bold" else "plain",
    family = "serif",
    color  = colour,
    angle  = angle
  )
}

# --- Colorblind-friendly DIVERGING palette presets --------------------------
# All built from Okabe-Ito colors. Pass the name to `diverging_palette` in
# heatmap/correlation functions. low <-> mid (white) <-> high.
.mbm_div_palettes <- list(
  "BuOr" = c("#0072B2", "white", "#E69F00"),   # blue <-> orange   (default)
  "BuVm" = c("#0072B2", "white", "#D55E00"),   # blue <-> vermillion
  "BuPk" = c("#0072B2", "white", "#CC79A7"),   # blue <-> pink
  "GnPk" = c("#009E73", "white", "#CC79A7"),   # green <-> pink
  "PuYl" = c("#440154", "white", "#FDE725")    # purple <-> yellow (viridis endpoints)
)

# --- Hill-number facet strip labels ------------------------------------------
# Renders q0/q1/q2 facet strips as bold "q=0", "q=1", "q=2" (plotmath), with
# only the "q" itself italic (standard math-variable convention) and the
# "=0"/"=1"/"=2" part upright.
# Use with `ggplot2::as_labeller(.mbm_q_labels, default = ggplot2::label_parsed)`.
.mbm_q_labels <- c(
  q0 = 'bolditalic(q)*bold("=0")',
  q1 = 'bolditalic(q)*bold("=1")',
  q2 = 'bolditalic(q)*bold("=2")'
)

# Plain-weight counterpart, used when `strip_text_bold = FALSE` (the package
# default). The plotmath above hard-codes bold(), so element_text(face =
# "plain") alone cannot unbold it - the expression itself has to change.
.mbm_q_labels_plain <- c(
  q0 = 'italic(q)*"=0"',
  q1 = 'italic(q)*"=1"',
  q2 = 'italic(q)*"=2"'
)

# Returns the q0/q1/q2 facet labeller at the requested weight.
.mbm_q_labeller <- function(bold = FALSE) {
  ggplot2::as_labeller(
    if (bold) .mbm_q_labels else .mbm_q_labels_plain,
    default = ggplot2::label_parsed
  )
}

# Re-sorts each "A_vs_B" string alphabetically (pmin/pmax, like the
# condition1_group/compar_condition1/2 columns built internally in
# beta_dissimilarity_plot()/beta_turnover_plot()), so a comparison_condition
# vector the caller wrote in whichever order ("Rhizosphere_vs_Bulk soil")
# still matches the internally-normalized column ("Bulk soil_vs_Rhizosphere")
# instead of requiring the caller to already know that internal convention.
.mbm_normalize_pair <- function(x) {
  parts <- strsplit(x, "_vs_", fixed = TRUE)
  vapply(parts, function(p) paste0(pmin(p[1], p[2]), "_vs_", pmax(p[1], p[2])), character(1))
}

# --- Colorblind-friendly palette (Okabe-Ito, 8 colors) ----------------------
# Safe for deuteranopia, protanopia and tritanopia.
# Used as the default qualitative palette across all functions except
# abundance barplots (which handle their own large palettes). Kept in the
# canonical Okabe-Ito publication order - some functions (e.g.
# aldex_heatmap_plot) index into specific positions of this exact ordering,
# so reordering it here would silently change those functions' colors too.
.mbm_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#000000")

# --- Dedicated 2-group color pair --------------------------------------------
# The orange/blue pairing used as THE 2-group comparison default across the
# package (e.g. aldex_volcano_plot's col_inf/col_sup, ancombc_plot's
# bar_colors). Deliberately a separate constant from `.mbm_colors` (rather
# than relying on `.mbm_colors[1:2]`) so it stays fixed regardless of how the
# 8-color qualitative palette above is ordered - use this explicitly whenever
# a function needs a default color for an exactly-2-group comparison.
# Matches `.mbm_colors[1:2]` exactly (same orange, same light blue) so a
# group's color doesn't change shade depending on whether it's being plotted
# alongside 1 or 3+ other groups (e.g. a "Rhizosphere" boxplot shouldn't be
# dark navy in a 2-group plot and light sky-blue in a 4-group one).
.mbm_colors_2group <- c("#E69F00", "#56B4E9")

# --- Second qualitative palette ("Safe", colorblind-friendly) ----------------
# From CARTOColors' "Safe" palette (based on Paul Tol's), reordered so its
# first colors (rose, indigo, olive) contrast with the Okabe-Ito orange, blue,
# green and yellow of `.mbm_colors`. Used for a second grouping variable shown
# next to the first one (e.g. treatments next to soil types in a heatmap).
.mbm_colors_safe <- c("#CC6677", "#332288", "#999933", "#882255",
                      "#44AA99", "#AA4499", "#661100", "#88CCEE",
                      "#DDCC77", "#117733", "#6699CC", "#888888")

# Named color vector for the values of a categorical variable: a named
# `colors` is matched by name (unmatched values grey, with a warning), an
# unnamed one is assigned in order.
.mbm_match_colors <- function(colors, values, arg = "colors") {
  if (is.null(names(colors))) {
    return(stats::setNames(colors[seq_along(values)], values))
  }
  .mbm_check_color_names(colors, values, arg)
  out <- stats::setNames(unname(colors[values]), values)
  out[is.na(out)] <- "grey70"
  out
}

# SILVA composite genus names of 3+ genera joined by "-" (e.g.
# "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium") become "<last genus>
# group" ("Rhizobium group") wherever they appear in a label, keeping any
# prefix or suffix ("ASV8_", "other ", " sp."). 2-genus names (e.g.
# "Escherichia-Shigella") and names whose parts aren't all genus-like
# (e.g. "Gitt-GS-136") are left as they are. Vectorized; NA stays NA.
.mbm_composite_genus <- function(x) {
  gsub("(?<![A-Za-z-])(?:[A-Z][a-z]+-){2,}([A-Z][a-z]+)(?![A-Za-z-])",
       "\\1 group", x, perl = TRUE)
}

# Row labels for heatmaps: composite genus names shortened (see
# .mbm_composite_genus(), if `composite`), and other labels longer than
# `max_length` cut with an ellipsis. Returns unique labels.
.mbm_shorten_labels <- function(x, max_length = 35, composite = TRUE) {
  if (composite) x <- .mbm_composite_genus(x)
  if (!is.null(max_length)) {
    long <- nchar(x) > max_length
    x[long] <- paste0(substr(x[long], 1, max_length - 1), "\u2026")
  }
  make.unique(x, sep = "_")
}

# Suppress R CMD check NOTEs for column names used in dplyr/ggplot2 NSE
# (no visible binding for global variable)
utils::globalVariables(c(
  # pipe
  "%>%",
  # common column names across functions
  ".",
  ".data",
  ":=",
  "taxonomy",
  "taxonomy2",
  "taxonomy_original",
  "OTUID",
  "OTU_ID",
  "SampleID",
  "Taxon",
  "Phylum",
  "taxon",
  "taxRank",
  # abundance / stats columns
  "Abundance",
  "RelativeAbundance",
  "MeanAbundance",
  "SumAbundance",
  "MeanAbund",
  "max_abund",
  "abun",
  "abund",
  # ordination columns
  "PC1",
  "PC2",
  "x",
  "y",
  "cntr.x",
  "cntr.y",
  "label",
  "ids",
  # beta diversity
  "TD_beta",
  "site1",
  "site2",
  "value",
  "name",
  "comparison",
  "grupo",
  # aldex columns
  "diff.btw",
  "effect",
  "wi.eBH",
  "wi.ep",
  "log_pvalue",
  ".p",
  ".tag",
  "significant",
  # global-test label in .mbm_stat_layer(): kept short, so after_stat() is
  # not namespaced there
  "after_stat",
  "p",
  "direction",
  # ancombc columns
  "lfc",
  "direct",
  # ratio plot
  "Ratio",
  "Dominant",
  "Condition",
  # random forest
  "MeanDecreaseGini",
  # corr plot
  "Correlation",
  "Environmental",
  # beta boxplot
  "compar_condition1",
  "compar_condition2",
  # taxonomy parsing
  "k", "p", "o", "f", "g", "s",
  # abundance heatmap
  "asv",
  "phylum",
  # misc
  "Group",
  "condition1_group",
  "across",
  "all_of",
  "where",
  "everything",
  "desc",
  "vars",
  "unit",
  # base R functions flagged
  "setNames",
  "reorder",
  "head",
  "as.dist",
  "as.formula",
  "cmdscale",
  "cor.test",
  "dist",
  "prcomp",
  "sd",
  "var",
  # alpha_decay_plot / beta_decay_plot corner-annotation layout columns
  ".grp_idx",
  "side",
  "slot",
  "hill",
  "x_pos",
  "y_pos",
  "hjust",
  "vjust",
  "geo_km",
  "similarity",
  "group",
  # alpha_diversity_plot facet formula
  "Index",
  # cca_rda_biplot vector labels
  "Variable",
  # random_forest_lollipop_plot row ids
  "row_id"
))
