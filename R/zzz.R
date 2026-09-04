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
  "significant",
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
  "var"
))
