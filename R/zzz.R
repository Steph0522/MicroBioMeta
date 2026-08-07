# ── Shared p-value formatter ───────────────────────────────────────────────────
# Formats a p-value to 3 decimal places, rendering values below `accuracy` as
# "<0.001". Matches the formatting used for group-comparison p-values drawn by
# ggpubr::stat_compare_means() (which apply scales::label_pvalue() inline in
# their aes(label = ...) mapping), so every p-value in the package reads the
# same way.
.mbm_format_pval <- function(p, accuracy = 0.001) {
  scales::label_pvalue(accuracy = accuracy)(p)
}

# ── Shared theme ──────────────────────────────────────────────────────────────
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

# ── Colorblind-friendly DIVERGING palette presets ─────────────────────────────
# All built from Okabe-Ito colors. Pass the name to `diverging_palette` in
# heatmap/correlation functions. low ↔ mid (white) ↔ high.
.mbm_div_palettes <- list(
  "BuOr" = c("#0072B2", "white", "#E69F00"),   # blue <-> orange   (default)
  "BuVm" = c("#0072B2", "white", "#D55E00"),   # blue <-> vermillion
  "BuPk" = c("#0072B2", "white", "#CC79A7"),   # blue <-> pink
  "GnPk" = c("#009E73", "white", "#CC79A7")    # green <-> pink
)

# ── Hill-number facet strip labels ────────────────────────────────────────────
# Renders q0/q1/q2 facet strips as bold italic "q=0", "q=1", "q=2" (plotmath).
# Use with `ggplot2::as_labeller(.mbm_q_labels, default = ggplot2::label_parsed)`.
.mbm_q_labels <- c(
  q0 = 'bolditalic("q=0")',
  q1 = 'bolditalic("q=1")',
  q2 = 'bolditalic("q=2")'
)

# ── Colorblind-friendly palette (Okabe-Ito, 8 colors) ─────────────────────────
# Safe for deuteranopia, protanopia and tritanopia.
# Used as the default qualitative palette across all functions except
# abundance barplots (which handle their own large palettes).
.mbm_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#000000")

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
  "install.packages",
  "as.dist",
  "as.formula",
  "cmdscale",
  "cor.test",
  "dist",
  "prcomp",
  "sd",
  "var"
))
