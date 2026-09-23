#' Box plot of beta diversity
#'
#' @param table table Data frame where columns are samples and rows are ASVs or taxa.
#' @param metadata  A data frame with sample metadata. The first column must match sample names in "table".
#' @param comparison_condition1 Vector of \code{"A_vs_B"} group-pair labels to
#'   keep (matched against \code{condition1.x}/\code{condition1.y}). Matching
#'   ignores order - listing \code{"A_vs_B"} also matches pairs the pairwise
#'   self-join happened to record as \code{"B_vs_A"}, so each pair only needs
#'   to be listed once. Each matched pair becomes one x-axis/fill group,
#'   labelled exactly as written (e.g. \code{"Rhizosphere_vs_Roots"}), same
#'   as \code{condition1_group} in \code{beta_dissimilarity_plot()}.
#' @param comparison_condition2 Optional. Same as \code{comparison_condition1},
#'   for the \code{condition2.x}/\code{condition2.y} pair (e.g. \code{"TC_vs_TC"}
#'   to keep only within-group pairs of a given treatment). Only used to filter
#'   which pairs are kept - not shown on the plot, so term order doesn't
#'   matter here. Default \code{NULL}: every pair matching
#'   \code{comparison_condition1} is kept regardless of this second condition.
#' @param condition1.x Metadata column (as it reads after the pairwise
#'   self-join, e.g. \code{"Type_of_soil.x"}) used - together with
#'   \code{condition1.y} - to build \code{comparison_condition1}'s group
#'   pairs.
#' @param condition1.y Metadata column (self-join suffix \code{.y}) paired
#'   with \code{condition1.x}.
#' @param condition2.x Optional. Metadata column (self-join suffix \code{.x})
#'   used with \code{condition2.y} to build \code{comparison_condition2}'s
#'   group pairs. Only needed when \code{comparison_condition2} is set.
#' @param condition2.y Optional. Metadata column (self-join suffix \code{.y})
#'   paired with \code{condition2.x}.
#' @param facet_by Optional. A metadata column (as it reads after the
#'   pairwise self-join, e.g. \code{"Treatment.x"}) to facet the plot by -
#'   purely visual, like \code{beta_dissimilarity_plot()}'s
#'   \code{condition2_col}: every pair is shown, faceted by this column's
#'   value, with no filtering. Unrelated to \code{comparison_condition2},
#'   which filters instead of faceting and doesn't require picking a side
#'   (\code{.x} vs \code{.y}) since it compares both.
#' @param color_facets_x Optional color vector for facet strips. Defaults to
#'   a neutral \code{"grey85"} background for each facet, like
#'   \code{beta_dissimilarity_plot()}'s \code{facet_colors}.
#' @param color_axis_x Optional named color vector for x-axis groups.
#'   Defaults to the package's colorblind-friendly Okabe-Ito palette
#'   (\code{.mbm_colors}, orange/blue first), like
#'   \code{beta_dissimilarity_plot()}'s \code{group_colors}.
#' @param x_axis_title Title for the x-axis. Default \code{"Section"}.
#' @param show_x_labels Logical. If \code{TRUE}, x-axis tick labels are shown.
#'   Default \code{FALSE}, since the same groups are already named in the legend.
#' @param x_label_angle Numeric. Rotation (in degrees) of the x-axis tick labels
#'   when \code{show_x_labels = TRUE}. Default \code{0} (horizontal).
#' @param strip_text_bold Logical. If \code{TRUE}, facet strip labels are bold.
#'   Default \code{FALSE} (plain).
#' @param strip_text_color Color of the top (x) facet strip labels, which sit on
#'   the \code{color_facets_x} backgrounds. Default \code{"white"}, unless
#'   \code{color_facets_x} is left at its own light \code{"grey85"} default,
#'   in which case this defaults to \code{"black"} instead so it stays legible.
#' @param aspect_ratio Numeric. Aspect ratio (height/width) of each panel.
#'   Default \code{NULL} (automatic).
#' @param stat Character or \code{NULL}. Statistical test to compare groups
#'   within each panel, passed to \code{ggpubr::stat_compare_means()} (e.g.
#'   \code{"wilcox.test"}, \code{"kruskal.test"}, \code{"anova"}). Default
#'   \code{NULL} (no test shown).
#' @param save_table Logical. If \code{TRUE}, saves the underlying turnover
#'   table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"betadiv_turnover.txt"}.
#'
#' @return A ggplot2 figure with beta diversity partitions across conditions.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "OTUID"
#'
#' # comparison_condition1/2 match group pairs regardless of order, so listing
#' # "Rhizosphere_vs_Roots" also catches pairs the self-join recorded the
#' # other way around ("Roots_vs_Rhizosphere") - no need to list both.
#' beta_turnover_plot(
#'   table                 = table,
#'   metadata              = metadata,
#'   comparison_condition1 = c("Rhizosphere_vs_Roots", "Rhizosphere_vs_Rhizosphere"),
#'   comparison_condition2 = c("Control_vs_Control", "Moderate_drought_vs_Moderate_drought"),
#'   condition1.x          = "Location.x",
#'   condition1.y          = "Location.y",
#'   condition2.x          = "Treatment.x",
#'   condition2.y          = "Treatment.y",
#'   color_facets_x        = c("#5D478B", "#8B668B"),
#'   color_axis_x          = c("Roots" = "#56B4E9", "Rhizosphere" = "#E69F00")
#' )

beta_turnover_plot <- function(table, 
                      metadata, 
                      comparison_condition1,
                      comparison_condition2 = NULL,
                      condition1.x,
                      condition1.y,
                      condition2.x = NULL,
                      condition2.y = NULL,
                      facet_by = NULL,
                      color_facets_x = NULL,
                      color_axis_x = NULL,
                      x_axis_title = "Section",
                      show_x_labels = FALSE,
                      x_label_angle = 0,
                      strip_text_bold = FALSE,
                      strip_text_color = "white",
                      aspect_ratio = NULL,
                      stat = NULL,
                      save_table = FALSE,
                      table_filename = "betadiv_turnover.txt") {

  # Treat metadata's first column as the sample ID regardless of its original name
  colnames(metadata)[1] <- "OTUID"

  # Drop a taxonomy column if present, so it isn't treated as a sample
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if (length(tax_col) == 1) table <- table[, -tax_col, drop = FALSE]

  # Paso 1: Filtrar y transponer tabla OTU
  otu_filter <- table %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric))) != 0)
  
  if (nrow(otu_filter) == 0) stop("OTU table is empty after filtering rows with sum=0.")

  otu_filter_t <- as.data.frame(t(otu_filter))

  # hill_taxa_parti_pairwise() (below) computes every pairwise sample
  # distance, which grows quadratically with sample count, and it's run once
  # per q (0/1/2) - so restrict to only the samples whose groups are actually
  # named in comparison_condition1/2 before that, instead of computing pairs
  # for every sample in the table and discarding most of them afterwards in
  # the compar_condition1/2 filter below.
  raw_col1 <- sub("\\.[xy]$", "", condition1.x)
  groups1_needed <- unique(unlist(strsplit(comparison_condition1, "_vs_")))
  keep_samples <- metadata$OTUID[metadata[[raw_col1]] %in% groups1_needed]

  if (!is.null(comparison_condition2)) {
    raw_col2 <- sub("\\.[xy]$", "", condition2.x)
    groups2_needed <- unique(unlist(strsplit(comparison_condition2, "_vs_")))
    keep_samples <- intersect(keep_samples, metadata$OTUID[metadata[[raw_col2]] %in% groups2_needed])
  }

  otu_filter_t <- otu_filter_t[rownames(otu_filter_t) %in% keep_samples, , drop = FALSE]
  if (nrow(otu_filter_t) == 0) {
    stop("No samples left after restricting to the groups named in comparison_condition1/2.")
  }

  # Paso 2: Calcular diversidad beta (Hill numbers)
  beta_q_list <- list()
  for (q in c(0, 1, 2)) {
    beta_res <- tryCatch({
      # hill_taxa_parti_pairwise() prints a raw pairwise-comparison progress
      # bar straight to stdout (not via message()/warning()), which
      # suppressMessages() can't catch - capture.output() discards it while
      # still returning the function's actual result.
      result <- NULL
      utils::capture.output(
        result <- hillR::hill_taxa_parti_pairwise(comm = otu_filter_t, q = q)
      )
      result %>%
        dplyr::mutate(Recambio = TD_beta - 1, q = q)
    }, error = function(e) {
      warning("Could not compute hill_taxa_parti_pairwise with q = ", q, ": ", e$message)
      NULL
    })
    if (!is.null(beta_res)) beta_q_list[[as.character(q)]] <- beta_res
  }
  
  beta_total <- dplyr::bind_rows(beta_q_list)
  
  if (nrow(beta_total) == 0) stop("Could not compute beta partitions (empty table).")


  # Paso 3: Unir con metadata
  beta_formato <- beta_total %>%
    dplyr::inner_join(metadata, by = c("site1" = "OTUID")) %>%
    dplyr::inner_join(metadata, by = c("site2" = "OTUID"))
  
  if (nrow(beta_formato) == 0) stop("Could not join beta_total with metadata (empty table).")


  # Paso 4: Crear comparaciones y filtrar
  # compar_condition1/2 (used only for matching) are order-independent
  # (pmin/pmax on the character values - as.character() first, because these
  # columns are typically factors, and pmin/pmax's underlying `>`/`<` on
  # unordered factors always returns NA - the exact bug this fixed in
  # beta_dissimilarity_plot). This means the caller only has to list a pair
  # once (e.g. "A_vs_B", not also "B_vs_A"), matched regardless of which
  # sample happened to land as site1 vs site2 in the pairwise self-join.
  c1_x <- as.character(beta_formato[[condition1.x]])
  c1_y <- as.character(beta_formato[[condition1.y]])
  beta_formato$compar_condition1 <- paste0(pmin(c1_x, c1_y), "_vs_", pmax(c1_x, c1_y))

  # comparison_condition2/condition2.x/condition2.y are optional - when left
  # NULL, every pair matching comparison_condition1 is kept regardless of
  # the second condition (e.g. Treatment), instead of requiring one.
  has_condition2 <- !is.null(comparison_condition2)
  if (has_condition2) {
    c2_x <- as.character(beta_formato[[condition2.x]])
    c2_y <- as.character(beta_formato[[condition2.y]])
    beta_formato$compar_condition2 <- paste0(pmin(c2_x, c2_y), "_vs_", pmax(c2_x, c2_y))
  }

  # Combined into one "A_vs_B" label - both the x-axis and fill group, same
  # as condition1_group in beta_dissimilarity_plot() - instead of splitting
  # the reference group off into its own facet. Uses the pair exactly as the
  # caller wrote it (not the alphabetically-normalized pair_norm), so e.g.
  # "Rhizosphere_vs_Roots" stays that way instead of becoming "Roots_vs_..." .
  pair_norm     <- .mbm_normalize_pair(comparison_condition1)
  label_by_norm <- stats::setNames(comparison_condition1, pair_norm)
  # A factor, ordered by the order comparison_condition1 was written in -
  # left as plain character, ggplot/ggpubr would order the x-axis/legend
  # alphabetically instead (or by first appearance in the data, neither of
  # which the caller controls), same fix as condition1_group's factor in
  # beta_dissimilarity_plot().
  beta_formato$.compar_label <- factor(
    unname(label_by_norm[beta_formato$compar_condition1]),
    levels = unique(comparison_condition1)
  )

  # facet_by is purely visual (no filtering) - every pair is kept and simply
  # faceted by this column's value, same role as condition2_col in
  # beta_dissimilarity_plot().
  has_facet_by <- !is.null(facet_by)
  if (has_facet_by) {
    beta_formato$.facet2_col <- as.character(beta_formato[[facet_by]])
  }

  beta_final <- beta_formato %>%
    dplyr::filter(compar_condition1 %in% pair_norm)

  if (has_condition2) {
    beta_final <- beta_final %>%
      dplyr::filter(compar_condition2 %in% .mbm_normalize_pair(comparison_condition2))
  }

  beta_final <- beta_final %>%
    dplyr::mutate(orden = dplyr::case_when(
      q == 0 ~ "q0", q == 1 ~ "q1", q == 2 ~ "q2", TRUE ~ as.character(q)
    ))
  
  if (nrow(beta_final) == 0) stop("After filtering comparisons, the table is empty.")

  # --- Default palettes if missing (same pattern as beta_dissimilarity_plot) ---
  if (is.null(color_axis_x)) {
    n_groups <- length(unique(beta_final$.compar_label))
    color_axis_x <- if (n_groups == 2) .mbm_colors_2group else rep_len(.mbm_colors, n_groups)
    names(color_axis_x) <- unique(beta_final$.compar_label)
  }
  if (has_facet_by && is.null(color_facets_x)) {
    # strip_text_color's own "white" default assumes the caller's (often
    # dark) color_facets_x - illegible against this light default, so switch
    # to black here instead, unless the caller asked for white explicitly.
    if (missing(strip_text_color)) strip_text_color <- "black"
    n_facets <- length(unique(beta_final$.facet2_col))
    color_facets_x <- rep("grey85", n_facets)
  }

  if (save_table) {
    utils::write.table(beta_final, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message("Table saved as: ", table_filename)
  }

  # Paso 5: Crear grafico
  q_labeller <- .mbm_q_labeller(strip_text_bold)

  # The x-axis groups are already named in the legend, so their tick labels are
  # hidden by default; `show_x_labels = TRUE` brings them back (rotated by
  # `x_label_angle`, as in the other bar/boxplot functions).
  x_text  <- if (show_x_labels) .mbm_x_text(x_label_angle) else ggplot2::element_blank()
  x_ticks <- if (show_x_labels) ggplot2::element_line(colour = "black") else ggplot2::element_blank()

  # facet_by (purely visual - see its docs) adds a facet column, e.g.
  # Treatment; unrelated to comparison_condition2, which only filters and
  # stays invisible here. Without it, the only facet dimension is q (orden).
  facet_spec <- if (has_facet_by) {
    ggh4x::facet_grid2(orden ~ .facet2_col,
                       scales = "free_x",
                       labeller = ggplot2::labeller(orden = q_labeller),
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x)))
  } else {
    ggh4x::facet_grid2(orden ~ .,
                       scales = "free_x",
                       labeller = ggplot2::labeller(orden = q_labeller))
  }

  figura <- beta_final %>%
    ggpubr::ggboxplot(x = ".compar_label", y = "Recambio", fill = ".compar_label") +
    # Generic "features" instead of "ASVs" - the table can just as well hold
    # OTUs, species, or any other feature type.
    ggplot2::ylab("Proportion of feature turnover") +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggplot2::labs(fill = "Comparison") +
    facet_spec +
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        panel.grid   = ggplot2::element_blank(),
        axis.text.x  = x_text,
        axis.ticks.x = x_ticks,
        # Plain weight, size 14, matching .mbm_theme()'s own axis.title.y
        # default (and every other plot function's axis titles) - only the
        # margin needed overriding here.
        axis.title.y = ggplot2::element_text(size = 14, color = "black",
                                             margin = ggplot2::margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        # x strips sit on the user-supplied `color_facets_x` backgrounds (often
        # dark), so their text color is exposed separately from the y strips,
        # which keep the package's plain grey-strip look (matching the
        # q0/q1/q2 strips in alpha_hill_plot/alpha_diversity_plot) instead of
        # `.mbm_theme()`'s white default.
        strip.background.y = ggplot2::element_rect(fill = "grey", color = "black"),
        strip.text.x = .mbm_strip_text(strip_text_bold, colour = strip_text_color),
        strip.text.y = .mbm_strip_text(strip_text_bold),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    ) +
    ggplot2::xlab(x_axis_title)

  if (!is.null(aspect_ratio)) {
    figura <- figura + ggplot2::theme(aspect.ratio = aspect_ratio)
  }

  if (!is.null(stat)) {
    # stat_compare_means() places its label a fixed fraction above the data's
    # max, which the default 5% top expansion doesn't leave room for - the
    # label gets clipped by the panel border. Widen the top expansion instead
    # of leaving that headroom out only for this stat-annotated case.
    figura <- figura +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
      ggpubr::stat_compare_means(
        method = stat,
        mapping = ggplot2::aes(
          label = paste0("p = ", scales::label_pvalue(accuracy = 0.001)(ggplot2::after_stat(p)))
        ),
        size = 3.5, family = "serif", hide.ns = TRUE
      )
  }

  # Returns the ggplot object directly (not wrapped in a list) so it drops
  # straight into cowplot::plot_grid()/patchwork etc. like every other plot
  # function in the package (beta_dissimilarity_plot(), alpha_hill_plot(), ...).
  figura
}
