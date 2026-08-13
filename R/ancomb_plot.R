#' ANCOMBC2 differential abundance bar/heatmap plot
#'
#' Runs ANCOMBC2 on a counts table and metadata, then visualizes differentially
#' abundant taxa. For two-group comparisons a bar plot of log-fold changes is
#' returned; for three or more groups a heatmap is returned; for a continuous
#' \code{col_cond} (e.g. \code{"dist_km"}) a bar plot of the effect size per
#' unit increase is returned instead, colored by the direction of the effect.
#'
#' @param table Data frame with taxa as rows and samples as columns. The last
#'   column must contain taxonomy strings (named "taxonomy", "Taxonomy",
#'   "taxon", "taxa", "Taxa", or "Taxon").
#' @param metadata Data frame with samples as rows. The first column must
#'   contain sample IDs that match the column names of \code{table}.
#' @param col_cond Character. Name of the column in \code{metadata} that
#'   defines the grouping variable. If this column is numeric (a continuous
#'   variable), it's treated as a covariate instead of a group: ANCOMBC2's
#'   \code{group}/structural-zero machinery (which requires discrete groups)
#'   is disabled, and the resulting plot shows the effect size per unit
#'   increase rather than a group-vs-group comparison.
#' @param tax_level Character or \code{NULL}. Taxonomic level to agglomerate
#'   to before running \code{ancombc2} (e.g. \code{"Genus"}, \code{"Family"}).
#'   Default \code{"Genus"}. Pass \code{NULL} to skip agglomeration and run
#'   ANCOMBC2 directly on the ASV/OTU-level table (rows of \code{table}, as-is).
#' @param prv_cut Numeric. Prevalence cut-off passed to \code{ancombc2}
#'   (default \code{0.1}). Lower values retain more taxa.
#' @param p_adj_method Character. Multiple-testing correction method passed to
#'   \code{ancombc2} (default \code{"holm"}). Use \code{"BH"} for a less
#'   strict correction when sample sizes are small.
#' @param formula Character. Right-hand side of the fixed-effects formula
#'   passed to \code{ancombc2} (e.g. \code{"group + age"}). If \code{NULL}
#'   (default) \code{col_cond} is used alone.
#' @param rand_formula Character. Random-effects formula passed to
#'   \code{ancombc2} for mixed models (e.g. \code{"~ 1 | subject_id"}).
#'   Default \code{NULL}.
#' @param ref_level Character. Reference level for \code{col_cond}. If
#'   \code{NULL} (default) the first factor level is used as reference.
#'   Use this to change which group appears as the baseline in comparisons
#'   (e.g. \code{ref_level = "P2"} to compare all other groups against P2).
#' @param diverging_palette Character. Name of the colorblind-friendly
#'   diverging palette used for the 3+-group heatmap's log-fold-change fill
#'   scale. One of \code{"BuOr"} (blue-orange, default), \code{"BuVm"}
#'   (blue-vermillion), \code{"BuPk"} (blue-pink), or \code{"GnPk"}
#'   (green-pink). Ignored for the 2-group / continuous bar plot, which uses
#'   \code{bar_colors} instead.
#' @param bar_colors Character vector of (at least) 2 colors used for the
#'   bar plot (2-group or continuous \code{col_cond}). First color is the
#'   "positive" direction (the non-reference group / increases with the
#'   variable); second color is the "negative" direction (the reference
#'   group / decreases with the variable). Default \code{c("#56B4E9",
#'   "#E69F00")} (colorblind-friendly blue/yellow). Ignored for the
#'   3+-group heatmap, which uses \code{diverging_palette} instead.
#' @param save_table Logical. If \code{TRUE}, saves the full ANCOMBC2 results
#'   table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"ancombc_results.txt"}.
#'
#' @return A \code{ggplot2} object: a bar plot (2 groups) or a heatmap
#'   (\eqn{\geq}3 groups).
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
#' # p_adj_method = "BH" (less strict than the "holm" default), and a lower
#' # prv_cut than the 0.1 default so rare-but-compartment-specific taxa
#' # (e.g. taxa found almost only in roots) aren't excluded before testing
#' ancombc_plot(
#'   table        = table,
#'   metadata     = metadata,
#'   col_cond     = "Location",
#'   prv_cut      = 0.01,
#'   p_adj_method = "BH"
#' )
#' }

ancombc_plot <- function(table,
                         metadata,
                         col_cond,
                         tax_level         = "Genus",
                         prv_cut           = 0.1,
                         p_adj_method      = "holm",
                         formula           = NULL,
                         rand_formula      = NULL,
                         ref_level         = NULL,
                         diverging_palette = "BuOr",
                         bar_colors        = c("#56B4E9", "#E69F00"),
                         save_table        = FALSE,
                         table_filename    = "ancombc_results.txt") {

  if (length(bar_colors) < 2)
    stop("`bar_colors` must have at least 2 colors.")

  # --- 0. package checks ---------------------------------------------------
  # ANCOMBC and phyloseq are on Bioconductor and can be auto-installed.
  for (pkg in c("ANCOMBC", "phyloseq")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    }
  }
  # qiime2R lives on GitHub only; it cannot be auto-installed via BiocManager.
  if (!requireNamespace("qiime2R", quietly = TRUE)) {
    stop(
      "Package 'qiime2R' is required but not installed.\n",
      "Install it with: remotes::install_github(\"jbisanz/qiime2R\")",
      call. = FALSE
    )
  }

  # --- 1. validate inputs --------------------------------------------------
  if (!col_cond %in% colnames(metadata))
    stop(paste("Column", col_cond, "not found in metadata."))

  if (!is.null(ref_level)) {
    if (!ref_level %in% metadata[[col_cond]])
      stop(paste("ref_level '", ref_level, "' not found in column '", col_cond, "'.", sep = ""))
    metadata[[col_cond]] <- relevel(factor(metadata[[col_cond]]), ref = ref_level)
  }

  tax_col <- grep("taxonomy|taxon|taxa", names(table), ignore.case = TRUE)
  if (length(tax_col) != 1)
    stop("Exactly one taxonomy column expected in 'table'.")

  names(table)[tax_col] <- "taxonomy"

  # --- 2. build phyloseq object ---------------------------------------------
  table_counts <- table %>% dplyr::select(-taxonomy)
  otumat  <- as.matrix(table_counts)

  taxa_df <- table %>%
    dplyr::select(Taxon = taxonomy) %>%
    tibble::rownames_to_column(var = "Feature.ID")
  taxmat  <- qiime2R::parse_taxonomy(taxa_df) %>% as.matrix()

  sample_id_col <- colnames(metadata)[1]

  # Build sample_data: use base R to avoid tibble's row-name restriction
  meta_df <- as.data.frame(metadata)
  rownames(meta_df) <- meta_df[[sample_id_col]]
  meta_df[[sample_id_col]] <- NULL

  OTU        <- phyloseq::otu_table(otumat, taxa_are_rows = TRUE)
  TAX        <- phyloseq::tax_table(taxmat)
  sampledata <- phyloseq::sample_data(meta_df)

  physeq      <- phyloseq::phyloseq(OTU, TAX, sampledata)
  physeq_filt <- phyloseq::prune_taxa(
    apply(phyloseq::otu_table(physeq), 1, var) > 0, physeq
  )
  dat <- physeq_filt

  # --- 3. run ANCOMBC2 -------------------------------------------------------
  fix_formula <- if (is.null(formula)) col_cond else formula

  # `group` (and the struc_zero/neg_lb machinery tied to it) only make sense
  # for a categorical grouping variable. If col_cond is numeric (e.g. a
  # continuous gradient like "dist_km"), treat it as a plain covariate:
  # no group, no structural-zero detection.
  is_continuous_cond <- is.numeric(meta_df[[col_cond]])
  if (is_continuous_cond) {
    message("`col_cond` ('", col_cond, "') is numeric; treating it as a ",
            "continuous covariate (group/struc_zero disabled).")
  }

  ancombc_res <- ANCOMBC::ancombc2(
    data          = dat,
    assay_name    = "counts",
    tax_level     = tax_level,
    fix_formula   = fix_formula,
    rand_formula  = rand_formula,
    p_adj_method  = p_adj_method,
    pseudo_sens   = TRUE,
    prv_cut       = prv_cut,
    lib_cut       = 1000,
    s0_perc       = 0.05,
    group         = if (is_continuous_cond) NULL else col_cond,
    struc_zero    = !is_continuous_cond,
    neg_lb        = !is_continuous_cond
  )

  res_prim <- ancombc_res$res

  if (save_table) {
    utils::write.table(res_prim, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # --- 4. parse formula terms to plot ----------------------------------------
  # Build list of terms: variables + interactions (e.g. "A*B" -> "A","B","A:B")
  fix_formula_str <- if (is.null(formula)) col_cond else formula
  simple_terms    <- trimws(unlist(strsplit(gsub("\\*", "+", fix_formula_str), "\\+")))
  simple_terms    <- unique(simple_terms[nchar(simple_terms) > 0])

  # add interaction terms if * was used
  all_terms <- simple_terms
  if (grepl("\\*", fix_formula_str)) {
    inter_vars  <- trimws(unlist(strsplit(
      grep("\\*", unlist(strsplit(fix_formula_str, "\\+")), value = TRUE),
      "\\*"
    )))
    inter_label <- paste(inter_vars, collapse = ":")
    all_terms   <- c(all_terms, inter_label)
  }

  # --- 5. helper functions --------------------------------------------------

  .make_barplot <- function(res, lfc_col, diff_col, se_col, term, meta_df) {
    cmp_group  <- sub(paste0("^lfc_", term), "", lfc_col)
    all_groups <- unique(as.character(meta_df[[term]]))
    ref_group  <- all_groups[!all_groups %in% cmp_group][1]
    if (is.na(ref_group)) ref_group <- "reference"

    df <- res %>%
      dplyr::select(taxon, dplyr::all_of(c(lfc_col, diff_col, se_col))) %>%
      dplyr::filter(.data[[diff_col]] %in% TRUE) %>%
      dplyr::arrange(.data[[lfc_col]]) %>%
      dplyr::mutate(
        taxon  = factor(taxon, levels = taxon),
        direct = factor(
          ifelse(.data[[lfc_col]] > 0, cmp_group, ref_group),
          levels = c(cmp_group, ref_group)
        )
      )

    if (nrow(df) == 0) return(NULL)

    fill_colors <- stats::setNames(bar_colors[1:2], c(cmp_group, ref_group))
    legend_labs <- stats::setNames(
      paste0("Higher in ", c(cmp_group, ref_group)),
      c(cmp_group, ref_group)
    )

    df %>%
      ggplot2::ggplot(ggplot2::aes(
        x = .data[[lfc_col]], y = taxon, fill = direct
      )) +
      ggplot2::geom_col(width = 0.7, color = "black", linewidth = 0.3) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          xmin = .data[[lfc_col]] - .data[[se_col]],
          xmax = .data[[lfc_col]] + .data[[se_col]]
        ),
        width = 0.3, color = "black", linewidth = 0.4, orientation = "y"
      ) +
      ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = fill_colors, name = NULL, labels = legend_labs) +
      ggplot2::labs(
        x     = paste0("Log fold change (", cmp_group, " vs ", ref_group, ")"),
        y     = NULL,
        title = paste0("Differential taxa by ", term)
      ) +
      .mbm_theme(
        legend_position = "bottom",
        extra = ggplot2::theme(
          axis.text.y        = ggplot2::element_text(size = 12, face = "italic",
                                                     color = "black"),
          panel.grid.major.y = ggplot2::element_blank()
        )
      )
  }

  # Continuous covariate (e.g. "dist_km"): there is no second group to
  # compare against, so bars are colored by the sign of the effect instead.
  .make_continuous_barplot <- function(res, lfc_col, diff_col, se_col, term) {
    df <- res %>%
      dplyr::select(taxon, dplyr::all_of(c(lfc_col, diff_col, se_col))) %>%
      dplyr::filter(.data[[diff_col]] %in% TRUE) %>%
      dplyr::arrange(.data[[lfc_col]]) %>%
      dplyr::mutate(
        taxon  = factor(taxon, levels = taxon),
        direct = factor(
          ifelse(.data[[lfc_col]] > 0, "Increases", "Decreases"),
          levels = c("Increases", "Decreases")
        )
      )

    if (nrow(df) == 0) return(NULL)

    fill_colors <- c(Increases = bar_colors[1], Decreases = bar_colors[2])
    legend_labs <- c(
      Increases = paste0("Increases with ", term),
      Decreases = paste0("Decreases with ", term)
    )

    df %>%
      ggplot2::ggplot(ggplot2::aes(
        x = .data[[lfc_col]], y = taxon, fill = direct
      )) +
      ggplot2::geom_col(width = 0.7, color = "black", linewidth = 0.3) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          xmin = .data[[lfc_col]] - .data[[se_col]],
          xmax = .data[[lfc_col]] + .data[[se_col]]
        ),
        width = 0.3, color = "black", linewidth = 0.4, orientation = "y"
      ) +
      ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = fill_colors, name = NULL, labels = legend_labs) +
      ggplot2::labs(
        x     = paste0("Log fold change (", term, ")"),
        y     = NULL,
        title = paste0("Differential taxa by ", term)
      ) +
      .mbm_theme(
        legend_position = "bottom",
        extra = ggplot2::theme(
          axis.text.y        = ggplot2::element_text(size = 12, face = "italic",
                                                     color = "black"),
          panel.grid.major.y = ggplot2::element_blank()
        )
      )
  }

  .make_heatmap <- function(res, lfc_cols, diff_cols, term) {
    df_cond <- res %>%
      dplyr::select(taxon, dplyr::all_of(c(diff_cols, lfc_cols)))

    df_sig <- df_cond %>%
      dplyr::filter(dplyr::if_any(dplyr::all_of(diff_cols), ~ . %in% TRUE))

    if (nrow(df_sig) == 0) return(NULL)

    df_long <- df_sig %>%
      dplyr::select(taxon, dplyr::all_of(lfc_cols)) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(lfc_cols), ~ round(.x, 2))) %>%
      tidyr::pivot_longer(
        cols      = dplyr::all_of(lfc_cols),
        names_to  = "comparison",
        values_to = "lfc"
      ) %>%
      dplyr::mutate(
        comparison = sub(paste0("^lfc_", term), "", comparison)
      )

    lo  <- floor(min(df_long$lfc,  na.rm = TRUE))
    up  <- ceiling(max(df_long$lfc, na.rm = TRUE))
    mid <- (lo + up) / 2

    df_long %>%
      ggplot2::ggplot(ggplot2::aes(x = comparison, y = taxon, fill = lfc)) +
      ggplot2::geom_tile(color = "black") +
      ggplot2::scale_fill_gradient2(
        low      = .mbm_div_palettes[[diverging_palette]][1],
        mid      = "white",
        high     = .mbm_div_palettes[[diverging_palette]][3],
        na.value = "white", midpoint = mid, limit = c(lo, up), name = "LFC"
      ) +
      ggplot2::geom_text(ggplot2::aes(label = lfc), size = 3.5) +
      ggplot2::labs(
        x     = NULL,
        y     = NULL,
        title = paste("Log fold changes -", term)
      ) +
      .mbm_theme(
        legend_position = "right",
        extra = ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 12, face = "italic",
                                              color = "black"),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                              size = 12, color = "black")
        )
      )
  }

  # --- 6. build one plot per term -------------------------------------------
  plots <- list()

  for (term in all_terms) {
    lfc_t  <- grep(paste0("^lfc_",  term), names(res_prim), value = TRUE)
    diff_t <- grep(paste0("^diff_", term), names(res_prim), value = TRUE)
    se_t   <- grep(paste0("^se_",   term), names(res_prim), value = TRUE)

    if (length(lfc_t) == 0) next   # term not found in results (e.g. naming mismatch) -> skip

    if (length(lfc_t) == 1) {
      p <- if (is.numeric(meta_df[[term]])) {
        .make_continuous_barplot(res_prim, lfc_t, diff_t, se_t, term)
      } else {
        .make_barplot(res_prim, lfc_t, diff_t, se_t, term, meta_df)
      }
    } else {
      p <- .make_heatmap(res_prim, lfc_t, diff_t, term)
    }

    if (!is.null(p)) {
      plots[[term]] <- p
      # also print diagnostic
      n_sig <- if (length(lfc_t) == 1)
        sum(res_prim[[diff_t]] %in% TRUE)
      else
        sum(apply(res_prim[, diff_t, drop = FALSE], 1,
                  function(r) any(r %in% TRUE)))
      message(sprintf("Term '%s': %d significant taxa", term, n_sig))
    } else {
      message(sprintf("Term '%s': no significant taxa found - skipped", term))
    }
  }

  if (length(plots) == 0)
    stop(paste0(
      "No significantly different taxa found for any term.\n",
      "Try p_adj_method = 'BH' or lower prv_cut (current = ", prv_cut, ")."
    ))

  # return single plot directly, or named list if multiple terms
  if (length(plots) == 1) return(plots[[1]])
  return(plots)
}
