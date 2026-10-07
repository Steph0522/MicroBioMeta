#' Beta Diversity Boxplot
#'
#' This function calculates beta diversity (shared species, turnover, or nestedness)
#' and generates a ggplot2 boxplot with facets.
#'
#' @param table A data frame with taxa in rows and samples in columns. 
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. 
#' Must include a `SAMPLEID` column matching sample names in `table`.
#' @param condition1_col Column name in metadata for the first condition.
#' @param condition2_col Optional column name in metadata for the second condition, used as facet. 
#' @param facet_colors Optional vector of colors for facet strips. Defaults to a neutral \code{"grey85"}.
#' @param x_axis_title Title for the x-axis.
#' @param y_axis_title Title for the y-axis. Default \code{NULL}: it is built automatically.
#' @param partition Type of beta diversity to compute: `"shared"` (default),
#'   `"turnover"`, or `"nestedness"`. Case-insensitive.
#' @param family Dissimilarity family for the turnover/nestedness partition:
#'   `"sorensen"` (default) or `"jaccard"`. Case-insensitive.
#' @param stat Character or \code{NULL}. Statistical test to compare the
#'   boxes within each panel/facet. \code{"wilcox.test"} or \code{"t.test"}
#'   compare every pair, (\code{ggpubr::stat_pwc()}); \code{"kruskal.test"} or \code{"anova"}
#'   give one global p-value for more than 3 factors
#'   Default \code{NULL} (no test shown).
#' @param p_adjust_method Multiple-comparison correction for the pairwise
#'   any method of \code{stats::p.adjust()}. Default
#'   \code{"holm"}; \code{"none"} shows the raw p-values.
#' @param show_x_labels Logical. If \code{TRUE} (default), x-axis tick labels
#'   are shown. 
#' @param x_label_angle Numeric. Rotation (in degrees) of the x-axis tick
#'   labels when \code{show_x_labels = TRUE}. Default \code{0} (horizontal);
#'   use e.g. \code{45} or \code{90} when comparison names are long enough to
#'   overlap.
#' @param strip_text_bold Logical. If \code{TRUE}, facet strip labels are bold. Default \code{FALSE}.
#' @param strip_text_color Color of the facet strip labels, which sit on the
#'   \code{facet_colors} backgrounds. Default \code{"black"} (the default
#'   \code{facet_colors} is a light \code{"grey85"}; pass \code{"white"} if
#'   you supply darker \code{facet_colors}).
#' @param aspect_ratio Numeric. Sets the aspect ratio (height/width) of each panel. Default \code{NULL} (automatic).
#' @param save_table Logical. If \code{TRUE}, saves the beta diversity table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table. Default \code{"betadiv_table.txt"}.
#'
#' @return A ggplot2 figure object.
#' @export
#' 
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' beta_dissimilarity_plot(
#'   table                = table,
#'   metadata             = metadata,
#'   comparison_condition1 = c("Rhizosphere_vs_Roots"),
#'   condition1_col       = "Location",
#'   condition2_col       = "Treatment",
#'   x_axis_title         = "Samples",
#'   partition            = "shared",
#'   family               = "sorensen"
#' )

beta_dissimilarity_plot <- function(
    table, metadata,
    comparison_condition1 = NULL,
    condition1_col,
    condition2_col = NULL,
    facet_colors = NULL,
    group_colors = NULL,
    x_axis_title = "Condition",
    y_axis_title = NULL,
    partition = c("shared","turnover","nestedness"),
    family = c("sorensen","jaccard"),
    stat = NULL,
    p_adjust_method = "holm",
    show_x_labels = TRUE,
    x_label_angle = 0,
    strip_text_bold = FALSE,
    strip_text_color = "black",
    aspect_ratio = NULL,
    save_table = FALSE,
    table_filename = "betadiv_table.txt"
) {
  
#prepare data
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) == 1) table <- table[ , -tax_col]
  
  partition <- match.arg(tolower(partition), c("shared", "turnover", "nestedness"))
  family <- match.arg(tolower(family), c("sorensen", "jaccard"))
  
  table[table>0] <- 1
  table_filtered <- table[rowSums(table)>1,]
  table_t <- as.data.frame(t(table_filtered))
  table_t[] <- lapply(table_t, as.numeric)
  rownames(table_t) <- colnames(table_filtered)
  
  if(partition=="shared") {
    beta_mat <- betapart::betapart.core(table_t)$shared
  } else {
    beta_pair <- betapart::beta.pair(table_t, index.family = family)
    if(family=="jaccard") {
      beta_mat <- if(partition=="turnover") beta_pair$beta.jtu else beta_pair$beta.jne
    } else { # sorensen
      beta_mat <- if(partition=="turnover") beta_pair$beta.sim else beta_pair$beta.sne
    }
  }

  beta_mat <- as.matrix(beta_mat)
  idx <- which(lower.tri(beta_mat), arr.ind = TRUE)
  beta_df <- data.frame(site1 = rownames(beta_mat)[idx[, "row"]],
                        site2 = colnames(beta_mat)[idx[, "col"]],
                        value = beta_mat[idx],
                        stringsAsFactors = FALSE)
  beta_df <- dplyr::filter(beta_df, !is.na(value))
  
  beta_df <- dplyr::left_join(beta_df, metadata, by=c("site1"=names(metadata)[1])) %>%
    dplyr::left_join(metadata, by=c("site2"=names(metadata)[1]))
  
  # --- Construct comparison labels ---
  c1_x <- as.character(beta_df[[paste0(condition1_col,".x")]])
  c1_y <- as.character(beta_df[[paste0(condition1_col,".y")]])
  beta_df <- dplyr::mutate(beta_df,
                           condition1_group = paste0(pmin(c1_x, c1_y),
                                                     "_vs_",
                                                     pmax(c1_x, c1_y)))
  
  if(!is.null(comparison_condition1)) {
    norm1 <- .mbm_normalize_pair(comparison_condition1)
    beta_df <- dplyr::filter(beta_df, condition1_group %in% norm1)
    beta_df$condition1_group <- factor(beta_df$condition1_group, levels = unique(norm1))
  }

  if (!is.null(condition2_col)) {
    c2_x <- as.character(beta_df[[paste0(condition2_col, ".x")]])
    c2_y <- as.character(beta_df[[paste0(condition2_col, ".y")]])
    beta_df <- beta_df[!is.na(c2_x) & c2_x == c2_y, , drop = FALSE]
    if (nrow(beta_df) == 0) {
      stop("No pairs of samples share the same value of condition2_col '",
           condition2_col, "'.")
    }
  }

  
  if (save_table) {
    utils::write.table(
      beta_df,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message("Table saved as: ", table_filename)
  }
  n_groups <- length(unique(beta_df$condition1_group))
  if(is.null(group_colors)) {
    group_colors <- if (n_groups == 2) .mbm_colors_2group else rep_len(.mbm_colors, n_groups)
  }
  if(is.null(names(group_colors))) names(group_colors) <- unique(beta_df$condition1_group)
  .mbm_check_color_names(group_colors, as.character(unique(beta_df$condition1_group)))

  if(!is.null(condition2_col)) {
    n_facets <- length(unique(beta_df[[paste0(condition2_col,".x")]]))
    if(is.null(facet_colors)) facet_colors <- rep("grey85", n_facets)
  }
  

  x_text  <- if (show_x_labels) .mbm_x_text(x_label_angle) else ggplot2::element_blank()
  x_ticks <- if (show_x_labels) ggplot2::element_line(colour = "black") else ggplot2::element_blank()

  base_theme <- .mbm_theme(
    legend_position = "right",
    extra = ggplot2::theme(
      panel.grid   = ggplot2::element_blank(),
      strip.text   = .mbm_strip_text(strip_text_bold, colour = strip_text_color),
      axis.text.x  = x_text,
      axis.ticks.x = x_ticks
    )
  )
  if (!is.null(aspect_ratio)) {
    base_theme <- base_theme + ggplot2::theme(aspect.ratio = aspect_ratio)
  }

  # --- Create figure ---
  # "shared" alone doesn't say shared *what* (samples could be ASVs, OTUs,
  # species...); "turnover"/"nestedness" are already standard beta-diversity
  # terms on their own, so only "shared" needs the extra word.
  partition_label <- if (partition == "shared") "shared features" else partition
  y_lab <- if (!is.null(y_axis_title)) y_axis_title else paste0("Beta diversity (", partition_label, ")")

  if(!is.null(condition2_col)) {
    figura <- ggpubr::ggboxplot(
      beta_df, x="condition1_group", y="value", fill="condition1_group"
    ) +
      ggplot2::ylab(y_lab) +
      ggplot2::scale_fill_manual(values=group_colors) +
      ggplot2::labs(fill = "Comparison") +
      base_theme +
      ggh4x::facet_grid2(
        stats::as.formula(paste(". ~", paste0(condition2_col,".x"))),
        scales="free_x",
        strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill=facet_colors))
      ) +
      ggplot2::xlab(x_axis_title)
  } else {
    figura <- ggpubr::ggboxplot(
      beta_df, x="condition1_group", y="value", fill="condition1_group"
    ) +
      ggplot2::ylab(y_lab) +
      ggplot2::scale_fill_manual(values=group_colors) +
      ggplot2::labs(fill = "Comparison") +
      ggplot2::xlab(x_axis_title) +
      base_theme
  }

  if (!is.null(stat)) {
    figura <- figura + .mbm_stat_layer(stat, p_adjust_method)
  }

  return(figura)
}
