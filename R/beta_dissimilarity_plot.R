#' Beta Diversity Boxplot
#'
#' This function calculates beta diversity (shared species, turnover, or nestedness)
#' and generates a ggplot2 boxplot with facets and custom coloring. The user only
#' needs to specify the metadata column(s) to be used for comparisons; the function
#' constructs the comparison pairs internally and removes duplicates (A_vs_B = B_vs_A).
#'
#' @param table Abundance matrix (samples in columns, species/features in rows).
#' @param metadata Data frame with sample metadata. First column must match sample names in table.
#' @param comparison_condition1 Optional vector of comparison labels for the first condition.
#' @param condition1_col Column name in metadata for the first condition.
#' @param condition2_col Optional column name in metadata for the second condition (used as facet).
#' @param facet_colors Optional vector of colors for facet strips. Defaults to a neutral \code{"grey85"} background.
#' @param group_colors Optional named vector of colors for x-axis groups. Defaults to the package's
#'   colorblind-friendly Okabe-Ito palette (\code{.mbm_colors}, orange/blue first).
#' @param x_axis_title Title for the x-axis.
#' @param partition Type of beta diversity to compute: "shared", "turnover", or "nestedness".
#' @param family Family for turnover/nestedness calculation: "sorensen" or "jaccard".
#' @param stat Character or \code{NULL}. Statistical test to compare groups,
#'   passed to \code{ggpubr::stat_compare_means()} (e.g. \code{"wilcox.test"},
#'   \code{"kruskal.test"}, \code{"anova"}). Default \code{NULL} (no test shown).
#' @param x_label_angle Numeric. Rotation (in degrees) of the x-axis tick
#'   labels. Default \code{0} (horizontal); use e.g. \code{45} or \code{90} when
#'   comparison names are long enough to overlap.
#' @param strip_text_bold Logical. If \code{TRUE}, facet strip labels are bold. Default \code{FALSE}.
#' @param aspect_ratio Numeric. Sets the aspect ratio (height/width) of each panel. Default \code{NULL} (automatic).
#' @param save_table Logical. If \code{TRUE}, saves the beta diversity table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table. Default \code{"betadiv_table.txt"}.
#'
#' @return A ggplot2 figure object.
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
#' }

beta_dissimilarity_plot <- function(
    table, metadata,
    comparison_condition1 = NULL,
    condition1_col,
    condition2_col = NULL,
    facet_colors = NULL,
    group_colors = NULL,
    x_axis_title = "Condition",
    partition = c("shared","turnover","nestedness"),
    family = c("sorensen","jaccard"),
    stat = NULL,
    x_label_angle = 0,
    strip_text_bold = FALSE,
    aspect_ratio = NULL,
    save_table = FALSE,
    table_filename = "betadiv_table.txt"
) {
  
  # --- Remove taxonomy column if present ---
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) == 1) table <- table[ , -tax_col]
  
  partition <- match.arg(partition)
  family <- match.arg(family)
  
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
  
  # --- Convert to long format ---
  beta_df <- reshape2::melt(as.matrix(beta_mat), varnames=c("site1","site2"))
  beta_df <- dplyr::filter(beta_df, !is.na(value) & value != 0)
  
  beta_df <- dplyr::left_join(beta_df, metadata, by=c("site1"=names(metadata)[1])) %>%
    dplyr::left_join(metadata, by=c("site2"=names(metadata)[1]))
  
  # --- Construct comparison labels ---
  beta_df <- dplyr::mutate(beta_df,
                           condition1_group = paste0(pmin(.data[[paste0(condition1_col,".x")]], .data[[paste0(condition1_col,".y")]]),
                                                     "_vs_",
                                                     pmax(.data[[paste0(condition1_col,".x")]], .data[[paste0(condition1_col,".y")]])))
  beta_df <- dplyr::distinct(beta_df, site1, site2, .keep_all = TRUE) # remove duplicates
  
  # Filter if comparison_condition1 specified
  if(!is.null(comparison_condition1)) {
    beta_df <- dplyr::filter(beta_df, condition1_group %in% comparison_condition1)
  }
  
  
  if (save_table) {
    utils::write.table(
      beta_df,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    message(paste("Table saved as:", table_filename))
  }
  # --- Default palettes if missing ---
  n_groups <- length(unique(beta_df$condition1_group))
  if(is.null(group_colors)) {
    group_colors <- if (n_groups == 2) .mbm_colors_2group else rep_len(.mbm_colors, n_groups)
  }
  if(!is.null(names(group_colors)==FALSE)) names(group_colors) <- unique(beta_df$condition1_group)

  if(!is.null(condition2_col)) {
    n_facets <- length(unique(beta_df[[paste0(condition2_col,".x")]]))
    if(is.null(facet_colors)) facet_colors <- rep("grey85", n_facets)
  }
  

  base_theme <- .mbm_theme(
    legend_position = "right",
    extra = ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      strip.text  = .mbm_strip_text(strip_text_bold),
      axis.text.x = .mbm_x_text(x_label_angle)
    )
  )
  if (!is.null(aspect_ratio)) {
    base_theme <- base_theme + ggplot2::theme(aspect.ratio = aspect_ratio)
  }

  # --- Create figure ---
  if(!is.null(condition2_col)) {
    figura <- ggpubr::ggboxplot(
      beta_df, x="condition1_group", y="value", fill="condition1_group"
    ) +
      ggplot2::ylab(paste0("Beta diversity (",partition,")")) +
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
      ggplot2::ylab(paste0("Beta diversity (",partition,")")) +
      ggplot2::scale_fill_manual(values=group_colors) +
      ggplot2::xlab(x_axis_title) +
      base_theme
  }

  if (!is.null(stat)) {
    figura <- figura + ggpubr::stat_compare_means(
      method = stat,
      mapping = ggplot2::aes(
        label = paste0("p = ", scales::label_pvalue(accuracy = 0.001)(ggplot2::after_stat(p)))
      ),
      size = 3.5, family = "serif", hide.ns = TRUE
    )
  }

  return(figura)
}
