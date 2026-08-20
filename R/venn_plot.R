#' Generate a Venn diagram of taxa shared between sample groups
#'
#' This function creates a Venn diagram using either the `ggVennDiagram` or `ggenn` package
#' based on the user's preference. It includes additional customization options like
#' minimum prevalence filtering and custom group color scales (manual).
#'
#' @param table A data frame containing taxonomic abundance data with a column named `taxonomy` and subsequent columns as sample IDs.
#' @param metadata A data frame containing metadata with a column named `SAMPLEID` that matches the sample columns in `table`.
#' @param merge_by A character string specifying the metadata column by which to group and merge samples.
#' @param selected_samples Optional character vector specifying a subset of sample IDs to include in the analysis.
#' @param min_prevalence Optional numeric value (0-1) to filter taxa based on minimum prevalence across groups.
#' @param title Optional character string for the title of the plot.
#' @param method Character string: `ggvenn` (default) or `ggVennDiagram`, specifying the package to use for Venn diagram generation.
#' @param group_colors Optional vector of colors for the groups. If NULL, a default `distiller` scale with `Set3` palette will be used.
#' @param save_table Logical. If \code{TRUE}, saves a long-format table of taxa
#'   membership per group to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"venn_taxa_sets.txt"}.

#' @return A ggplot object or other plot depending on the method.
#' @importFrom ggvenn ggvenn
#' @importFrom ggVennDiagram ggVennDiagram
#' @export
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SAMPLEID"
#'
#' venn_plot(
#'   table          = table,
#'   metadata       = metadata,
#'   merge_by       = "Location",
#'   min_prevalence = 0
#' )
#'
#' ## Filtering taxa by prevalence
#' venn_plot(
#'   table          = table,
#'   metadata       = metadata,
#'   merge_by       = "Location",
#'   min_prevalence = 0.2
#' )
#'
#' ## Custom colors
#' venn_plot(
#'   table          = table,
#'   metadata       = metadata,
#'   merge_by       = "Location",
#'   min_prevalence = 0,
#'   group_colors   = c("#1B9E77", "#D95F02")
#' )
#' }
venn_plot <- function(table, metadata, merge_by = NULL,
                              selected_samples = NULL, min_prevalence = 0,
                              title = NULL, method = "ggvenn",
                              group_colors = NULL,
                              save_table = FALSE,
                              table_filename = "venn_taxa_sets.txt") {
  table <- as.data.frame(table)
  metadata <- as.data.frame(metadata)
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  if (!merge_by %in% colnames(metadata)) stop("Group or merge column is not in the metadata file.")
  
  common_samples <- intersect(colnames(table), metadata[[1]])
  table <- table[, c(common_samples), drop = FALSE]
  metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
  rownames(metadata) <- NULL
  
  if (!is.null(selected_samples)) {
    common_samples <- intersect(selected_samples, colnames(table))
    table <- table[, c(common_samples), drop = FALSE]
    metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
  }
  
  metadata_split <- split(metadata[[1]], metadata[[merge_by]])
  metadata_split <- metadata_split[sapply(metadata_split, length) > 0]
  num_groups <- length(metadata_split)
  
  # Default palette - Okabe-Ito (colorblind-friendly)
  use_manual <- !is.null(group_colors)
  if (!use_manual) {
    group_colors <- if (num_groups == 2) {
      .mbm_colors_2group
    } else {
      rep_len(.mbm_colors, num_groups)
    }
  } else {
    group_colors <- rep(group_colors, length.out = num_groups)
  }
  
  # Build the list of sets to compare
  lista <- lapply(metadata_split, function(samps) {
    subset <- table[, samps, drop = FALSE]
    subset_core <- if (min_prevalence > 0) {
      subset_core <- subset[rowMeans(subset > 0) >= min_prevalence, ]
      
    } else {
      subset_core <- subset[rowSums(subset) != 0, ]

    }
    rownames(subset_core)
  })
  names(lista) <- names(metadata_split)

  if (save_table) {
    venn_table <- utils::stack(lista)
    colnames(venn_table) <- c("taxon_id", "group")
    utils::write.table(venn_table, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # Method selection
  if (tolower(method) == "ggvenndiagram") {
    venn_plot <- ggVennDiagram::ggVennDiagram(
      lista,
      label_alpha = 0,
      set_color   = group_colors,
      edge_size   = 1
    ) +
      ggplot2::scale_fill_gradient(low = "white", high = "grey60",
                                   na.value = NA, name = "Count")

  } else if (tolower(method) == "ggvenn") {
    venn_plot <- ggvenn::ggvenn(lista,
                                fill_color   = group_colors,
                                fill_alpha   = 0.5,
                                stroke_color = "grey30",
                                set_name_size = 5,
                                text_size     = 4)
  } else {
    stop("Method must be 'ggvenn' or 'ggVennDiagram'.")
  }

  # Base theme for Venn (theme_void keeps circles clean)
  venn_plot <- venn_plot +
    ggplot2::theme_void(base_family = "serif") +
    ggplot2::theme(
      legend.position = "right",
      legend.text     = ggplot2::element_text(size = 12, color = "black"),
      legend.title    = ggplot2::element_text(size = 14, face = "bold",
                                              color = "black")
    )

  # Add title only when provided
  if (!is.null(title)) {
    venn_plot <- venn_plot +
      ggplot2::labs(title = title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                           size = 14, color = "black")
      )
  }
  return(venn_plot)
}

