#' Generate a Lollipop plot from random forest results
#'
#' @param table Data frame where columns are samples and rows are ASVs or taxa.
#' @param metadata Data frame containing sample metadata.
#' @param variable_to_predict Variable to predict from the metadata.
#' @param top_n Number of top features to plot (default = 15).
#' @param group_colors Custom color palette (optional).
#' @param title Main title for the figure.
#' @param size the size of the point of the lollipop.
#' @param save_table Logical. If \code{TRUE}, saves the feature-importance
#'   table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"randomforest_importance.txt"}.
#'
#' @return A lollipop plot showing top important features from random forest analysis.
#' @importFrom randomForest randomForest importance
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
#' random_forest_lollipop_plot(
#'   table               = table,
#'   metadata            = metadata,
#'   variable_to_predict = "Location",
#'   top_n               = 20,
#'   size                = 6
#' )
#' }
random_forest_lollipop_plot <- function(table,
                                  metadata,
                                  top_n = 15,
                                  size =8,
                                  variable_to_predict,
                                  group_colors = NULL,
                                  title = NULL,
                                  save_table = FALSE,
                                  table_filename = "randomforest_importance.txt") {

  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")

  # Set default figure title if not provided
  if (is.null(title)) {
    title <- sprintf("Top %d most important features (Random Forest)", top_n)
  }
  
  # Extract taxonomy column for later use
  ncols <- ncol(table)
  table_numeric <- table[-ncols]
  taxonomy <- table[ncols]
  
  # Transpose if samples are in rows
  if (ncol(table_numeric) < nrow(table_numeric)) {
    table_numeric <- data.frame(t(table_numeric), check.names = FALSE)
  }
  
  # Find common samples between table and metadata
  common_samples <- intersect(rownames(table_numeric), metadata[[1]])
  
  # Filter both datasets to include only common samples
  otu_filtered <- table_numeric[common_samples, , drop = FALSE]
  metadata_filtered <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
  
  # Convert response variable to factor
  response <- as.factor(metadata_filtered[[variable_to_predict]])
  
  # Validate response variable has at least 2 classes
  if (length(unique(response)) < 2) {
    stop("The variable to predict must contain at least two classes.")
  }
  
  # Run Random Forest analysis
  rf_model <- randomForest::randomForest(
    x = otu_filtered,
    y = response,
    importance = TRUE,
    ntree = 500
  )
  
  # Extract feature importance
  importance_df <- randomForest::importance(rf_model)
  importance_df <- as.data.frame(importance_df)
  
  # Get top features by MeanDecreaseGini
  importance_df$ASV <- rownames(importance_df)
  top_asvs <- importance_df %>% 
    dplyr::arrange(desc(MeanDecreaseGini)) %>% 
    head(top_n)
  
  # Merge with original taxonomy
  top_asvs <- dplyr::left_join(
    top_asvs,
    data.frame(ASV = colnames(as.data.frame(table_numeric)), taxonomy = taxonomy),
    by = "ASV"
  ) %>%
    dplyr::mutate(taxonomy_original = taxonomy)
  
  # Parse taxonomy into different levels
  top_asvs <- top_asvs %>%
    tidyr::separate(
      col = taxonomy_original,
      into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = ";",
      fill = "right" 
    ) %>%
    dplyr::mutate(across(everything(), ~ trimws(.)))
  
  # Simplify taxonomy labels for plotting
  top_asvs.modified <- top_asvs %>%
    dplyr::mutate(
      taxonomy = dplyr::case_when(
        grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ 
          sub(".*g__([^;]*).*", "\\1", taxonomy),
        grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ 
          paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
        grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ 
          paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
        grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ 
          paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
        grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ 
          paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
        TRUE ~ "Unclassified"
      ),
      Phylum = dplyr::case_when(
        Phylum == "p__Proteobacteria" ~ "Pseudomonadota",
        Phylum == "p__Actinobacteriota" ~ "Actinomycetota",
        Phylum == "p__Bacteroidetes" ~ "Bacteroidota",
        Phylum == "p__Firmicutes" ~ "Bacillota",
        TRUE ~ sub("^p__", "", as.character(Phylum))
      ),
      # A unique row id (not the displayed label) keeps each bar distinct on
      # the axis even when two features share the same simplified taxonomy
      # name; the axis is then labeled with the plain taxonomy text via
      # scale_x_discrete(labels = ...) below.
      row_id = paste0("row", dplyr::row_number())
    )
  # Generate warning if any updated phyla are present
  updated_phyla <- c("p__Proteobacteria", "p__Actinobacteriota", 
                     "p__Bacteroidetes", "p__Firmicutes")
  present_phyla <- intersect(updated_phyla, top_asvs$Phylum)
  
  if(length(present_phyla) > 0) {
    old_new <- data.frame(
      Old = c("Proteobacteria", "Actinobacteriota", "Bacteroidetes", "Firmicutes"),
      New = c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Bacillota")
    )
    
    changes <- old_new[old_new$Old %in% sub("^p__", "", present_phyla), ]
    
    warning("Note: Some bacterial phylum names have been updated to match NCBI's revised taxonomy:\n",
            paste(sprintf("  - '%s' changed to '%s'", changes$Old, changes$New), collapse = "\n"),
            "\nReference: https://ncbiinsights.ncbi.nlm.nih.gov/2021/12/10/ncbi-taxonomy-prokaryote-phyla-added/")
  }
  # Use custom palette if provided, otherwise use Okabe-Ito default
  n_phyla <- length(unique(top_asvs.modified$Phylum))
  if (!is.null(group_colors)) {
    if (length(group_colors) < n_phyla) {
      warning(sprintf("Custom palette has %d colors but %d are needed. Recycling palette.",
                      length(group_colors), n_phyla))
      group_colors <- rep_len(group_colors, n_phyla)
    }
    fill_colors <- group_colors
  } else {
    fill_colors <- rep_len(.mbm_colors, n_phyla)
  }
  
  # Convert importance values to numeric
  top_asvs.modified$MeanDecreaseGini <- as.numeric(top_asvs.modified$MeanDecreaseGini)

  if (save_table) {
    utils::write.table(top_asvs.modified, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # Create lollipop plot
  row_labels <- setNames(top_asvs.modified$taxonomy, top_asvs.modified$row_id)
  lollipop <- ggplot2::ggplot(
    top_asvs.modified,
    ggplot2::aes(x = reorder(row_id, MeanDecreaseGini),
        y = MeanDecreaseGini,
        fill = Phylum)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = reorder(row_id, MeanDecreaseGini),
        xend = reorder(row_id, MeanDecreaseGini),
        y = 0,
        yend = MeanDecreaseGini
      ),
      color = "black",
      linewidth = 1
    ) +
    ggplot2::ylab("Feature importance (MeanDecreaseGini)") +
    ggplot2::geom_point(size = size, shape = 21, color = "black") +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::scale_x_discrete(labels = row_labels) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title) +
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        panel.grid   = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y  = ggplot2::element_text(size = 12, color = "black",
                                              face = "italic")
      )
    )
  
  return(lollipop) 
}
