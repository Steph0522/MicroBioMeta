#' Generate a Lollipop plot from random forest results
#'
#' @param table Data frame where columns are samples and rows are ASVs or taxa.
#' @param metadata Data frame containing sample metadata.
#' @param variable_to_predict Variable to predict from the metadata.
#' @param top_n Number of top features to plot (default = 15).
#' @param col_palette Custom color palette (optional).
#' @param legend_figure Main title for the figure.
#' @param size the size of the point of the lollipop.
#' 
#' @return A lollipop plot showing top important features from random forest analysis.
#' @export
#'
#' @examples 
#' randomf_lollipop_plot(table = otu_table, 
#'                      metadata = sample_metadata,
#'                      variable_to_predict = "season",
#'                      legend_figure = "Top important features (Random Forest)")
randomf_lollipop_plot <- function(table,
                                  metadata,
                                  top_n = 15,
                                  size =8,
                                  variable_to_predict,
                                  col_palette = NULL,
                                  legend_figure = NULL) {
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  # Set default figure title if not provided
  if (is.null(legend_figure)) {
    legend_figure <- sprintf("Top %d most important features (Random Forest)", top_n)
  }
  
  # Extract taxonomy column for later use
  ncols <- ncol(table)
  table_numeric <- table[-ncols]
  taxonomy <- table[ncols]
  
  # Transpose if samples are in rows
  if (ncol(table_numeric) < nrow(table_numeric)) {
    table_numeric <- data.frame(t(table_numeric), check.names = F)
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
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Please install the randomForest package to use this function.")
  }
  
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
      Phylum = case_when(
        Phylum == "p__Proteobacteria" ~ "Pseudomonadota",
        Phylum == "p__Actinobacteriota" ~ "Actinomycetota",
        Phylum == "p__Bacteroidetes" ~ "Bacteroidota",
        Phylum == "p__Firmicutes" ~ "Bacillota",
        TRUE ~ sub("^p__", "", as.character(Phylum))  
      ),
      taxonomy2 = paste0(letters[1:n()], ".", taxonomy) # Add letter labels
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
  # Default color palette (colorblind-friendly)
  default_palette <- c(
    "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482",
    "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C",
    "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26"
  )
  
  # Use custom palette if provided, otherwise use default
  if (!is.null(col_palette)) {
    # Check if palette has enough colors
    n_phyla <- length(unique(top_asvs.modified$Phylum))
    if (length(col_palette) < n_phyla) {
      warning(sprintf("Custom palette has %d colors but %d are needed. Using default palette to complete.", 
                      length(col_palette), n_phyla))
      col_palette <- c(col_palette, default_palette)[1:n_phyla]
    }
    fill_colors <- col_palette
  } else {
    fill_colors <- default_palette
  }
  
  # Convert importance values to numeric
  top_asvs.modified$MeanDecreaseGini <- as.numeric(top_asvs.modified$MeanDecreaseGini)
  
  # Create lollipop plot
  lollipop <- ggplot2::ggplot(
    top_asvs.modified, 
    aes(x = reorder(taxonomy2, MeanDecreaseGini), 
        y = MeanDecreaseGini, 
        fill = Phylum)
  ) +
    ggplot2::geom_segment(
      aes(
        x = reorder(taxonomy2, MeanDecreaseGini), 
        xend = reorder(taxonomy2, MeanDecreaseGini), 
        y = 0, 
        yend = MeanDecreaseGini
      ),
      color = "black", 
      linewidth = 2
    ) +
    ggplot2::ylab("Feature importance (MeanDecreaseGini)") +
    ggplot2::geom_point(size = size, shape = 21, color = "black") +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::theme_classic() +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(legend_figure) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(size = 11, color = "black"),
      axis.text.y = ggplot2::element_text(size = 14, face = "italic", colour = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 20)
    )
  
  return(lollipop) 
}