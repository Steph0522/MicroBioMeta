#####RANDOM FOREST WITH LOLIPOP####
#'
#' @param table Data frame, where, the columns are the samples and rows are ASV's or taxa.
#' @param metadata Data frame of characteristics or important information of the samples.
#' @param variable_to_predict Variable to predict from the metadata of the samples analized.
#' @param col_pallete Blind-friendly color palette.
#' @param legend_figure Principal title of the figure.
#'
#' @return A lollipop plot with Top 15 most important ASV´s of random forest analysis .
#' @export
#'
#' @examples random_forest(table = table, 
#'              metadata = metadata,
#'              variable_to_predict = condition (e.g."season", "environment", "soil")
#'              legend_figure =  "Top 15 most important ASVs (Random Forest)")
#' 
#1. Define function 
random_forest <- function(table, metadata, variable_to_predict, col_pallete= NULL, legend_figure)
  
{#Eliminate taxonomy column for numeric analysis
  taxonomy <- table$taxonomy
  table_numeric <- table %>% dplyr::select(-taxonomy) 
  
  if (ncol(table_numeric) < nrow(table_numeric)) {
    table_numeric <- t(table_numeric)
  }
  
  #Verify that samples names samples match
  common_samples <- intersect(rownames(table_numeric), rownames(metadata))
  
  # Filter both datasets
  otu_filtered <- table_numeric[common_samples, ]
  metadata_filtered <- metadata[common_samples, ]
  
  # Define response variable 
  response <- as.factor(metadata_filtered[[variable_to_predict]])
  
  # Model Random Forest
  modelo_rf <- randomForest::randomForest(x = otu_filtered, y = response, importance = TRUE, ntree = 500)
  
  # Obtaine importance
  importance_df <- randomForest::importance(modelo_rf)
  importance_df <- as.data.frame(importance_df)
  
  # Order MeanDecreaseGini and get the 15 most important ASVs
  importance_df$ASV <- rownames(importance_df)
  top_asvs <- importance_df %>% dplyr::arrange(desc(MeanDecreaseGini)) %>% head(15)
  
  #Duplicate taxonomy column for modify file
  top_asvs <- dplyr::left_join(top_asvs, data.frame(ASV = colnames(table_numeric), taxonomy = taxonomy), by = "ASV") %>%
    dplyr::mutate(taxonomy_original = taxonomy)
  
  #Separate taxonomy column for agregate Phylum legend
  top_asvs <- top_asvs %>%
    tidyr::separate(col= taxonomy_original, into = c("Dominio","Phylum","Class","Orden","Family","Genus","Specie"), sep = ";") 
  
  #Remove all white space
  top_asvs<- top_asvs %>%
    dplyr::mutate(across(everything(), ~ trimws(.)))
  
  #Modify table
  top_asvs.modificada <- top_asvs %>%
    dplyr::mutate(taxonomy = dplyr::case_when(
      grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
      grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
      grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
      grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
      grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
      TRUE ~ "Unclassified")) %>%
    dplyr::mutate(Phylum = case_when(Phylum=="p__Proteobacteria" ~ "Pseudomonadata",
                                     Phylum=="p__Firmicutes" ~ "Bacillota",  TRUE~as.character(Phylum))) %>%
    mutate(taxonomy2=paste0(LETTERS[1:n()], ".", taxonomy)) #add letters
  
  
  paleta_colores <- c("#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482",
                      "#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C",
                      "#DCD300","#882D17","#8DB600","#654522","#E25822","#2B3D26")
  
  #Convert MeanDecreaseGini column to numeric
  top_asvs.modificada$MeanDecreaseGini <- as.numeric(top_asvs.modificada$MeanDecreaseGini)
  
  #FigurE
  lollipop <- ggplot2::ggplot(top_asvs.modificada, aes(x = reorder(taxonomy2, MeanDecreaseGini), y = MeanDecreaseGini, fill = Phylum)) +
    ggplot2::geom_segment(aes(x = reorder(taxonomy2, MeanDecreaseGini), 
                              xend = reorder(taxonomy2, MeanDecreaseGini), 
                              y = 0, yend = MeanDecreaseGini), 
                          color = "grey30", lwd = 2) +
    ggplot2::ylab("Feature importance") +
    ggplot2::geom_point(size = 9, pch = 21, col = "grey30") +
    ggplot2::scale_fill_manual(values = paleta_colores) + 
    ggplot2::theme_classic() +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(legend_figure) +
    ggplot2::theme(axis.title.y = element_blank(),
                   axis.title.x = element_text(size = 13, face = "bold"),
                   axis.text.x = element_text(size = 11),
                   axis.text.y = element_text(size = 14, face = "bold"),
                   legend.position = "bottom",
                   legend.title = element_text(size = 11, face = "bold"),
                   legend.text = element_text(size = 10),
                   plot.title = element_text(size = 20))
  
  return(lollipop) 
}

