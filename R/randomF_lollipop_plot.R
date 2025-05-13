#' Generate a Lollipot of a random forest.
#'
#' @param table Data frame, where, the columns are the samples and rows are ASV's or taxa.
#' @param metadata Data frame of characteristics or important information of the samples.
#' @param variable_to_predict Variable to predict from the metadata of the samples analized.
#' @param top_n Number of features to plot.
#' @param col_pallete Blind-friendly color palette.
#' @param legend_figure Principal title of the figure.
#' 
#'
#' @return A lollipop plot with Top 15 most important ASV´s of random forest analysis.
#' @export
#'
#' @examples randomf_lollipop_plot(table = table, 
#'              metadata = metadata,
#'              variable_to_predict = "season",
#'              legend_figure = "Top 15 most important ASVs (Random Forest)")
#' 

randomf_lollipop_plot <- function(table,
                                  metadata,
                                  top_n = 15,
                                  variable_to_predict,
                                  col_pallete = NULL,
                                  legend_figure = NULL) {
  
  # Asignar valor por defecto a `legend_figure`
  if (is.null(legend_figure)) {
    legend_figure <- sprintf("Top %d most important ASVs (Random Forest)", top_n)
  }
  
  # Eliminar columna de taxonomy para análisis numérico
  taxonomy <- table$taxonomy
  table_numeric <- table %>% dplyr::select(-taxonomy)
  
  # Transponer si es necesario
  if (ncol(table_numeric) < nrow(table_numeric)) {
    table_numeric <- t(table_numeric)
  }
  
  # Verificar que los nombres de las muestras coincidan
  common_samples <- intersect(rownames(table_numeric), metadata[[1]])
  
  # Filtrar ambos datasets
  otu_filtered <- table_numeric[common_samples, , drop = FALSE]
  metadata_filtered <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
  
  # Definir variable respuesta
  response <- as.factor(metadata_filtered[[variable_to_predict]])
  
  # Verificar que haya al menos dos clases en la variable a predecir
  if (length(unique(response)) < 2) {
    stop("La variable a predecir debe contener al menos dos clases.")
  }
  
  # Modelado Random Forest
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("Install randomForest package.")
  modelo_rf <- randomForest::randomForest(
    x = otu_filtered,
    y = response,
    importance = TRUE,
    ntree = 500
  )
  
  # Obtener importancia
  importance_df <- randomForest::importance(modelo_rf)
  importance_df <- as.data.frame(importance_df)
  
  # Ordenar MeanDecreaseGini y obtener los 15 ASVs más importantes
  importance_df$ASV <- rownames(importance_df)
  top_asvs <- importance_df %>% dplyr::arrange(desc(MeanDecreaseGini)) %>% head(top_n)
  
  # Unir con taxonomía original
  top_asvs <- dplyr::left_join(
    top_asvs,
    data.frame(ASV = colnames(table_numeric), taxonomy = taxonomy),
    by = "ASV"
  ) %>%
    dplyr::mutate(taxonomy_original = taxonomy)
  
  # Separar taxonomía y crear leyenda por Phylum
  top_asvs <- top_asvs %>%
    tidyr::separate(
      col = taxonomy_original,
      into = c("Dominio", "Phylum", "Class", "Orden", "Family", "Genus", "Specie"),
      sep = ";"
    ) %>%
    dplyr::mutate(across(everything(), ~ trimws(.)))
  
  # Modificar la tabla para crear etiquetas simplificadas
  top_asvs.modificada <- top_asvs %>%
    dplyr::mutate(
      taxonomy = dplyr::case_when(
        grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
        grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
        grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
        grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
        grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
        TRUE ~ "Unclassified"
      ),
      Phylum = case_when(
        Phylum == "p__Proteobacteria" ~ "Pseudomonadata",
        Phylum == "p__Firmicutes" ~ "Bacillota",
        TRUE ~ as.character(Phylum)
      ),
      taxonomy2 = paste0(letters[1:n()], ".", taxonomy) # Agregar letras como etiquetas
    )
  
  # Paleta de colores
  paleta_colores <- c(
    "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482",
    "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C",
    "#DCD300", "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26"
  )
  
  # Convertir MeanDecreaseGini a numérico
  top_asvs.modificada$MeanDecreaseGini <- as.numeric(top_asvs.modificada$MeanDecreaseGini)
  
  # Crear gráfico Lollipop
  lollipop <- ggplot2::ggplot(top_asvs.modificada, aes(x = reorder(taxonomy2, MeanDecreaseGini), y = MeanDecreaseGini, fill = Phylum)) +
    ggplot2::geom_segment(
      aes(
        x = reorder(taxonomy2, MeanDecreaseGini), 
        xend = reorder(taxonomy2, MeanDecreaseGini), 
        y = 0, yend = MeanDecreaseGini
      ),
      color = "black", lwd = 2
    ) +
    ggplot2::ylab("Feature importance") +
    ggplot2::geom_point(size = 9, pch = 21, col = "black") +
    ggplot2::scale_fill_manual(values = paleta_colores) + 
    ggplot2::theme_classic() +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(legend_figure) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 13, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 11,  color = "black"),
      axis.text.y = ggplot2::element_text(size = 14, face = "bold.italic", colour = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 20)
    )
  
  return(lollipop) 
}
