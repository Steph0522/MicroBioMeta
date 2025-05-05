#####RANDOM FOREST CON LOLIPOP####

 #Definir la funcion 
#' Title
#'
#' @param table Data frame, where, the columns are the samples and rows are ASV's or taxa.
#' @param metadata Data frame of characteristics or important information of the samples.
#' @param taxonomy Data frame that contains the ASVs assigned to a taxonomic group.
#' @param variable_to_predict Variable to predict from the metadata of the samples analized.
#' @param col_pallete Blind-friendly color palette.
#' @param legend_figure Principal title of the figure.
#'
#' @return A lollipop plot with Top 15 most important ASV´s of random forest analysis .
#' @export
#'
<<<<<<< HEAD
#' @examples
random_forest(table = table, 
              metadata = metadata,
              taxonomy = taxonomia,
              variable_to_predict = "temporada",
              legend_figure =  "Top 15 most important ASVs (Random Forest)")
=======
#' @examples random_forest(table = gestacion.recto.filtrada, 
#'              metadata = metadata.gestacion.recto,
#'              taxonomy = taxonomia_0.99,
#'              variable_to_predict = "temporada",
#'              legend_figure =  "Top 15 most important ASVs (Random Forest)")
>>>>>>> 96efc1f3b7ffc43cc107ac8b590d67c23c5b8467
#' 
#' 
#' 
random_forest <- function(table, metadata, taxonomy, variable_to_predict, col_pallete= NULL, legend_figure)
  
  #cargar librerias
{library(tidyverse)
  library(ggplot2)
  library(randomForest)
  library(dplyr)
  
  # Asegurar que las muestras coincidan
  if (ncol(table) < nrow(table)) {
    table <- t(table)
  }
  
  # Verifica si los nombres de muestra coinciden
  common_samples <- intersect(rownames(table), rownames(metadata))
  
  # Filtrar ambos datasets
  otu_filtered <- table[common_samples, ]
  metadata_filtered <- metadata[common_samples, ]
  
  # Ejemplo: predecir la variable
  #variable_to_predict <- as.factor(metadata_filtered[[variable_to_predict]])
  if (!(variable_to_predict %in% colnames(metadata_filtered))) {
    stop("La variable a predecir no está en los metadatos.")
  }
  
  # Combinar variable temporada con datos OTU
  # datos_rf <- data.frame(variable_to_predict, otu_filtered)
  datos_rf <- data.frame(variable = metadata_filtered[[variable_to_predict]], otu_filtered)
  colnames(datos_rf)[1] <- variable_to_predict
  
  #modelo Random Forest
  modelo_rf <- randomForest(as.formula(paste(variable_to_predict, "~ .")), data = datos_rf, importance = TRUE, ntree = 500)
  
  # Obtener la importancia de las variables
  importancia <- importance(modelo_rf)
  
  # Convertir a data frame con nombres de OTUs
  importancia_df <- data.frame(OTU = rownames(importancia), importancia)
  
  # Asegúrate de que la taxonomía tenga un campo OTU si rownames no están como columna
  taxonomy$OTU <- rownames(taxonomy)
  
  # Hacemos merge para unir importancia + taxonomía
  importancia_taxa <- merge(importancia_df, taxonomy, by = "OTU")
  
  # Ordenar de mayor a menor importancia
  importancia_taxa_ordenada <- importancia_taxa %>% arrange(desc(MeanDecreaseGini))
  
  # Mostrar top 10
  #head(importancia_taxa_ordenada, 15)
  
  # Crear una columna con taxonomía simplificada, por ejemplo usando Género
  importancia_taxa_ordenada$Taxon <- ifelse(
    is.na(importancia_taxa_ordenada$Genus) | importancia_taxa_ordenada$Genus == "",
    importancia_taxa_ordenada$Family,  # usa Familia si no hay Género
    importancia_taxa_ordenada$Genus)
  
  #Seleccionar las 15 más importantes
  top15 <- importancia_taxa_ordenada %>% 
    slice_max(order_by = MeanDecreaseGini, n = 15)
  
  #tabla
  top15_modificada <- top15 %>%
    mutate(Taxon = gsub("^g__|^f__|;", "", Taxon)) %>%
    #mutate_at(c("Taxon"), funs(Taxon=case_when(Taxon=="g__Incertae_Sedis;" ~ "Ruminococcaceae", TRUE~as.character(Taxon)))) %>%
    mutate_at(c("Phylum"), funs(Phylum=case_when(Phylum=="p__Proteobacteria;" ~ "Pseudomonadata",
                                                 Phylum=="p__Firmicutes;" ~ "Bacillota",  TRUE~as.character(Phylum)))) %>%
    mutate(Taxon2=paste0(LETTERS[1:n()], ".", Taxon)) #agregar letras al inicio del nombre en la columna nueva llamda Taxon2
  
  
  # Paleta por defecto si no se proporciona
  if (is.null(col_pallete)) {
    col_pallete <- c("#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482",
                     "#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C",
                     "#DCD300","#882D17","#8DB600","#654522","#E25822","#2B3D26")
  }
  
  
  #Plot
  lollipop_randomR <- ggplot(top15_modificada, aes(x = reorder(Taxon2, + MeanDecreaseGini), y = MeanDecreaseGini, fill = Phylum)) +
    geom_segment(aes(x = reorder(Taxon2, + MeanDecreaseGini), 
                     xend = reorder(Taxon2, + MeanDecreaseGini), 
                     y = 0, yend = MeanDecreaseGini), 
                 color = "grey30", lwd = 2) +
    ylab("Feature importance") +
    geom_point(size = 9, pch = 21, col = "grey30") +
    scale_fill_manual(values = col_pallete) + 
    theme_classic() +
    coord_flip() +
    ggtitle(legend_figure) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 13, face = "bold"),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 14, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 22, face = "bold"),
          legend.background = element_rect(colour="black"))
  
  
  return(lollipop_randomR) 
}


