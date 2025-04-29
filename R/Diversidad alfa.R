##### FUNCIÓN PARA DIVERSIDAD ALFA ###


#1. Definir la funcion 
#' Diversidad alfa
#'
#' @param table: Data frame, where, the columns are the samples and rows are ASV's or taxa.
#'               In addition, the first column should be called "OTUID", other names, e.g.,OTU-ID,
#'               Sample-ID, sample-ids, etc., will flag an error.
#' @param metadata Data frame of characteristics or important information of the samples.
#' @param x_col Variable that defines the axis x
#' @param y_col Variable that defines the axis y
#' @param fill_col Variable that separates the samples for color 
#' @param facet_x Variable that separates the samples in horizontal facets
#' @param facet_y Variable that separates the samples in vertical facets
#' @param col_pallete Blind-friendly color palette
#' @param legend_title Title of the principal legend 
#' @param legend_figure Principal title of the figure
#' @param axis_y_title Title of the axis y
#'
#' @return
#' @export
#'
#' @examples
diversidad_alfa <- function(table, metadata, x_col, y_col, fill_col, facet_x, facet_y, col_pallete, legend_title,
                            legend_figure, axis_y_title)
  
  
  #Obtener indices de diversidad por cada orden
  #cargar librerias
{library(tidyverse)
  library(hilldiv)
  library(ggplot2)
  library(ggpubr)
  #El nombre de la primera columna de la otu debe ser OTUID, el nombre del metadata se cambia en el siguiente paso. 
  otu.q0<- hilldiv::hill_div(otu,0) %>% as.data.frame() %>% mutate(orden="q=0") %>% rownames_to_column(var = "OTUID")
  otu.q1<- hilldiv::hill_div(otu,1) %>% as.data.frame() %>% mutate(orden="q=1") %>% rownames_to_column(var = "OTUID")
  otu.q2<- hilldiv::hill_div(otu,2) %>% as.data.frame() %>% mutate(orden="q=2") %>% rownames_to_column(var = "OTUID")
  
  #renombrar nombre del metadata
  colnames(metadata)[1]<- "OTUID"
  
  #Unir archivos de cada orden en una sola tabla y darle nombre a la segunda columna del archivo.
  otu.completa <- rbind(otu.q0, otu.q1, otu.q2) %>% inner_join(metadata, by = "OTUID")
  colnames(otu.completa)[2]<- "Numero efectivo de ASV´s"
  
  #Especificar el orden del eje x en el que quiero la grafica final.
  #otu.completa$variable.ordenada<- factor(otu.completa$variable.a.ordenar, levels=niveles.ordenados)
  
  
  paleta_colores <- c("#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482",
                      "#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C",
                      "#DCD300","#882D17","#8DB600","#654522","#E25822","#2B3D26")
  
  #Figura completa
  figura_completa<- otu.completa %>% 
    ggpubr::ggboxplot(x = x_col, y = y_col, fill = fill_col)+
    ggh4x::facet_grid2(as.formula(paste(V1,"~", V2)) , space = "fixed", scales = "free") +
    ylab(axis_y_title)+
    xlab(NULL)+
    scale_fill_manual(values = paleta_colores) + 
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(fill=legend_title) + #Modificar titulo de la leyenda
    ggtitle(legend_figure)+
    guides(color = guide_legend(override.aes = list(size = 5))) + #tamaño del key de la primer leyenda (season)
    theme(panel.border = element_rect(fill = "transparent", # Necesario para agregar el borde
                                      color = "black", linewidth = 0.5),
          strip.text.x = element_text(face = "bold", color = "white", size = 14),
          strip.text.y = element_text(face = "bold", color = "white", size = 14),
          strip.background.x = element_rect(fill = "#676778", linetype = "solid",
                                            color = "black", linewidth = 0.5),
          strip.background.y = element_rect(fill = "gray30", linetype = "solid",
                                            color = "black", linewidth = 0.5),
          axis.title.y = element_text(size=14, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11))
  
  return(figura_completa)
}

#2. Correr la funcion 
diversidad_alfa(table = table.qiime, 
                metadata = sample_metadata,
                x_col = "body_site",
                y_col="Numero efectivo de ASV´s",
                fill_col = "day",
                facet_x="subject",
                facet_y="orden", 
                col_pallete = paleta_colores,
                legend_title="Days",
                legend_figure="Alfa diversity", 
                axis_y_title= "Effective number of ASVs")


