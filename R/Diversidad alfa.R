##### FUNCIÓN PARA DIVERSIDAD ALFA ###


#1. Definir la funcion 
#' Diversidad alfa
#'
#' @param otu: Otu table, en donde las columnas son las muestras y las filas ASV's. Ademas, la primer columna debe decir "OTUID", otros nombres marcara error (ejemplo OTU-ID, Sample-ID).
#' @param metadata Mapa de características que contiene información util de las muestras.
#' @param x Nombre de la variable del eje x
#' @param y Nombre de la variable del eje y
#' @param fill Nombre de la variable que separará las muestras por color 
#' @param V1 Variable que dividirá la figura en facets verticales
#' @param V2 Varible que separara la figura en facets horizontales
#' @param paleta_colores Paleta de colores blind-fliendly 
#' @param Titulo.leyenda Título de la leyenda principal 
#' @param Titulo.figura Título principal de la figura
#' @param Titulo.eje.y Título del eje y
#'
#' @return
#' @export
#'
#' @examples
diversidad_alfa <- function(otu, metadata, x, y, fill, V1, V2, paleta_colores, Titulo.leyenda,
                            Titulo.figura, Titulo.eje.y)
  
  
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
    ggpubr::ggboxplot(x = x, y = y, fill = fill)+
    ggh4x::facet_grid2(as.formula(paste(V1,"~", V2)) , space = "fixed", scales = "free") +
    ylab(Titulo.eje.y)+
    xlab(NULL)+
    scale_fill_manual(values = paleta_colores) + 
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(fill=Titulo.leyenda) + #Modificar titulo de la leyenda
    ggtitle(Titulo.figura)+
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
diversidad_alfa(otu = table.qiime, 
                metadata = sample_metadata,
                x = "body_site",
                y="Numero efectivo de ASV´s",
                fill = "day",
                V1="orden",
                V2="subject", 
                paleta_colores = paleta_colores,
                Titulo.leyenda="Days",
                Titulo.figura="Alfa diversity", 
                Titulo.eje.y= "Effective number of ASVs")


