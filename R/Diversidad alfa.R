##### FUNCIÓN PARA DIVERSIDAD ALFA ###


#1. Definir la funcion 

diversidad_alfa <- function(otu, metadata, x, y, fill, V1, V2, color1, color2, color3, color4, color5, Titulo.leyenda,
                            Titulo.figura, Titulo.eje.y)
  
  #Obtener indices de diversidad por cada orden
  #especificar el nombre de la columna con las muestras (OTUID o sample-id y verificar al unir con el metadata)
{otu.q0<- hilldiv::hill_div(otu,0) %>% as.data.frame() %>% mutate(orden="q=0") %>% rownames_to_column(var = "OTUID")
otu.q1<- hilldiv::hill_div(otu,1) %>% as.data.frame() %>% mutate(orden="q=1") %>% rownames_to_column(var = "OTUID")
otu.q2<- hilldiv::hill_div(otu,2) %>% as.data.frame() %>% mutate(orden="q=2") %>% rownames_to_column(var = "OTUID")

#renombrar nombre del metadata
colnames(metadata)[1]<- "OTUID"

#Unir archivos de cada orden en una sola tabla y darle nombre a la segunda columna del archivo.
otu.completa <- rbind(otu.q0, otu.q1, otu.q2) %>% inner_join(metadata, by = "OTUID")
colnames(otu.completa)[2]<- "Numero efectivo de ASV´s"

#Especificar el orden del eje x en el que quiero la grafica final.
#otu.completa$variable.ordenada<- factor(otu.completa$variable.a.ordenar, levels=niveles.ordenados)


#Figura completa
figura_completa<- otu.completa %>% 
  ggpubr::ggboxplot(x = x, y = y, fill = fill)+
  ggh4x::facet_grid2(as.formula(paste(V1,"~", V2)) , space = "fixed", scales = "free") +
  ylab(Titulo.eje.y)+
  xlab(NULL)+
  scale_fill_manual(values = c(color1, color2, color3, color4, color5))+ 
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
                color1="#d55e00",
                color2="#cc79a7",
                color3="#0072b2",
                color4="#f0e442",
                color5="#009e73",
                Titulo.leyenda="Days",
                Titulo.figura="Alfa diversity", 
                Titulo.eje.y= "Effective number of ASVs")


