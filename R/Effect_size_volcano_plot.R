#' Effect size volcano plot
#'
#' @param tabla Tabla de conteos donde las columnas son las muestras y las filas los 'features'
#' @param metadata Tabla de metadatos con información de las muestras
#' @param columna_condiciones Columna dentro de la tabla de metadatos que 
#' se tomará como vector de condiciones para hacer la transofmración con aldex.clr
#' @param color_menor Color para los puntos por debajo del valor de lim_inf
#' @param color_mayor Color para los puntos por encima del valor de lim_sup
#' @param lim_inf Valor negativo donde cruza el eje x 
#' @param lim_sup Valor positivo donde cruza el eje x
#' @param condicion1 Nombre de la condición que se ubica primero en la tabla de conteos
#' @param mostrar_etiquetas Opción de incluir o no la etiqueta de mayor y menor
#'
#' @return
#' @export
#'
#' @examples
#
#' 
#' 
#' 
effect_size_plot <-
  function(tabla,
           metadata,
           columna_condiciones,
           color_menor,
           color_mayor,
           threshold_lower,
           threshold_upper,
           condicion1,
           mostrar_etiquetas = TRUE) {  
    library(ALDEx2)
    library(ggplot2)
    library(scales)
    library(ggtext)
    
    condiciones <- metadata[[columna_condiciones]]
    
    aldex_clr <-
      aldex.clr(tabla, condiciones, mc.samples = 128, denom = "all")
    
    effect_size <- aldex.effect(
      aldex_clr,
      verbose = TRUE,
      include.sample.summary = FALSE,
      useMC = FALSE,
      CI = FALSE
    )
    
    KW <- aldex.kw(aldex_clr, useMC = FALSE, verbose = FALSE)
    
    resultado <- cbind(effect_size, KW)
    
    # Clasificar en 3 grupos autom?ticamente
    resultado$grupo <- ifelse(
      resultado$effect <= threshold_lower,
      "Menor en condicion 1",
      ifelse(
        resultado$effect >= threshold_upper,
        "Mayor en condicion 1",
        "Normal"
      )
    )
    
    p <- ggplot(resultado, aes(x = effect, y = kw.ep, color = grupo)) +
      geom_point(size = 3.5) +
      scale_color_manual(
        values = c(
          "Menor en condicion 1" = color_menor,
          "Normal" = "gray",
          "Mayor en condicion 1" = color_mayor
        )
      ) +
      geom_vline(
        xintercept = c(threshold_lower, threshold_upper),
        linetype = 2,
        color = "black"
      ) +
      geom_hline(yintercept = 0.05, linetype = 2, color = "black") +
      labs(
        x = "Effect size",
        y = bquote(italic("p") ~ " value"),
        color = "Grupo"  # T?tulo de la leyenda
      ) +
      ylim(c(0, 1)) +
      theme_classic() +
      theme(
        axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12)
      ) +
      theme(legend.position = "none") 
    
    # Agregar etiquetas solo si mostrar_etiquetas = TRUE
    if (mostrar_etiquetas) {
      p <- p +
        geom_richtext(
          aes(x = threshold_lower, y = 0.95, 
              label = paste("<b>Menor en", condicion1, "</b>")),
          color = color_menor,
          size = 5,
          hjust = 1,
          vjust = 1,
          inherit.aes = FALSE
        ) +
        geom_richtext(
          aes(x = threshold_upper, y = 0.95, 
              label = paste("<b>Mayor en", condicion1, "</b>")),
          color = color_mayor,
          size = 5,
          hjust = 0,
          vjust = 1,
          inherit.aes = FALSE
        )
    }
    
    return(p)
  }

