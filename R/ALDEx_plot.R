
#' Heatmap ALDEx
#'
#' @param table, Data frame with taxonomy, where, the columns are the samples and rows are ASV's or taxa.
#' @param conditions, Vector that defines the categories or classes to compare. 
#' @param effect, Effect size (default: >= 0.8)
#' @param pvalue_BH, Value p-ajusted (optional) 
#'
#' @return A plot with the deferentially abundant taxonomic groups between two categories or groups of samples.
#' @export
#'
#' @examples aldex_plot(table = table, 
#'                      conditions = conditions,
#'                      effect = 0.8,
#'                      pvalue_BH = NULL)
#' 


  aldex_plot <- function(table,
                       conditions,
                       effect = 0.8, 
                       pvalue_BH = NULL) {
  
  
  #Verify columns
  if (!"OTUID" %in% colnames(table)) stop("La tabla debe contener una columna 'OTUID'.")
  if (!"taxonomy" %in% colnames(table)) stop("La tabla debe contener una columna 'taxonomy'.")
  
  
  #Delete taxonomy
  table_counts <- table %>%
    dplyr::select(-taxonomy) %>%
    tibble::column_to_rownames("OTUID")
  
  #Separate taxonomy
  table_tax <- table %>% dplyr::select(OTUID, taxonomy)
  
  #Verify conditions 
  if (length(conditions) != ncol(table_counts)) {
    stop("El número de condiciones no coincide con el número de columnas de abundancia.")
  }
  
  #Run ALDEX
  aldex_results <- ALDEx2::aldex(reads = table_counts,
                                 conditions = conditions,
                                 mc.samples = 1000,
                                 effect = TRUE, 
                                 test = "t",
                                 verbose = TRUE,
                                 denom = "all",
                                 include.sample.summary = FALSE)
  
  
  #Filter by effect size (default)
  aldex.filtered.effect <- aldex_results %>% dplyr::filter(abs(effect) >= effect)
  
  #Filter table by value p if specified
  if (!is.null(pvalue_BH)) {
    aldex.filtered.effect <- aldex.filtered.effect %>% dplyr::filter(wi.eBH <= pvalue_BH)
  }
  
  #Declare and change for the plot 
  aldex.plot<- aldex.filtered.effect %>% 
    tibble::rownames_to_column("OTUID") %>%
    dplyr::left_join(table_tax, by = "OTUID") %>%
    dplyr::mutate(seccion = dplyr::case_when( #declarar condiciones de acuerdo al valor negativo o positivo de diff.btw
      diff.btw < 0 ~ "NonReproductive.female",
      diff.btw > 0 ~ "Reproductive.female"),
      taxonomy = dplyr::case_when(
        grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
        grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
        grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
        grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
        grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
        TRUE ~ "Unclassified"),
      taxonomy = stringr::str_trim(taxonomy),  # Elimina espacios al inicio/final
      taxonomy = make.unique(taxonomy), 
      p.value = dplyr::case_when( #declarar cambio de valor p (wi.eBH) por rangos establecidos 
        wi.eBH <= 0.001 ~ "<0.001",
        wi.eBH <= 0.01 ~ "<0.01",
        wi.eBH <  0.05 ~ "<0.05", TRUE ~ ">0.05")) %>%
    dplyr::arrange(diff.btw)   #ordenar de mayor a menor 
  
  #Heatmap table
  heat.data <- aldex.plot %>%
    dplyr::select(taxonomy, "Reproductive.female"= rab.win.Reproductive.female,
                            "NonReproductive.female"= rab.win.NonReproductive.female) %>% #Cambia el nombre de las columnas rab.win.ALTA y rab.win.BAJA
    tibble::column_to_rownames(var = "taxonomy") %>%
    as.matrix()  
  
  #Declare heatmap size
  heatmap_width <- grid::unit(ncol(heat.data) * 7, "mm")
  heatmap_height <- grid::unit(nrow(heat.data) * 6, "mm")
  
  #Figure annotation
  #Bar colors
  treatment.colors = c("Reproductive.female" = "#808000", "NonReproductive.female" = "#1B5E20")
  
  # Bar annotation 
  barpl = ComplexHeatmap::rowAnnotation("difference \nbetween groups" = ComplexHeatmap::anno_barplot(aldex.plot$diff.btw,
                                                                                                     which = "row",
                                                                                                     gp = grid::gpar(fill = treatment.colors[aldex.plot$seccion]),
                                                                                                     width = unit(4, "cm")),  
                                                                                                     show_annotation_name = T,
                                                                                                     annotation_name_gp = grid::gpar(fontsize = 0),
                                                                                                     annotation_name_rot = 0)
  
  #value p colors
  pvalue.colors <- list('p-value' = c(
    "<0.001" = '#C70039',
    "<0.01"='#FF5733',
    "<0.05"="#FFC300",
    ">0.05"="#F9E79F"))
  
  #valor p annotation
  annP = ComplexHeatmap::rowAnnotation("p-value" =aldex.plot$p.value, 
                                       simple_anno_size = unit(0.45, "cm"),
                                       annotation_name_gp = grid::gpar(fontsize = 8, fontface="bold"),
                                       annotation_legend_param = list(
                                         title_gp = grid::gpar(fontsize = 8, fontface="bold"),
                                         labels_gp = grid::gpar(fontsize = 8),
                                         direction ="vertical"),
                                       col = pvalue.colors,
                                       show_legend = TRUE,
                                       gp = grid::gpar(col = "white"), 
                                       show_annotation_name = TRUE)
  
  
  #effect colors
  effect.colors = circlize::colorRamp2(c(-1.5, 0, 1.5), c("lightsalmon4", "white", "lightseagreen"))
  
  #Effect size annotation
  annE = ComplexHeatmap::rowAnnotation("Effect size" = aldex.plot$effect, 
                                       col = list("Effect size" = effect.colors),
                                       simple_anno_size = unit(0.45, "cm"),
                                       annotation_name_gp = grid::gpar(fontsize = 8, fontface="bold"), 
                                       annotation_legend_param = list(title_gp = grid::gpar(fontsize = 8,fontface="bold"),
                                                                      labels_gp = grid::gpar(fontsize = 8),
                                                                      direction ="vertical"),
                                       show_legend = TRUE, 
                                       gp = grid::gpar(col = "white"),
                                       show_annotation_name = TRUE)
  
  
  #heatmap colors
  heatmap.colors = circlize::colorRamp2(seq(min(heat.data), max(heat.data), length= 3),
                                        c("#cbdbcd", "#759c8a", "#00544d"))
  
  
  #Final heatmap
  heatmap <- ComplexHeatmap::Heatmap(heat.data, 
                                     cluster_rows = F,
                                     cluster_columns = F,
                                     width = heatmap_width, 
                                     height = heatmap_height,
                                     column_names_rot = 90,
                                     rect_gp = grid::gpar(col = "white", lwd = 2),
                                     left_annotation = annE + annP,
                                     right_annotation = barpl,
                                     name = "Median clr value", 
                                     heatmap_legend_param = list(direction = "vertical" , 
                                                                 labels_gp = grid::gpar(fontsize = 8),
                                                                 title_gp = grid::gpar(fontsize = 8,  fontface="bold"),
                                                                 legend_height = unit(1.4, "cm")), 
                                     column_names_gp = grid::gpar(fontsize = 8, fontface="bold", direction="vertical"),
                                     col = heatmap.colors,
                                     row_names_gp = grid::gpar(fontsize = 8),
                                     show_heatmap_legend = T)
  
  return(heatmap) 
  
}





