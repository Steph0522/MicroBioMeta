#' Box plot of beta diversity
#'
#' @param table table Data frame where columns are samples and rows are ASVs or taxa.
#' @param metadata  A data frame with sample metadata. The first column must match sample names in "table".
#' @param comparison_condition1 Vector of comparisons 1
#' @param comparison_condition2 Vector of comparisons 2
#' @param condition1.x Condition 1 of axis-x
#' @param condition1.y Condition 1 of axis-y
#' @param condition2.x Condition 2 of axis-x
#' @param condition2.y Condition 2 of axis-y
#' @param color_facets_x Color vector for facet
#' @param color_axis_x Color and conditions vector for axis-x
#' @param title_axis_x  Axis-x title 
#' @param partition Component of beta diversity, e.g. "turnover", "nestedness" or "shared"
#' @param family  Family for turnover o nestednes, e.g. "sorensen" or "jaccard".
#'
#' @return A plot with the values of beta diversity. 
#' @export
#'
#' @examples    res <- beta_plot_flexible(
#'                        table = table,
#'                        metadata = metadata,
#'                        comparison_condition1 = c("Boca_vs_L.amniotico","Boca_vs_Membrana","Boca_vs_Yema","Boca_vs_Tracto.embrionario",
#'                                                  "Cloaca_vs_L.amniotico","Cloaca_vs_Membrana","Cloaca_vs_Yema","Cloaca_vs_Tracto.embrionario",
#'                                                  "Ileon_vs_L.amniotico","Ileon_vs_Membrana","Ileon_vs_Yema","Ileon_vs_Tracto.embrionario",
#'                                                  "Dorso_vs_L.amniotico", "Dorso_vs_Membrana", "Dorso_vs_Yema", "Dorso_vs_Tracto.embrionario"),
#'                        comparison_condition2 = c("3_vs_3", "7_vs_7", "16_vs_16", "8_vs_8", "12_vs_12"),
#'                        condition1.x = "Seccion.x",
#'                        condition1.y = "Seccion.y",
#'                        condition2.x = "ID.identificador.x",
#'                        condition2.y = "ID.identificador.y",
#'                        color_facets_x = c("#5D478B", "#8B668B", "#CDB5CD"),
#'                        color_axis_x = c("L.amniotico"="#2F4F4F", "Tracto.embrionario"="#698B69", "Membrana"="#458B74", "Yema"="#B4EEB4"),
#'                        partition = "turnover",
#'                        family = "jaccard",
#'                        title_axis_x = "Maternal samples"
#'                              )
#' 
#'                   # Show plot
#'                     res$plot
#' 
#' 
#' #FUNCION RECAMBIO
beta_plot <- function(table, 
                      metadata, 
                      comparison_condition1, 
                      comparison_condition2,
                      condition1.x,
                      condition1.y,
                      condition2.x,
                      condition2.y,
                      color_facets_x,
                      color_axis_x) {
  
  
  # Paso 1: Filtrar y transponer tabla OTU
  otu_filter <- table %>%
    dplyr::filter(rowSums(across(where(is.numeric))) != 0)
  
  if (nrow(otu_filter) == 0) stop("Error: La tabla OTU quedó vacía después de filtrar filas con suma=0.")
  
  otu_filter_t <- as.data.frame(t(otu_filter))
  message("✅ Paso 1: OTU filtrada y transpuesta - dimensión: ", paste(dim(otu_filter_t), collapse = " x "))
  
  # Paso 2: Calcular diversidad beta (Hill numbers)
  beta_q_list <- list()
  for (q in c(0, 1, 2)) {
    beta_res <- tryCatch({
      hillR::hill_taxa_parti_pairwise(comm = otu_filter_t, q = q) %>%
        dplyr::mutate(Recambio = TD_beta - 1, q = q)
    }, error = function(e) {
      warning(paste("⚠️ Error calculando hill_taxa_parti_pairwise con q =", q, ":", e$message))
      NULL
    })
    if (!is.null(beta_res)) beta_q_list[[as.character(q)]] <- beta_res
  }
  
  beta_total <- bind_rows(beta_q_list)
  
  if (nrow(beta_total) == 0) stop("Error: No se pudieron calcular las particiones beta (tabla vacía).")
  
  message("✅ Paso 2: Beta_total calculada - filas: ", nrow(beta_total))
  
  # Paso 3: Unir con metadata
  beta_formato <- beta_total %>%
    dplyr::inner_join(metadata, by = c("site1" = "OTUID")) %>%
    dplyr::inner_join(metadata, by = c("site2" = "OTUID"))
  
  if (nrow(beta_formato) == 0) stop("Error: No se pudo unir beta_total con metadata (tabla vacía).")
  
  message("✅ Paso 3: Beta_formato unido con metadata - filas: ", nrow(beta_formato))
  
  # Paso 4: Crear comparaciones y filtrar
  beta_final <- beta_formato %>%
    tidyr::unite("compar_condition1", dplyr::all_of(c(condition1.x, condition1.y)), sep = "_vs_", remove = FALSE) %>%
    tidyr::unite("compar_condition2", dplyr::all_of(c(condition2.x, condition2.y)), sep = "_vs_", remove = FALSE) %>%
    dplyr::filter(compar_condition1 %in% comparison_condition1)  %>%
    dplyr::filter(compar_condition2 %in% comparison_condition2) %>%
    dplyr::mutate(orden = dplyr::case_when(
      q == 0 ~ "q=0", q == 1 ~ "q=1", q == 2 ~ "q=2", TRUE ~ as.character(q)
    ))
  
  if (nrow(beta_final) == 0) stop("Error: Después de filtrar comparaciones, la tabla quedó vacía.")
  
  message("✅ Paso 4: Beta_final lista - filas: ", nrow(beta_final))
  
  # Paso 5: Crear gráfico
  figura <- beta_final %>%
    ggpubr::ggboxplot(x = condition1.y, y = "Recambio", fill = condition1.y) +
    ggplot2::ylab("Proportion of ASVs turnover") +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggh4x::facet_grid2(stats::as.formula(paste("orden ~", condition1.x)), 
                       scales = "free_x", 
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x))) +  
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 8),
                   axis.title.y = ggplot2::element_text(size = 12, face = "bold", 
                                                        margin = ggplot2::margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
                   strip.text.x = ggplot2::element_text(size = 7, face = "bold", color = "white"),
                   strip.text.y = ggplot2::element_text(face = "italic"),
                   legend.text = ggplot2::element_text(size = 5),
                   legend.title = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", fill=NA, size=0.5)) + # <<--- añade esto
    ggplot2::xlab("Section")
  
  # Devuelve plot
  return(list(plot = figura))
}


###FUNCION SHARED 

shared_plot <- function(table, 
                        metadata, 
                        comparison_condition1, 
                        comparison_condition2,
                        condition1.x,
                        condition1.y,
                        condition2.x,
                        condition2.y,
                        color_facets_x,
                        color_axis_x) 
{
  
  #Obtener base de datos sin singletons por muestra
  asv_table <- table 
  asv_table[asv_table>0]=1 
  asv_no_single<-asv_table %>% 
    dplyr::filter(rowSums(across(where(is.numeric)))>1) %>%
    t() %>%
    as.data.frame() 
  
  #Obtener el objeto core
  core.beta<- betapart::betapart.core (asv_no_single)
  shared<-core.beta$shared
  
  beta.shared <- as.matrix(shared) %>% 
    reshape2::melt(varnames = c("site1", "site2")) %>%
    tidyr::drop_na() %>%
    dplyr::filter(!value==0)
  
  #Unir indices de diversidad beta con metadata
  beta.shared.formato<- beta.shared %>%
    dplyr::inner_join(metadata, by = c("site1"="OTUID")) %>%
    dplyr::inner_join(metadata, by = c("site2"="OTUID")) 
  
  
  #Subset para obtener comparaciones anteriores
  beta.shared.final<-beta.shared.formato %>%
    tidyr::unite("compar_condition1", all_of(c(condition1.x, condition1.y)), sep="_vs_", remove=F) %>%
    tidyr::unite("compar_condition2", all_of(c(condition2.x, condition2.y)), sep="_vs_", remove=F) %>%
    dplyr::filter(compar_condition1 %in% comparison_condition1) %>%
    dplyr::filter(compar_condition2 %in% comparison_condition2)
  
  #Calcúlo con base de datos sin singletones
  asv_no_single$sum <- rowSums(asv_no_single) 
  asv_emb<- asv_no_single %>% tibble::rownames_to_column(var="site2")
  # Toma la columna 'sum' de la última columna
  asv_emb <- asv_emb[, c("site2", "sum")]
  beta.shared.final2 <- dplyr::left_join(beta.shared.final, asv_emb, by="site2")
  beta.shared.final2$overlap<- (beta.shared.final2$value / beta.shared.final2$sum)*100   
  
  
  #Boxplot
  
  beta.shared<-beta.shared.final2 %>% 
    ggpubr::ggboxplot(x = condition1.y, y="overlap", fill=condition1.y, facet.by = condition1.x) +
    ggplot2::ylab("Percentage of shared ASVs") +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggh4x::facet_grid2(stats::as.formula(paste(". ~", condition1.x)), 
                       scales = "free_x", 
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x))) +  
    ggplot2::theme(axis.text.x = element_blank(),
                   axis.ticks.x =element_blank(),
                   axis.text.y = element_text(size = 10),
                   axis.title.y = element_text(size = 12, face = "bold", margin =margin(t=0, r=0.5, b=0, l=0, "cm")),
                   strip.text.x = element_text(size = 10, face = "bold", color = "white"),
                   legend.title = element_blank())+
    ggplot2::xlab(element_blank())+
    ggpubr::stat_compare_means(label="p.format", label.x = 2, size=2.8)
  
  
  
  return(beta.shared)
  
}



###FUNCIÓN COMPLETA

beta_plot_flexible <- function(table, 
                               metadata, 
                               comparison_condition1, 
                               comparison_condition2,
                               condition1.x, 
                               condition1.y, 
                               condition2.x, 
                               condition2.y,
                               color_facets_x, 
                               color_axis_x,
                               title_axis_x,
                               partition = c("shared", "turnover", "nestedness"),
                               family = c("sorensen", "jaccard")) {
  
  # Validar argumentos
  partition <- match.arg(partition)
  family <- match.arg(family)
  
  message("✅ Filtrando tabla (quitar singletons)...")
  # Quitar singletons por muestra
  asv_table <- table
  asv_table[asv_table > 0] <- 1
  asv_no_single <- asv_table %>%
    filter(rowSums(across(where(is.numeric))) > 1) %>%
    t() %>%
    as.data.frame()
  
  message("✅ Calculando partición de beta diversity: ", partition, " (family = ", family, ")")
  
  if (partition == "shared") {
    # Betapart.core da matriz shared
    core.beta <- betapart::betapart.core(asv_no_single)
    shared <- core.beta$shared
    beta_mat <- shared
  } else {
    # beta.pair da list con turnover, nestedness o beta
    beta_pair <- betapart::beta.pair(asv_no_single, index.family = family)
    if (partition == "turnover") {
      beta_mat <- beta_pair[[1]] # turnover
    } else if (partition == "nestedness") {
      beta_mat <- beta_pair[[2]] # nestedness
    }
  }
  
  # Convertir a long format
  beta_df <- as.matrix(beta_mat) %>%
    reshape2::melt(varnames = c("site1", "site2")) %>%
    drop_na() %>%
    dplyr::filter(!value == 0)
  
  message("✅ Uniendo con metadata...")
  beta_format <- beta_df %>%
    dplyr::inner_join(metadata, by = c("site1" = "OTUID")) %>%
    dplyr::inner_join(metadata, by = c("site2" = "OTUID"))
  
  # Crear comparaciones
  beta_final <- beta_format %>%
    tidyr::unite("compar_condition1", all_of(c(condition1.x, condition1.y)), sep = "_vs_", remove=FALSE) %>%
    tidyr::unite("compar_condition2", all_of(c(condition2.x, condition2.y)), sep = "_vs_", remove=FALSE) %>%
    dplyr::filter(compar_condition1 %in% comparison_condition1) %>%
    dplyr::filter(compar_condition2 %in% comparison_condition2)
  
  message("✅ Filtrado: ", nrow(beta_final), " filas finales.")
  
  # Crear plot
  figura <- beta_final %>%
    ggpubr::ggboxplot(x = condition1.y, y = "value", fill = condition1.y, facet.by = condition1.x) +
    ggplot2::ylab(paste("Beta diversity (", partition, ")", sep="")) +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggh4x::facet_grid2(
      stats::as.formula(paste(". ~", condition1.x)),
      scales = "free_x",
      strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x))) +
    ggplot2::theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.y = element_text(size = 8),
                   axis.title.y = element_text(size = 12, face = "bold", 
                                               margin = margin(t=0, r=0.5, b=0, l=0, "cm")),
                   strip.text.x = element_text(size = 7, face = "bold", color = "white"),
                   legend.text = element_text(size = 5),
                   legend.title = element_blank(),
                   panel.border = element_rect(color = "black", fill=NA, size=0.5)) +
    ggplot2::xlab(title_axis_x)
  
  # Devolver plot + tabla
  return(list(plot = figura, data = beta_final))
}
