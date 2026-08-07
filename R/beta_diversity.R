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
#' @param color_facets_x Color vector for facet strips.
#' @param color_axis_x Named color vector for x-axis groups.
#' @param save_table Logical. If \code{TRUE}, saves the underlying turnover
#'   table to disk. Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"betadiv_turnover.txt"}.
#'
#' @return A ggplot2 figure with beta diversity partitions across conditions.
#' @export
#'
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "table_with_taxonomy.tsv", package = "MicroBioMeta")
#' table <- read.delim(table_path, skip = 1, comment.char = "", check.names = FALSE, row.names = 1)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE, comment.char = "")
#' colnames(metadata)[1] <- "OTUID"
#'
#' # comparison_condition1/2 match "value.x_vs_value.y" pairs created by the
#' # pairwise self-join (condition.x = site1's value, condition.y = site2's value)
#' beta_plot(
#'   table                 = table,
#'   metadata              = metadata,
#'   comparison_condition1 = c("Rizosphere_vs_Roots", "Rizosphere_vs_Non-Rizospheric"),
#'   comparison_condition2 = c("1_vs_1", "2_vs_2"),
#'   condition1.x          = "Type_of_soil.x",
#'   condition1.y          = "Type_of_soil.y",
#'   condition2.x          = "Treatment.x",
#'   condition2.y          = "Treatment.y",
#'   color_facets_x        = c("#5D478B", "#8B668B"),
#'   color_axis_x          = c("Roots" = "#2F4F4F", "Non-Rizospheric" = "#698B69")
#' )
#' }

beta_plot <- function(table, 
                      metadata, 
                      comparison_condition1, 
                      comparison_condition2,
                      condition1.x,
                      condition1.y,
                      condition2.x,
                      condition2.y,
                      color_facets_x,
                      color_axis_x,
                      save_table = FALSE,
                      table_filename = "betadiv_turnover.txt") {

  # Treat metadata's first column as the sample ID regardless of its original name
  colnames(metadata)[1] <- "OTUID"

  # Paso 1: Filtrar y transponer tabla OTU
  otu_filter <- table %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric))) != 0)
  
  if (nrow(otu_filter) == 0) stop("Error: OTU table is empty after filtering rows with sum=0.")

  otu_filter_t <- as.data.frame(t(otu_filter))
  message("Step 1: OTU table filtered and transposed - dimension: ", paste(dim(otu_filter_t), collapse = " x "))
  
  # Paso 2: Calcular diversidad beta (Hill numbers)
  beta_q_list <- list()
  for (q in c(0, 1, 2)) {
    beta_res <- tryCatch({
      hillR::hill_taxa_parti_pairwise(comm = otu_filter_t, q = q) %>%
        dplyr::mutate(Recambio = TD_beta - 1, q = q)
    }, error = function(e) {
      warning(paste("WARNING: Error computing hill_taxa_parti_pairwise with q =", q, ":", e$message))
      NULL
    })
    if (!is.null(beta_res)) beta_q_list[[as.character(q)]] <- beta_res
  }
  
  beta_total <- dplyr::bind_rows(beta_q_list)
  
  if (nrow(beta_total) == 0) stop("Error: Could not compute beta partitions (empty table).")

  message("Step 2: beta_total computed - rows: ", nrow(beta_total))
  
  # Paso 3: Unir con metadata
  beta_formato <- beta_total %>%
    dplyr::inner_join(metadata, by = c("site1" = "OTUID")) %>%
    dplyr::inner_join(metadata, by = c("site2" = "OTUID"))
  
  if (nrow(beta_formato) == 0) stop("Error: Could not join beta_total with metadata (empty table).")

  message("Step 3: beta_formato joined with metadata - rows: ", nrow(beta_formato))
  
  # Paso 4: Crear comparaciones y filtrar
  beta_final <- beta_formato %>%
    tidyr::unite("compar_condition1", dplyr::all_of(c(condition1.x, condition1.y)), sep = "_vs_", remove = FALSE) %>%
    tidyr::unite("compar_condition2", dplyr::all_of(c(condition2.x, condition2.y)), sep = "_vs_", remove = FALSE) %>%
    dplyr::filter(compar_condition1 %in% comparison_condition1)  %>%
    dplyr::filter(compar_condition2 %in% comparison_condition2) %>%
    dplyr::mutate(orden = dplyr::case_when(
      q == 0 ~ "q=0", q == 1 ~ "q=1", q == 2 ~ "q=2", TRUE ~ as.character(q)
    ))
  
  if (nrow(beta_final) == 0) stop("Error: After filtering comparisons, the table is empty.")

  message("Step 4: beta_final ready - rows: ", nrow(beta_final))

  if (save_table) {
    utils::write.table(beta_final, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # Paso 5: Crear grafico
  figura <- beta_final %>%
    ggpubr::ggboxplot(x = condition1.y, y = "Recambio", fill = condition1.y) +
    ggplot2::ylab("Proportion of ASVs turnover") +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggh4x::facet_grid2(stats::as.formula(paste("orden ~", condition1.x)), 
                       scales = "free_x", 
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x))) +  
    ggplot2::theme_bw(base_family = "serif") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 12, color = "black"),
                   axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black",
                                                        margin = ggplot2::margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
                   strip.text.x = ggplot2::element_text(size = 12, face = "bold", color = "white"),
                   strip.text.y = ggplot2::element_text(size = 12, face = "italic"),
                   legend.text = ggplot2::element_text(size = 12, color = "black"),
                   legend.title = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", fill=NA, size=0.5)) +
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
                        color_axis_x,
                        save_table = FALSE,
                        table_filename = "betadiv_shared.txt")
{
  # Treat metadata's first column as the sample ID regardless of its original name
  colnames(metadata)[1] <- "OTUID"

  #Obtener base de datos sin singletons por muestra
  asv_table <- table 
  asv_table[asv_table>0]=1 
  asv_no_single<-asv_table %>% 
    dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric)))>1) %>%
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
    tidyr::unite("compar_condition1", dplyr::all_of(c(condition1.x, condition1.y)), sep="_vs_", remove=F) %>%
    tidyr::unite("compar_condition2", dplyr::all_of(c(condition2.x, condition2.y)), sep="_vs_", remove=F) %>%
    dplyr::filter(compar_condition1 %in% comparison_condition1) %>%
    dplyr::filter(compar_condition2 %in% comparison_condition2)
  
  #Calculo con base de datos sin singletones
  asv_no_single$sum <- rowSums(asv_no_single) 
  asv_emb<- asv_no_single %>% tibble::rownames_to_column(var="site2")
  # Toma la columna 'sum' de la ultima columna
  asv_emb <- asv_emb[, c("site2", "sum")]
  beta.shared.final2 <- dplyr::left_join(beta.shared.final, asv_emb, by="site2")
  beta.shared.final2$overlap<- (beta.shared.final2$value / beta.shared.final2$sum)*100

  if (save_table) {
    utils::write.table(beta.shared.final2, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  #Boxplot
  
  beta.shared<-beta.shared.final2 %>% 
    ggpubr::ggboxplot(x = condition1.y, y="overlap", fill=condition1.y, facet.by = condition1.x) +
    ggplot2::ylab("Percentage of shared ASVs") +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggh4x::facet_grid2(stats::as.formula(paste(". ~", condition1.x)), 
                       scales = "free_x", 
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x))) +  
    ggplot2::theme_bw(base_family = "serif") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 12, color = "black"),
                   axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black", margin = ggplot2::margin(t=0, r=0.5, b=0, l=0, "cm")),
                   strip.text.x = ggplot2::element_text(size = 12, face = "bold", color = "white"),
                   legend.text = ggplot2::element_text(size = 12, color = "black"),
                   legend.title = ggplot2::element_blank())+
    ggplot2::xlab(ggplot2::element_blank())+
    ggpubr::stat_compare_means(
      mapping = ggplot2::aes(
        label = scales::label_pvalue(accuracy = 0.001)(ggplot2::after_stat(p))
      ),
      label.x = 2, size = 3.5, family = "serif")
  
  
  
  return(beta.shared)
  
}



###FUNCIoN COMPLETA

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
                               family = c("sorensen", "jaccard"),
                               save_table = FALSE,
                               table_filename = "betadiv_partition.txt") {

  # Treat metadata's first column as the sample ID regardless of its original name
  colnames(metadata)[1] <- "OTUID"

  # Validar argumentos
  partition <- match.arg(partition)
  family <- match.arg(family)
  
  message("Filtering table (removing singletons)...")
  # Quitar singletons por muestra
  asv_table <- table
  asv_table[asv_table > 0] <- 1
  asv_no_single <- asv_table %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric))) > 1) %>%
    t() %>%
    as.data.frame()
  
  message("Computing beta diversity partition: ", partition, " (family = ", family, ")")
  
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
    tidyr::drop_na() %>%
    dplyr::filter(!value == 0)
  
  message("Joining with metadata...")
  beta_format <- beta_df %>%
    dplyr::inner_join(metadata, by = c("site1" = "OTUID")) %>%
    dplyr::inner_join(metadata, by = c("site2" = "OTUID"))
  
  # Crear comparaciones
  beta_final <- beta_format %>%
    tidyr::unite("compar_condition1", dplyr::all_of(c(condition1.x, condition1.y)), sep = "_vs_", remove=FALSE) %>%
    tidyr::unite("compar_condition2", dplyr::all_of(c(condition2.x, condition2.y)), sep = "_vs_", remove=FALSE) %>%
    dplyr::filter(compar_condition1 %in% comparison_condition1) %>%
    dplyr::filter(compar_condition2 %in% comparison_condition2)
  
  message("Filtered: ", nrow(beta_final), " final rows.")

  if (save_table) {
    utils::write.table(beta_final, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # Crear plot
  figura <- beta_final %>%
    ggpubr::ggboxplot(x = condition1.y, y = "value", fill = condition1.y, facet.by = condition1.x) +
    ggplot2::ylab(paste("Beta diversity (", partition, ")", sep="")) +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggh4x::facet_grid2(
      stats::as.formula(paste(". ~", condition1.x)),
      scales = "free_x",
      strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x))) +
    ggplot2::theme_bw(base_family = "serif") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 12, color = "black"),
                   axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black",
                                               margin = ggplot2::margin(t=0, r=0.5, b=0, l=0, "cm")),
                   strip.text.x = ggplot2::element_text(size = 12, face = "bold", color = "white"),
                   legend.text = ggplot2::element_text(size = 12, color = "black"),
                   legend.title = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = "black", fill=NA, size=0.5)) +
    ggplot2::xlab(title_axis_x)
  
  # Devolver plot + tabla
  return(list(plot = figura, data = beta_final))
}
