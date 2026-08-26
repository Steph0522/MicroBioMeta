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
#' @param x_axis_title Title for the x-axis. Default \code{"Section"}.
#' @param show_x_labels Logical. If \code{TRUE}, x-axis tick labels are shown.
#'   Default \code{FALSE}, since the same groups are already named in the legend.
#' @param x_label_angle Numeric. Rotation (in degrees) of the x-axis tick labels
#'   when \code{show_x_labels = TRUE}. Default \code{0} (horizontal).
#' @param strip_text_bold Logical. If \code{TRUE}, facet strip labels are bold.
#'   Default \code{FALSE} (plain).
#' @param strip_text_color Color of the top (x) facet strip labels, which sit on
#'   the \code{color_facets_x} backgrounds. Default \code{"white"}.
#' @param aspect_ratio Numeric. Aspect ratio (height/width) of each panel.
#'   Default \code{NULL} (automatic).
#' @param stat Character or \code{NULL}. Statistical test to compare groups
#'   within each panel, passed to \code{ggpubr::stat_compare_means()} (e.g.
#'   \code{"wilcox.test"}, \code{"kruskal.test"}, \code{"anova"}). Default
#'   \code{NULL} (no test shown).
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
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "OTUID"
#'
#' # comparison_condition1/2 match "value.x_vs_value.y" pairs created by the
#' # pairwise self-join (condition.x = site1's value, condition.y = site2's value)
#' beta_turnover_plot(
#'   table                 = table,
#'   metadata              = metadata,
#'   comparison_condition1 = c("Rhizosphere_vs_Roots", "Rhizosphere_vs_Rhizosphere"),
#'   comparison_condition2 = c("Control_vs_Control", "Moderate_drought_vs_Moderate_drought"),
#'   condition1.x          = "Location.x",
#'   condition1.y          = "Location.y",
#'   condition2.x          = "Treatment.x",
#'   condition2.y          = "Treatment.y",
#'   color_facets_x        = c("#5D478B", "#8B668B"),
#'   color_axis_x          = c("Roots" = "#56B4E9", "Rhizosphere" = "#E69F00")
#' )
#' }

beta_turnover_plot <- function(table, 
                      metadata, 
                      comparison_condition1, 
                      comparison_condition2,
                      condition1.x,
                      condition1.y,
                      condition2.x,
                      condition2.y,
                      color_facets_x,
                      color_axis_x,
                      x_axis_title = "Section",
                      show_x_labels = FALSE,
                      x_label_angle = 0,
                      strip_text_bold = FALSE,
                      strip_text_color = "white",
                      aspect_ratio = NULL,
                      stat = NULL,
                      save_table = FALSE,
                      table_filename = "betadiv_turnover.txt") {

  # Treat metadata's first column as the sample ID regardless of its original name
  colnames(metadata)[1] <- "OTUID"

  # Drop a taxonomy column if present, so it isn't treated as a sample
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if (length(tax_col) == 1) table <- table[, -tax_col, drop = FALSE]

  # Paso 1: Filtrar y transponer tabla OTU
  otu_filter <- table %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric))) != 0)
  
  if (nrow(otu_filter) == 0) stop("Error: OTU table is empty after filtering rows with sum=0.")

  otu_filter_t <- as.data.frame(t(otu_filter))

  # Paso 2: Calcular diversidad beta (Hill numbers)
  beta_q_list <- list()
  for (q in c(0, 1, 2)) {
    beta_res <- tryCatch({
      # hill_taxa_parti_pairwise() prints a raw pairwise-comparison progress
      # bar straight to stdout (not via message()/warning()), which
      # suppressMessages() can't catch - capture.output() discards it while
      # still returning the function's actual result.
      result <- NULL
      utils::capture.output(
        result <- hillR::hill_taxa_parti_pairwise(comm = otu_filter_t, q = q)
      )
      result %>%
        dplyr::mutate(Recambio = TD_beta - 1, q = q)
    }, error = function(e) {
      warning(paste("WARNING: Error computing hill_taxa_parti_pairwise with q =", q, ":", e$message))
      NULL
    })
    if (!is.null(beta_res)) beta_q_list[[as.character(q)]] <- beta_res
  }
  
  beta_total <- dplyr::bind_rows(beta_q_list)
  
  if (nrow(beta_total) == 0) stop("Error: Could not compute beta partitions (empty table).")


  # Paso 3: Unir con metadata
  beta_formato <- beta_total %>%
    dplyr::inner_join(metadata, by = c("site1" = "OTUID")) %>%
    dplyr::inner_join(metadata, by = c("site2" = "OTUID"))
  
  if (nrow(beta_formato) == 0) stop("Error: Could not join beta_total with metadata (empty table).")


  # Paso 4: Crear comparaciones y filtrar
  beta_final <- beta_formato %>%
    tidyr::unite("compar_condition1", dplyr::all_of(c(condition1.x, condition1.y)), sep = "_vs_", remove = FALSE) %>%
    tidyr::unite("compar_condition2", dplyr::all_of(c(condition2.x, condition2.y)), sep = "_vs_", remove = FALSE) %>%
    dplyr::filter(compar_condition1 %in% comparison_condition1)  %>%
    dplyr::filter(compar_condition2 %in% comparison_condition2) %>%
    dplyr::mutate(orden = dplyr::case_when(
      q == 0 ~ "q0", q == 1 ~ "q1", q == 2 ~ "q2", TRUE ~ as.character(q)
    ))
  
  if (nrow(beta_final) == 0) stop("Error: After filtering comparisons, the table is empty.")

  if (save_table) {
    utils::write.table(beta_final, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # Paso 5: Crear grafico
  q_labeller <- .mbm_q_labeller(strip_text_bold)

  # The x-axis groups are already named in the legend, so their tick labels are
  # hidden by default; `show_x_labels = TRUE` brings them back (rotated by
  # `x_label_angle`, as in the other bar/boxplot functions).
  x_text  <- if (show_x_labels) .mbm_x_text(x_label_angle) else ggplot2::element_blank()
  x_ticks <- if (show_x_labels) ggplot2::element_line(colour = "black") else ggplot2::element_blank()

  figura <- beta_final %>%
    ggpubr::ggboxplot(x = condition1.y, y = "Recambio", fill = condition1.y) +
    ggplot2::ylab("Proportion of ASVs turnover") +
    ggplot2::scale_fill_manual(values = color_axis_x) +
    ggh4x::facet_grid2(stats::as.formula(paste("orden ~", condition1.x)),
                       scales = "free_x",
                       labeller = ggplot2::labeller(orden = q_labeller),
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = color_facets_x))) +
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        panel.grid   = ggplot2::element_blank(),
        axis.text.x  = x_text,
        axis.ticks.x = x_ticks,
        axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black",
                                             margin = ggplot2::margin(t = 0, r = 0.5, b = 0, l = 0, "cm")),
        # x strips sit on the user-supplied `color_facets_x` backgrounds (often
        # dark), so their text color is exposed separately from the y strips,
        # which keep the package's plain black-on-white look.
        strip.text.x = .mbm_strip_text(strip_text_bold, colour = strip_text_color),
        strip.text.y = .mbm_strip_text(strip_text_bold),
        legend.title = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    ) +
    ggplot2::xlab(x_axis_title)

  if (!is.null(aspect_ratio)) {
    figura <- figura + ggplot2::theme(aspect.ratio = aspect_ratio)
  }

  if (!is.null(stat)) {
    figura <- figura + ggpubr::stat_compare_means(
      method = stat,
      mapping = ggplot2::aes(
        label = paste0("p = ", scales::label_pvalue(accuracy = 0.001)(ggplot2::after_stat(p)))
      ),
      size = 3.5, family = "serif", hide.ns = TRUE
    )
  }

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

  # Drop a taxonomy column if present, so it isn't treated as a sample
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if (length(tax_col) == 1) table <- table[, -tax_col, drop = FALSE]

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
    tidyr::unite("compar_condition1", dplyr::all_of(c(condition1.x, condition1.y)), sep="_vs_", remove=FALSE) %>%
    tidyr::unite("compar_condition2", dplyr::all_of(c(condition2.x, condition2.y)), sep="_vs_", remove=FALSE) %>%
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
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        panel.grid   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black", margin = ggplot2::margin(t=0, r=0.5, b=0, l=0, "cm")),
        strip.text.x = .mbm_strip_text(colour = "white"),
        legend.title = ggplot2::element_blank()
      )
    ) +
    ggplot2::xlab(ggplot2::element_blank())+
    ggpubr::stat_compare_means(
      mapping = ggplot2::aes(
        label = paste0("p = ", scales::label_pvalue(accuracy = 0.001)(ggplot2::after_stat(p)))
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

  # Drop a taxonomy column if present, so it isn't treated as a sample
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if (length(tax_col) == 1) table <- table[, -tax_col, drop = FALSE]

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
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        panel.grid   = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black",
                                             margin = ggplot2::margin(t=0, r=0.5, b=0, l=0, "cm")),
        strip.text.x = .mbm_strip_text(colour = "white"),
        legend.title = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
    ) +
    ggplot2::xlab(title_axis_x)
  
  # Devolver plot + tabla
  return(list(plot = figura, data = beta_final))
}
