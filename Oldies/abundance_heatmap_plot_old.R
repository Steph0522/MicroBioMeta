#' Heatmap of relative abundance
#'
#'This function create a heatmap to visualize the relative abundance of ASV's, features o bacterial groups.
#'
#' @param table Data frame with taxonomy, where, the columns are the samples and rows are ASV's or taxa.
#' @param metadata Data frame of characteristics or important information of the samples
#' @param condition1 Variable of the first horizontal annotation
#' @param condition2 Variable of the second horizontal annotation
#' @param condition3 Variable of the third horizontal annotation
#' @param colors_condition1 Color vector for condition 1 
#' @param colors_condition2 Color vector for condition 2
#' @param colors_condition3 Color vector for condition 3
#' @param name_legend_condition1 Title assigned to legend of condition 1 
#' @param name_legend_condition2 Title assigned to legend of condition 2
#' @param name_legend_condition3 Title assigned to legend of condition 3
#' @param arguments_condition1 Vector of specific arguments for condition 1
#' @param arguments_condition2  Vector of specific arguments for condition 2
#' @param arguments_condition3  Vector of specific arguments for condition 3
#' @param top_n Number of features to plot.
#'
#' @return A plot with the fifty (XX) taxonomic groups most abundant. 
#' @export
#'
#' @examples abundance_heatmap_plot(table = table_taxonomy, 
#'                        metadata = metadata.gestacion.recto,
#'                        condition1 = "poblacion",
#'                        condition2 = "temporada",
#'                        condition3 = "sexo",
#'                        top_n = 50,
#'                        name_legend_condition1 = "Altitude",
#'                        name_legend_condition2 = "Season",
#'                        name_legend_condition3 = "Sex",
#'                        colors_condition1 = c("#676778","#D9D9C2"),
#'                        colors_condition2 = c("#0E6251", "#1B5E20"),
#'                        colors_condition3 = c("#5D3277","#AF6502"),
#'                        arguments_condition1 = c("2600 masl" = "#676778", "4150 masl"= "#D9D9C2"),
#'                        arguments_condition2 = c("Reproductive" = "#0E6251", "No Reproductive"= "#1B5E20"),
#'                        arguments_condition3 = c("Female" = "#5D3277", "Male" ="#AF6502"))
#' 
#' 
abundance_heatmap_plot <- function(table,
                                   metadata,
                                   condition1 = NULL,
                                   condition2 = NULL,
                                   condition3 = NULL,
                                   colors_condition1 = NULL,
                                   colors_condition2 = NULL,
                                   colors_condition3 = NULL,
                                   name_legend_condition1 = NULL,
                                   name_legend_condition2 = NULL,
                                   name_legend_condition3 = NULL,
                                   arguments_condition1 = NULL,
                                   arguments_condition2 = NULL,
                                   arguments_condition3 = NULL,
                                   top_n) {
  
  # Delete taxonomy column
  taxonomy_col <- ncol(table)
  table_counts <- table[, -taxonomy_col, drop = FALSE]
  
  OTU_ID <- colnames(metadata)[1]
  
  # Separate taxonomy
  table_tax <- table %>%
    tibble::rownames_to_column("OTUID") %>%
    dplyr::select(OTUID, taxonomy)
  
  # Initial table
  phy.ra.complete <- t(t(table_counts)/colSums(table_counts)*100) %>%
    as.data.frame()
  
  # Editing table
  table_abundance <- phy.ra.complete %>%
    mutate(abun = rowMeans(.)) %>% 
    tibble::rownames_to_column(var="OTUID") %>% 
    dplyr::arrange(-abun) %>%
    dplyr::slice(1:top_n) %>%
    dplyr::left_join(table_tax, by = "OTUID")  %>%
    mutate(taxonomy2 = taxonomy) %>%
    tidyr::separate(taxonomy2, into = c("dominio","phylum","clase","orden","familia","genero","especie"), 
                    sep = ";", fill = "right", extra = "merge") %>% 
    dplyr::mutate(taxonomy = dplyr::case_when(
      grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
      grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
      grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
      grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
      grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
      TRUE ~ "Unclassified")) %>% 
    mutate(asv=paste0("ASV", dplyr::row_number())) %>%
    tidyr::unite("taxa", asv, taxonomy, remove = F) %>%
    mutate(across(everything(), ~ trimws(.))) %>%
    mutate(phylum = sub("^p__", "", phylum)) %>%
    dplyr::mutate(phylum = case_when(phylum == "Proteobacteria" ~ "Pseudomonadata",
                                     phylum=="Firmicutes" ~ "Bacillota",  
                                     TRUE~as.character(phylum)))
  
  # Join table with metadata
  heat <- table_abundance %>% 
    tibble::column_to_rownames(var = "taxa") %>%
    t() %>%
    as.data.frame() %>% 
    dplyr::slice(-1) %>% 
    tibble::rownames_to_column(var = OTU_ID) %>%
    dplyr::inner_join(metadata, by = OTU_ID) 
  
  # Transformation of ranges
  heatmap <- heat %>%
    tibble::column_to_rownames(var = OTU_ID) %>%
    dplyr::select(starts_with("ASV")) %>% 
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(across(everything(), ~ case_when(
      . <= 0.001 ~ 0,
      . >  0.001 & .  <= 0.005 ~ 1,
      . >  0.005 & .  <= 0.01 ~ 2,
      . >  0.01 & .  <= 0.10 ~ 3,
      . >  0.10 & .  <= 0.20 ~ 4,
      . >  0.20 & .  <= 1.00 ~ 5,
      . >  1.00 & .  <= 2.00 ~ 6,
      . >  2.00 & .  <= 5.00 ~ 7,
      . >  5.00 & .  <= 10.00 ~ 8,
      . >  10.00 & .  <= 25.00 ~ 9,
      . >  25.00 & .  <= 50.00 ~ 10,
      . >  50.00 & .  <= 75.00 ~ 11,
      . >  75.00 ~ 12))) 
  
  heatmap <- as.matrix(heatmap)
  
  # Annotation row 
  annotation_rows <- table_abundance %>%
    dplyr::select(OTUID, phylum) %>%
    tibble::column_to_rownames(var = "OTUID")
  rownames(annotation_rows) <- rownames(heatmap)
  
  # Dynamic selection of conditions
  conditions <- purrr::compact(list(condition1, condition2, condition3))
  
  # Create annotation_columns only with provided conditions
  if (length(conditions) > 0) {
    annotation_columns <- heat %>%
      dplyr::select(all_of(unlist(conditions))) 
    rownames(annotation_columns) <- colnames(heatmap)
  } else {
    annotation_columns <- data.frame()
  }
  
  # Color palette for heatmap
  my_palette <- viridis::viridis(n = 12, option = "C", direction = -1)
  
  # Color for phylum annotations
  unique_phyla <- unique(annotation_rows$phylum)
  c5.phylum <- RColorBrewer::brewer.pal(max(3, length(unique_phyla)), "Set2")[seq_along(unique_phyla)]
  cols_phyl <- list(Phylum = setNames(c5.phylum, unique_phyla))
  
  # Phylum annotation (always present)
  annphylum <- ComplexHeatmap::rowAnnotation(
    Phylum = annotation_rows$phylum, 
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
    gp = grid::gpar(col = "white"),
    col = cols_phyl
  )
  
  # Initialize lists for dynamic annotations and legends
  heatmap_annotations <- list()
  legend_list <- list(lgd1 = ComplexHeatmap::Legend(
    at = unique(annotation_rows$phylum), 
    legend_gp = grid::gpar(fill = c5.phylum), 
    title = "Phylum", 
    labels_gp = grid::gpar(fontsize=11)))
  
  # Process conditions if they exist
  if (!is.null(condition1)) {
    # Usa colors_condition1 si se proporciona, de lo contrario usa viridis
    if (!is.null(colors_condition1)) {
      unique_vals <- sort(unique(annotation_columns[[condition1]]))
      arguments_condition1 <- setNames(colors_condition1, unique_vals)
    } else if (is.null(arguments_condition1)) {
      unique_vals <- sort(unique(annotation_columns[[condition1]]))
      arguments_condition1 <- setNames(
        viridis::viridis(length(unique_vals)),
        unique_vals
      )
    }
    
    heatmap_annotations$ann1 <- ComplexHeatmap::HeatmapAnnotation(
      df = annotation_columns[condition1],
      which = "column",
      col = setNames(list(arguments_condition1), condition1),
      show_legend = FALSE,
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
      gp = grid::gpar(col = "white")
    )
    legend_list$lgd2 <- ComplexHeatmap::Legend(
      at = names(arguments_condition1),
      legend_gp = grid::gpar(fill = arguments_condition1),
      title = if (!is.null(name_legend_condition1)) name_legend_condition1 else condition1,
      labels_gp = grid::gpar(fontsize=11)
    )
  }
  
  if (!is.null(condition2)) {
    # Usa colors_condition2 si se proporciona
    if (!is.null(colors_condition2)) {
      unique_vals <- sort(unique(annotation_columns[[condition2]]))
      arguments_condition2 <- setNames(colors_condition2, unique_vals)
    } else if (is.null(arguments_condition2)) {
      unique_vals <- sort(unique(annotation_columns[[condition2]]))
      arguments_condition2 <- setNames(
        viridis::inferno(length(unique_vals)),
        unique_vals
      )
    }  
    
    heatmap_annotations$ann2 <- ComplexHeatmap::HeatmapAnnotation(
      df = annotation_columns[condition2],
      which = "column",
      col = setNames(list(arguments_condition2), condition2),
      show_legend = FALSE,
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
      gp = grid::gpar(col = "white")
    )
    legend_list$lgd3 <- ComplexHeatmap::Legend(
      at = names(arguments_condition2),
      legend_gp = grid::gpar(fill = arguments_condition2),
      title = if (!is.null(name_legend_condition2)) name_legend_condition2 else condition2,
      labels_gp = grid::gpar(fontsize=11)
    )
  } 
  
  if (!is.null(condition3)) {
    # Usa colors_condition3 si se proporciona
    if (!is.null(colors_condition3)) {
      unique_vals <- sort(unique(annotation_columns[[condition3]]))
      if(length(colors_condition3) < length(unique_vals)) {
        # Completa los colores faltantes con viridis
        colors_condition3 <- c(colors_condition3, 
                               viridis::viridis(length(unique_vals) - length(colors_condition3)))
      }
      arguments_condition3 <- setNames(colors_condition3[1:length(unique_vals)], unique_vals)
    } else if (is.null(arguments_condition3)) {
      unique_vals <- sort(unique(annotation_columns[[condition3]]))
      arguments_condition3 <- setNames(
        viridis::plasma(length(unique_vals)),
        unique_vals
      )
    } else {
      # Verifica que todos los niveles tengan colores
      unique_vals <- sort(unique(annotation_columns[[condition3]]))
      missing_levels <- setdiff(unique_vals, names(arguments_condition3))
      if(length(missing_levels) > 0) {
        arguments_condition3 <- c(arguments_condition3,
                                  setNames(viridis::magma(length(missing_levels)), missing_levels))
      }
    }
    
    heatmap_annotations$ann3 <- ComplexHeatmap::HeatmapAnnotation(
      df = annotation_columns[condition3],
      which = "column",
      col = setNames(list(arguments_condition3), condition3),  # Corregido aquí
      show_legend = FALSE,
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
      gp = grid::gpar(col = "white")
    )
    legend_list$lgd4 <- ComplexHeatmap::Legend(
      at = names(arguments_condition3),
      legend_gp = grid::gpar(fill = arguments_condition3),
      title = if (!is.null(name_legend_condition3)) name_legend_condition3 else condition3,
      labels_gp = grid::gpar(fontsize=11)
    )
  }
  
  # Combine all annotations
  top_annotation <- do.call(c, unname(heatmap_annotations))
  
  # Combine all legends
  pd.legends <- do.call(ComplexHeatmap::packLegend, legend_list)
  
  # Create heatmap
  heats <- ComplexHeatmap::Heatmap(
    heatmap, 
    col = my_palette,
    heatmap_legend_param = list(
      direction = "horizontal",
      labels_gp = grid::gpar(fontsize = 7),
      legend_gp = grid::gpar(fontsize = 9),
      title = "Relative abundance (%)",
      title_position = "topcenter",
      at = c(0,1,2,3,5,8,10,25, 50, 100),
      break_dist = 1
    ),
    rect_gp = grid::gpar(col = "black", lwd = 0.5),    
    row_names_gp = grid::gpar(fontsize=7),
    column_names_gp = grid::gpar(fontsize=7),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_column_names = TRUE,
    show_heatmap_legend = TRUE, 
    top_annotation = if (length(heatmap_annotations) > 0) top_annotation else NULL,
    left_annotation = annphylum
  )
  
  heatmap_output <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heats,
      heatmap_legend_side = "top",
      annotation_legend_side = "right", 
      merge_legend = FALSE,
      annotation_legend_list = pd.legends
    )
  )
  
  grid::grid.newpage()
  grid::grid.draw(heatmap_output)
  return(invisible(heatmap_output))  
}