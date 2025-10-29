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
#' @param top_n Number of features to plot.
#' @param cluster Logical indicating whether to cluster rows (TRUE) or order by abundance (FALSE)
#' @param show_column_names Logical indicating whether to show column names (TRUE) or not (FALSE)
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
#'                        cluster = TRUE,
#'                        show_column_names = FALSE,
#'                        name_legend_condition1 = "Altitude",
#'                        name_legend_condition2 = "Season",
#'                        name_legend_condition3 = "Sex",
#'                        colors_condition1 = c("#676778","#D9D9C2"),
#'                        colors_condition2 = c("#0E6251", "#1B5E20"),
#'                        colors_condition3 = c("#5D3277","#AF6502"))

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
                                   top_n,
                                   cluster = TRUE,
                                   show_column_names = TRUE) {
  
  #Check for taxonomy column
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  
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
    dplyr::arrange(-abun) %>%  # Esto ordena por abundancia descendente
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
                                     phylum == "Actinobacteriota" ~ "Actinomycetota",
                                     phylum == "Bacteroidetes" ~ "Bacteroidota",
                                     TRUE~as.character(phylum)))
  warning("Note: Some bacterial phylum names have been updated to match NCBI's revised taxonomy:\n",
          "\nReference: https://ncbiinsights.ncbi.nlm.nih.gov/2021/12/10/ncbi-taxonomy-prokaryote-phyla-added/")

  
  ordered_taxa <- table_abundance$taxa
  
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
    dplyr::select(all_of(ordered_taxa)) %>%   
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(across(everything(), ~ as.numeric(.))) %>%  
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
  rownames(heatmap) <- ordered_taxa
  
  if (!cluster) {
    row_order <- ordered_taxa  
  } else {
    row_order <- NULL
  }
  

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
    annotation_name_gp = grid::gpar(fontsize = 12, fontface="bold", fontfamily= "serif"),
    gp = grid::gpar(col = "white"),
    col = cols_phyl
  )
  
  # Initialize lists for dynamic annotations and legends
  heatmap_annotations <- list()
  legend_list <- list(lgd1 = ComplexHeatmap::Legend(
    at = unique(annotation_rows$phylum), 
    legend_gp = grid::gpar(fill = c5.phylum), 
    title = "Phylum", 
    labels_gp = grid::gpar(fontsize=12, fontfamily= "serif")))
  
  
  # Process condition1 if exists
  if (!is.null(condition1)) {
    unique_vals <- sort(unique(annotation_columns[[condition1]]))
    
    if (is.null(colors_condition1)) {
      colors_condition1 <- viridis::viridis(length(unique_vals))
    } else if (length(colors_condition1) < length(unique_vals)) {
      colors_condition1 <- c(colors_condition1, 
                             viridis::viridis(length(unique_vals) - length(colors_condition1)))
    }
    
    color_mapping <- setNames(colors_condition1[1:length(unique_vals)], unique_vals)
    
    heatmap_annotations$ann1 <- ComplexHeatmap::HeatmapAnnotation(
      df = annotation_columns[condition1],
      which = "column",
      col = setNames(list(color_mapping), condition1),
      show_legend = FALSE,
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 12, fontface="bold", fontfamily= "serif"),
      gp = grid::gpar(col = "white")
    )
    
    legend_list$lgd2 <- ComplexHeatmap::Legend(
      at = names(color_mapping),
      legend_gp = grid::gpar(fill = color_mapping),
      title = if (!is.null(name_legend_condition1)) name_legend_condition1 else condition1,
      labels_gp = grid::gpar(fontsize=12, fontfamily= "serif")
    )
  }
#Process condition2 if exists
if (!is.null(condition2)) {
  unique_vals <- sort(unique(annotation_columns[[condition2]]))
  
  if (is.null(colors_condition2)) {
    colors_condition2 <- viridis::inferno(length(unique_vals))
  } else if (length(colors_condition2) < length(unique_vals)) {
    colors_condition2 <- c(colors_condition2, 
                           viridis::inferno(length(unique_vals) - length(colors_condition2)))
  }
  
  color_mapping <- setNames(colors_condition2[1:length(unique_vals)], unique_vals)
  
  heatmap_annotations$ann2 <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_columns[condition2],
    which = "column",
    col = setNames(list(color_mapping), condition2),
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 12, fontface="bold", fontfamily= "serif"),
    gp = grid::gpar(col = "white")
  )
  
  legend_list$lgd3 <- ComplexHeatmap::Legend(
    at = names(color_mapping),
    legend_gp = grid::gpar(fill = color_mapping),
    title = if (!is.null(name_legend_condition2)) name_legend_condition2 else condition2,
    labels_gp = grid::gpar(fontsize=12, fontfamily= "serif")
  )
}

# Process condition3 if exists
if (!is.null(condition3)) {
  unique_vals <- sort(unique(annotation_columns[[condition3]]))
  
  if (is.null(colors_condition3)) {
    colors_condition3 <- viridis::plasma(length(unique_vals))
  } else if (length(colors_condition3) < length(unique_vals)) {
    colors_condition3 <- c(colors_condition3, 
                           viridis::plasma(length(unique_vals) - length(colors_condition3)))
  }
  
  color_mapping <- setNames(colors_condition3[1:length(unique_vals)], unique_vals)
  
  heatmap_annotations$ann3 <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_columns[condition3],
    which = "column",
    col = setNames(list(color_mapping), condition3),
    show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 12, fontface="bold", fontfamily= "serif"),
    gp = grid::gpar(col = "white")
  )
  
  legend_list$lgd4 <- ComplexHeatmap::Legend(
    at = names(color_mapping),
    legend_gp = grid::gpar(fill = color_mapping),
    title = if (!is.null(name_legend_condition3)) name_legend_condition3 else condition3,
    labels_gp = grid::gpar(fontsize=12, fontfamily= "serif")
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
      labels_gp = grid::gpar(fontsize = 8.4, fontfamily= "serif"),
      legend_gp = grid::gpar(fontsize = 10.4, fontfamily= "serif"),
      title = "Relative abundance (%)",
      title_position = "topcenter",
      at = c(0,1,2,3,5,8,10,25, 50, 100),
      break_dist = 1
    ),
    rect_gp = grid::gpar(col = "black", lwd = 0.5),    
    row_names_gp = grid::gpar(fontsize = 11, fontface = "italic", fontfamily= "serif"), 
    column_names_gp = grid::gpar(fontsize=11, fontfamily= "serif"),
    cluster_columns = FALSE,
    cluster_rows = cluster,
    row_order = row_order,  # Añadido para mantener orden cuando cluster=FALSE
    show_column_names = show_column_names,
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
