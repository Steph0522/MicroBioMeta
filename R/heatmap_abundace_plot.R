###FUNCTION####

heatmap_abundance_plot <- function(table,
                         metadata,
                         condition1,
                         condition2,
                         condition3,
                         colors_condition1,
                         colors_condition2,
                         colors_condition3,
                         name_legend_condition1,
                         name_legend_condition2,
                         name_legend_condition3,
                         arguments_condition1,
                         arguments_condition2,
                         arguments_condition3,
                         top_n) {
  
#Delete taxonomy
table_counts <- table %>%
    
dplyr::select(-taxonomy) 
  
#Separate taxonomy
table_tax <- table %>%
tibble::rownames_to_column("OTUID") %>%
dplyr::select(OTUID, taxonomy)
  
#Initial table
phy.ra.complete <- t(t(table_counts)/colSums(table_counts)*100) %>%
as.data.frame()
  
#Editing table
table_abundance <- phy.ra.complete %>%
mutate(abun = rowMeans(.)) %>% 
tibble::rownames_to_column(var="OTUID") %>% 
dplyr::arrange(-abun) %>%
dplyr::slice(1:top_n) %>%
dplyr::left_join(table_tax, by = "OTUID")  %>%
mutate(taxonomy2 = taxonomy) %>%
tidyr::separate(taxonomy2, into = c("dominio","phylum","clase","orden","familia","genero","especie"), sep = ";", fill = "right", extra = "merge") %>% 
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
                                     phylum=="Firmicutes" ~ "Bacillota",  TRUE~as.character(phylum)))
  
  
#Join table with metadata
heat <- table_abundance %>% 
tibble::column_to_rownames(var = "taxa") %>%
t() %>%
as.data.frame() %>%
tibble::rownames_to_column(var = "OTUID") %>%
dplyr::inner_join(metadata, by="OTUID") 
  
#Transformation of ranges
heatmap <- heat %>%
tibble::column_to_rownames(var = "OTUID") %>%
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
  
#annotation row 
annotation_rows <- table_abundance %>%
dplyr::select(OTUID, phylum) %>%
tibble::column_to_rownames(var = "OTUID")
rownames(annotation_rows) <- rownames(heatmap)
  
#annotations columns
annotation_columns <- heat %>%
dplyr::select(all_of(c(condition1, condition2, condition3))) 
rownames(annotation_columns) <- colnames(heatmap)
  
#Color palette
#Color palette to heatmap
my_palette <- viridis::viridis (n = 12, option = "C", direction = -1)
#Color for annotations
c5.phylum <- RColorBrewer::brewer.pal(8, "Set2")[1:6]
cols_phyl <- list(Phylum = setNames(c5.phylum, unique(annotation_rows$phylum)))
  
#arguments for annotations
cols_condition1 <- setNames(list(arguments_condition1), condition1)
cols_condition2 <- setNames(list(arguments_condition2), condition2)
cols_condition3 <- setNames(list(arguments_condition3), condition3)
  
#Annotations
annphylum = ComplexHeatmap::rowAnnotation(Phylum = annotation_rows$phylum, 
                                            show_legend = F,
                                            show_annotation_name = T,
                                            annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
                                            gp = grid::gpar(col = "white"),
                                            col = cols_phyl)
  
  
anncondition1 = ComplexHeatmap::HeatmapAnnotation(df = annotation_columns[condition1],
                                                    which = "column", 
                                                    col = cols_condition1,
                                                    show_legend = FALSE,
                                                    show_annotation_name = TRUE,
                                                    annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
                                                    gp = grid::gpar(col = "white"))
  
  
anncondition2 = ComplexHeatmap::HeatmapAnnotation(df = annotation_columns[condition2],
                                                    which = "column", 
                                                    col = cols_condition2,
                                                    show_legend = FALSE, 
                                                    show_annotation_name = TRUE,
                                                    annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
                                                    gp = grid::gpar(col = "white"))
  
  
  
anncondition3 = ComplexHeatmap::HeatmapAnnotation(df = annotation_columns[condition3],
                                                    which = "column", 
                                                    col = cols_condition3,
                                                    show_legend = FALSE, 
                                                    show_annotation_name = TRUE,
                                                    annotation_name_gp = grid::gpar(fontsize = 11, fontface="bold"),
                                                    gp = grid::gpar(col = "white"))
  
#Legends
lgd1 = ComplexHeatmap::Legend(at = unique(annotation_rows$phylum), legend_gp = grid::gpar(fill = c5.phylum), title = "Phylum", labels_gp = grid::gpar(fontsize=11))
lgd2 = ComplexHeatmap::Legend(at = sort(unique(annotation_columns[[condition1]])), legend_gp = grid::gpar(fill = colors_condition1), title = name_legend_condition1, labels_gp = grid::gpar(fontsize=11))
lgd3 = ComplexHeatmap::Legend(at = sort(unique(annotation_columns[[condition2]])), legend_gp = grid::gpar(fill = colors_condition2), title = name_legend_condition2, labels_gp = grid::gpar(fontsize=11))
lgd4 = ComplexHeatmap::Legend(at = sort(unique(annotation_columns[[condition3]])), legend_gp = grid::gpar(fill = colors_condition3), title = name_legend_condition3, labels_gp = grid::gpar(fontsize=11))
  
pd.legends <- ComplexHeatmap::packLegend(lgd1, lgd2, lgd3, lgd4)
  
  
#######HEATMAP#####
heats <- ComplexHeatmap::Heatmap(heatmap, col=my_palette,
                                   heatmap_legend_param = list(direction = "horizontal",
                                                               labels_gp = grid::gpar(fontsize = 7),
                                                               legend_gp = grid::gpar(fontsize = 9),
                                                               title = "Relative abundance (%)",
                                                               title_position = "topcenter",
                                                               at = c(0,1,2,3,5,8,10,25, 50, 100),
                                                               break_dist = 1),
                                   rect_gp = grid::gpar(col = "black", lwd = 0.5),    
                                   row_names_gp =  grid::gpar(fontsize=7),
                                   column_names_gp = grid::gpar(fontsize=7),
                                   cluster_columns = FALSE,
                                   cluster_rows = TRUE,
                                   show_column_names = TRUE,
                                   show_heatmap_legend = TRUE, 
                                   top_annotation = anncondition1 + anncondition2 + anncondition3, 
                                   left_annotation = annphylum)
  
  
heatmap_output <- grid.grabExpr(draw(heats,
                                       heatmap_legend_side = "top",
                                       annotation_legend_side = "right", 
                                       merge_legend=FALSE,
                                       annotation_legend_list = pd.legends))
  
  return(heatmap_output)
  
}
