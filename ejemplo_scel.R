# testing microbiometa with sceloporus data

#load data
library(qiime2R)
library(tidyverse)
devtools::load_all()

table <- read_qza("test_data/run_f250_r230_feature-table.qza")$data
metadata <- read.csv("test_data/Metadata.csv", sep = ";") %>% 
  dplyr::rename(SAMPLEID = "INDEX") %>%
  mutate_all(as.character) %>% 
  dplyr::select(SAMPLEID, everything())
new_names <- sub(".*\\.", "", colnames(table))
colnames(table) <- new_names
table<- as.data.frame(table)
taxonomy<- read_qza("test_data/run_f250_r230_taxa_gg2.qza")$data %>%
  column_to_rownames(var = "Feature.ID") %>%
  dplyr::select(-Confidence)

taxonomy2<- read_qza("test_data/run_f250_r230_taxa_silva.qza")$data %>%
  column_to_rownames(var = "Feature.ID") %>%
  dplyr::select(-Confidence)


table_taxa <- merge_feature_taxonomy(table, taxonomy)
table_taxa2 <- merge_feature_taxonomy(table, taxonomy2)


abundance_barplot(
  table = table_taxa,
  metadata = metadata,
  taxonomy_db = "silva",
  level = "genus",
  x_col = "MUESTRA",
  label = "Genus",
  facet_col = "SITIO",
  width_equal = FALSE,
  # group_var = "SAMPLEID",
  top_n_groups = 30,
  add_remained  = TRUE
)

abundance_barplot(
  table = table_taxa2,
  metadata = metadata,
  taxonomy_db = "silva",
  level = "genus",
  x_col = "MUESTRA",
  label = "Genus",
  facet_col = "SITIO",
  width_equal = FALSE,
  # group_var = "SAMPLEID",
  top_n_groups = 30,
  add_remained  = TRUE
)

abundance_heatmap_plot(
  table = table_taxa,
  metadata = metadata,
  condition1 = "SITIO",
  condition2 = "ID.CAM",
  top_n = 30,
  cluster = TRUE,
  colors_condition1 = c("red", "blue"),
  colors_condition2 = c("green", "orange"),
  show_column_names = FALSE
)

alpha_hill_corrplot(table = table_taxa, 
                    facet_orientation = "horizontal")

alpha_hill_plot(
  table = table_taxa,
  metadata = metadata %>% filter(!SAMPLEID=="30"),
  type = "boxplot",
  fill_col = "SITIO",
  x_col = "SITIO",
  facet_orientation = "horizontal",
  #facet_by = "edad",
  #     facet_by2 = "estructura",
  free_y = T,
  legend_position = "top",
  stat = "kruskal.test"
)

metadata1 = metadata %>% filter(!SAMPLEID=="30")
ordered = metadata1$SAMPLEID
table_taxa <- table_taxa %>% dplyr::select(ordered, taxonomy)

beta_div_plot(
  table = table_taxa,
  metadata = metadata1,
  distance = "aitchison",
  ordination = "NMDS",
  group_col  = "SITIO",
  n_taxa = 5,
  arrows = 100
)


randomf_lollipop_plot(
  table_taxa,
  metadata,
  variable_to_predict = "SITIO",
  col_palette = c("red", "blue", "green"),
  top_n = 20,
  size = 6
)



venn_diagram_plot(
  table_taxa,
  metadata,
  merge_by = "ID.CAM",
  min_prevalence = 0,
  #  denom = "all",
  group_colors = c("blue", "yellow"),
  method = "ggvenn"
)


abundance_sankey_plot(
  table_taxa2,
  output_file = "sankey_scel2.html",
  maxn = 20,
  taxonomy_db = "silva"
)


#summarice

metadata %>%  group_by(SITIO, ID.CAM) %>% count()
