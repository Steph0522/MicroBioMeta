# testing microbiometa with sceloporus data

#load data
library(qiime2R)
library(tidyverse)
devtools::load_all()

table <- read_qza("tests/data/run_f250_r230_feature-table.qza")$data %>% as.data.frame()
depth <- read.delim("tests/data/depth.csv", sep = ";") 
depth <- depth[match(colnames(table), depth$Sample.ID),]

metadata <- read.csv("tests/data/Metadata.csv", sep = ";") %>% 
  dplyr::rename(SAMPLEID = "INDEX") %>%
  mutate_all(as.character) %>% 
  dplyr::select(SAMPLEID, everything())

new_names <- sub(".*\\.", "", colnames(table))
colnames(table) <- new_names
depth$Sample.ID <- new_names
metadata_depth <- metadata %>%
  full_join(depth, by = c("SAMPLEID" = "Sample.ID")) %>%
  filter(!is.na(Frequency) & Frequency > 5000)



taxonomy_gg2<- read_qza("tests/data/run_f250_r230_taxa_gg2_scel.qza")$data %>%column_to_rownames(var = "Feature.ID") %>% dplyr::select(-Confidence)
taxonomy_silva<- read_qza("tests/data/run_f250_r230_taxa_silva.qza")$data %>%column_to_rownames(var = "Feature.ID") %>%dplyr::select(-Confidence)
taxonomy_gg2_weighted<- read_qza("tests/data/taxonomy_gg2_weighted.qza")$data %>%column_to_rownames(var = "Feature.ID") %>% dplyr::select(-Confidence)
taxonomy_silva_weighted<- read_qza("tests/data/taxonomy_silva_weighted.qza")$data %>%column_to_rownames(var = "Feature.ID") %>% dplyr::select(-Confidence)


parse_taxa_gg2 <- qiime2R::parse_taxonomy(read_qza("tests/data/run_f250_r230_taxa_gg2_scel.qza")$data)
parse_taxa_gg2_weighted <- qiime2R::parse_taxonomy(read_qza("tests/data/taxonomy_gg2_weighted.qza")$data)
parse_taxa_silva <- qiime2R::parse_taxonomy(read_qza("tests/data/run_f250_r230_taxa_silva.qza")$data)
parse_taxa_silva_weighted <- qiime2R::parse_taxonomy(read_qza("tests/data/taxonomy_silva_weighted.qza")$data)



table_taxa <- merge_feature_taxonomy(table, taxonomy_gg2)
table_taxa2 <- merge_feature_taxonomy(table, taxonomy_gg2_weighted)
table_taxa3 <- merge_feature_taxonomy(table, taxonomy_silva)
table_taxa4 <- merge_feature_taxonomy(table, taxonomy_silva_weighted)


table_taxa[,-55] <- sweep(table_taxa[,-55], 2, colSums(table_taxa[,-55], na.rm = TRUE), FUN = "/") * 100
table_taxa2[,-55] <- sweep(table_taxa2[,-55], 2, colSums(table_taxa2[,-55], na.rm = TRUE), FUN = "/") * 100
table_taxa3[,-55] <- sweep(table_taxa3[,-55], 2, colSums(table_taxa3[,-55], na.rm = TRUE), FUN = "/") * 100
table_taxa4[,-55] <- sweep(table_taxa4[,-55], 2, colSums(table_taxa4[,-55], na.rm = TRUE), FUN = "/") * 100

table_parse_taxa_gg2_genus <- table_taxa %>% 
  rownames_to_column(var = "Feature.ID") %>% 
  inner_join(parse_taxa_gg2 %>% rownames_to_column(var = "Feature.ID")) %>% 
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), mean), .groups = "drop") %>% 
  drop_na() %>% 
  mutate(sums = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(sums))

table_parse_taxa_gg2_w_genus <- table_taxa2 %>% 
  rownames_to_column(var = "Feature.ID") %>% 
  inner_join(parse_taxa_gg2_weighted %>% rownames_to_column(var = "Feature.ID")) %>% 
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), mean), .groups = "drop") %>% 
  drop_na() %>% 
  mutate(sums = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(sums))

table_parse_taxa_silva_genus <- table_taxa3 %>% 
  rownames_to_column(var = "Feature.ID") %>% 
  inner_join(parse_taxa_silva %>% rownames_to_column(var = "Feature.ID")) %>% 
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), mean), .groups = "drop") %>% 
  drop_na() %>% 
  mutate(sums = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(sums))

table_parse_taxa_silva_w_genus <- table_taxa4 %>% 
  rownames_to_column(var = "Feature.ID") %>% 
  inner_join(parse_taxa_silva_weighted %>% rownames_to_column(var = "Feature.ID")) %>% 
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), mean), .groups = "drop") %>% 
  drop_na() %>% 
  mutate(sums = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(sums))

dim(table_parse_taxa_gg2_genus)
dim(table_parse_taxa_gg2_w_genus)
dim(table_parse_taxa_silva_genus)
dim(table_parse_taxa_silva_w_genus)

table_taxa <- merge_feature_taxonomy(table, taxonomy_gg2)
table_taxa2 <- merge_feature_taxonomy(table, taxonomy_gg2_weighted)
table_taxa3 <- merge_feature_taxonomy(table, taxonomy_silva)
table_taxa4 <- merge_feature_taxonomy(table, taxonomy_silva_weighted)

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

abundance_barplot(
  table = table_taxa3,
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
  table = table_taxa4,
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

alpha_hill_corrplot(table = table_taxa2, 
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
  table_taxa3,
  metadata %>% drop_na(),
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
  table_taxa3,
  output_file <- file.path(getwd(), "sankey_scel_gg2.html"),
  maxn = 10,
  taxRanks = c("P", "C", "G", "S"),
  taxonomy_db = "silva"
)

#getwd()
#summarice

metadata %>%  group_by(SITIO, ID.CAM) %>% count()
devtools::load_all()




ancombc_plot(table = table_taxa, metadata = metadata,
             col_cond = "SITIO", p_adj_method = "BH", prv_cut = 0.08)
