getwd()

otu = read.delim("otutable_with_taxonomy.txt",
                 skip = 1,
                 row.names = 1)
library(tidyverse)
table = otu #%>% dplyr::select(taxonomy, everything())
metadata = read.delim("meta.txt", sep = "") %>% rename(SAMPLEID = "sample.id")
metadata = metadata[match(colnames(otu), metadata$SAMPLEID), ] %>% filter(!SAMPLEID ==
                                                                            "NA")
devtools::load_all()

aldex_heatmap_plot(table = table, 
                   metadata = metadata,
                   col_cond = "metodo",
                   heatmap_colors = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
                   treatment_colors = c("Higher" = "green", "Lower" = "yellow"))

aldex_volcano_plot(
  table = table,
  metadata = metadata,
  col_inf = "blue",
  col_sup = "red",
  col_cond = "metodo",
  type ="effect",
  threshold_lower = -1,
  threshold_upper = 1,
  cond = "kit",
  show_labels = TRUE
)



alpha_hill_corrplot(table = table[-1],
                    facet_orientation = "horizontal")

alpha_hill_plot(
  table = table,
  metadata = metadata,
  type = "barplot",
  fill_col = "metodo",
  x_col = "metodo",
  facet_orientation = "horizontal",
  #facet_by = "edad",
  #     facet_by2 = "estructura",
  free_y = T,
  legend_position = "top",
  stat = "t.test"
)

alpha_diversity_plot(
  table = table %>% remove_rownames(),
  metadata = metadata,
  type = "barplot",
  fill_col = "metodo",
  x_col = "metodo",
  facet_orientation = "horizontal",
  #   facet_by = "edad",
  #  facet_by2 = "estructura",
  free_y = T,
  legend_position = "top",
  stat = "t.test"
)


relative_abundance_plot(
  table = table,
  metadata = metadata,
  taxonomy_db = "silva",
  #checar con unite, gg2, gtdb, kraken
  level = "genus",
  x_col = "metodo",
  label = "Phylum",
  # facet_col = "edad",
  # group_var = "SAMPLEID",
  top_n_groups = 10,
  add_remained  = TRUE,
)


venn_diagram(
  table = table %>% tibble::remove_rownames(),
  metadata %>% remove_rownames(),
  merge_by = "metodo",
  method = "microeco"
)


venn_diagram_plot(
  table,
  #%>% remove_rownames(),
  metadata,
  #%>% remove_rownames(),
  merge_by = "metodo",
  min_prevalence = 0,
  #  denom = "all",
  #  group_colors = c("blue", "yellow"),
  method = "ggvenn"
)



#metadata$metodo <- as.factor(metadata$metodo)

randomf_lollipop_plot(table,
                      metadata,
                      variable_to_predict = "metodo",
                      top_n = 20)


beta_div_plot(
  table = table,
  metadata = metadata,
  distance = "compositional",
  ordination = "PCA",
  color_by = "metodo",
  n_taxa = 5
)

#data
data("varechem")
data("varespec")

#extracting species and env-data
species_data <- varespec[, 1:18]
env_data <- varechem[, 2:7]
metadata = data.frame(ids = 1:24, group = c(rep("A", 12), rep("B", 12)))

cca_biplot(
  table = t(species_data),
  env_data = env_data,
  metadata = metadata,
  group_col = "group"
)
