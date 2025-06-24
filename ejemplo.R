library(tidyverse)
getwd()

otu = read.delim("otutable_with_taxonomy.txt",
                 skip = 1,
                 row.names = 1) %>%
  dplyr::select(-taxonomy)

taxonomy = read.delim("taxonomy.tsv", row.names = 1) %>%
  dplyr::select(-Confidence)

table = merge_feature_taxonomy(otu, taxonomy)


library(tidyverse)
metadata = read.delim("meta.txt", sep = "") %>% rename(SAMPLEID = "sample.id")
metadata = metadata[match(colnames(otu), metadata$SAMPLEID), ] %>% filter(!SAMPLEID ==
                                                                            "NA")
devtools::load_all()



abundance_heatmap_plot(
  table = table,
  metadata = metadata,
  condition1 = "metodo",
  condition2 = "edad",
  condition3 = "estructura.metodo",
  top_n = 20,
  cluster = TRUE,
  colors_condition1 = c("red", "blue"),
  colors_condition2 = c("green", "orange"),
  show_column_names = FALSE
)

aldex_heatmap_plot(
  table = table,
  metadata = metadata,
  col_cond = "metodo",
  heatmap_colors = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
  treatment_colors = c("Higher" = "green", "Lower" = "yellow")
)

aldex_volcano_plot(
  table = table,
  metadata = metadata,
  col_inf = "blue",
  col_sup = "red",
  col_cond = "metodo",
  type = "effect",
  threshold_lower = -1,
  threshold_upper = 1,
  cond = "kit",
  show_labels = TRUE
)



alpha_hill_corrplot(table = table, facet_orientation = "horizontal")

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

beta_div_plot(
  table = table,
  metadata = metadata,
  distance = "compositional",
  ordination = "PCA",
  color_by = "metodo",
  n_taxa = 5,
  shape_by = "edad"
)

# let's make some example data for the next analysis

set.seed(123) # For reproducibility

# Create environmental data frame
env_data <- data.frame(
  SAMPLEID = c(
    "Fe1",
    "Fe2",
    "Fe3",
    "Fe4",
    "Fe5",
    "Fe6",
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "Q5",
    "Q6"
  ),
  pH = round(runif(12, 5.5, 7.5), 1),
  # pH values between 5.5-7.5
  Conductivity = round(rnorm(12, mean = 500, sd = 150)),
  # µS/cm
  Organic_Matter = round(rnorm(12, mean = 3, sd = 0.8), 1),
  # %
  Nitrogen = round(rnorm(12, mean = 0.15, sd = 0.05), 2),
  # %
  Phosphorus = round(rnorm(12, mean = 25, sd = 8)),
  # mg/kg
  Potassium = round(rnorm(12, mean = 150, sd = 40)),
  # mg/kg
  Calcium = round(rnorm(12, mean = 500, sd = 200)),
  # mg/kg
  row.names = "SAMPLEID"
)

# View the generated data
print(env_data)

cca_rda_biplot(
  table = table,
  env_data = env_data,
  metadata = metadata ,
  group_col = "metodo",
  analysis = "RDA",
  show_all_env_vectors = TRUE,
  legend_title = "Método",
  group_colors = c("red", "blue"),
  env_vars = c("pH", "Nitrogen", "Calcium"),
  title = "tittle"
  
  
)

randomf_lollipop_plot(
  table,
  metadata,
  variable_to_predict = "metodo",
  col_palette = c("red", "blue", "green"),
  top_n = 20,
  size = 6
)



relative_abundance_plot(
  table = table,
  metadata = metadata,
  taxonomy_db = "silva",
  #checar con unite, gg2, gtdb, kraken
  level = "genus",
  x_col = "metodo",
  label = "Genus",
  # facet_col = "edad",
  # group_var = "SAMPLEID",
  top_n_groups = 10,
  add_remained  = TRUE,
)


abundance_sankey_plot(
  table,
  output_file = "sankey_myxo.html",
  maxn = 20,
  taxonomy_db = "silva"
)


venn_diagram_plot(
  table,
  #%>% remove_rownames(),
  metadata,
  #%>% remove_rownames(),
  merge_by = "metodo",
  min_prevalence = 0,
  #  denom = "all",
  group_colors = c("blue", "yellow"),
  method = "ggvenn"
)

