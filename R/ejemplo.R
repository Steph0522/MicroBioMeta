
otu = read.delim("otutable_with_taxonomy.txt", skip = 1, row.names = 1) 
library(tidyverse)
table = otu %>% dplyr::select(taxonomy, everything())
metadata = read.delim("metadata.txt") %>% rename(SAMPLEID="sample.id")

alpha_hill_plot(table = table,
                metadata = metadata, 
                fill_col = "metodo", 
                x_col = "metodo", facet_orientation = "vertical",
                free_y = T,
                legend_position = "right",
                stat = "t.test")


relative_abundance_plot(table = table,
                        metadata = metadata,
                        taxonomy_db = "silva",
                        level = "phylum", 
                        x_col = "estructura",
                        label = "Phylum",
                        facet_col = "edad",
                        group_var = "SAMPLEID",
                        top_n_groups = 10)

effect_size_plot(
  table = table[-1],
  metadata= metadata,
  col_cond = "metodo",
  threshold_lower = -2,
  threshold_upper = 2,
  cond = kit,
  show_labels = TRUE
)
