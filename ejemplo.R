getwd()

otu = read.delim("otutable_with_taxonomy.txt", skip = 1, row.names = 1) 
library(tidyverse)
table = otu %>% dplyr::select(taxonomy, everything())
metadata = read.delim("meta.txt", sep = "")%>% rename(SAMPLEID="sample.id")
metadata = metadata[match(colnames(otu), metadata$SAMPLEID),] %>% filter(!SAMPLEID =="NA")

alpha_hill_corrplot(table = table[-1],
                    facet_orientation = "horizontal")

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
                       # facet_col = "edad",
                        group_var = "SAMPLEID",
                        top_n_groups = 10)

effect_size_plot(
  table = table[-1],
  metadata= metadata,
  col_inf = "blue",
  col_sup = "red",
  col_cond = "metodo",
  threshold_lower = -1,
  threshold_upper = 1,
  cond = "kit",
  show_labels = TRUE
)


