table<- read.delim("table_species_plus_composta.tsv", row.names = 1)
table<- table[,-9]

columnas<-colnames(table)
columnas<- gsub("\\.","_", columnas)
columnas<- gsub("\\_RO","_Ra", columnas)
columnas<- gsub("\\_Ra1","_Ra", columnas)
columnas<- gsub("X","", columnas)
colnames(table)<- columnas

# Filtrar solo filas que pertenezcan al reino Bacteria
table <- table %>%
  filter(str_detect(taxonomy, "^k__Bacteria"))

# Definir niveles comunes para TRATAMIENTO_SAMPLE_DAY
tratamiento_day_levels <- c( "Soil_m_Soil_Jul","Soil_Soil_Jul", "No.compost_Soil_Ago","Compost_Soil_Ago",
                             "Sterile.compost_Soil_Ago",  
                             "No.compost_Rhizosphere_Ago", 
                             "Compost_Rhizosphere_Ago", 
                             "Sterile.compost_Rhizosphere_Ago", 
                             "No.compost_Root_Ago",  "Compost_Root_Ago",
                             "Sterile.compost_Root_Ago")

tratamiento_day_labels <- c("CM","S", "Un_S", "C_S", "SC_S", 
                            "Un_R","C_R", "SC_R", 
                            "Un_Ro",  "C_Ro",  "SC_Ro")

metadata <- read.delim("sample_metadata_TODO.tsv", row.names = NULL) %>% 
  slice(1:33) %>%          # Tomar las primeras 33 filas
  slice(-10) %>% 
  slice(1:32) %>%            
  mutate(
    TRATAMIENTO_SAMPLE_DAY = paste(Treatment, Sample_type, Time.days, sep = "_"),
    TRATAMIENTO_SAMPLE_DAY = factor(TRATAMIENTO_SAMPLE_DAY, levels = tratamiento_day_levels, labels = tratamiento_day_labels)
  ) %>%
  arrange(Sample_type)

metadata$Sample_type = factor(metadata$Sample_type, levels = c("Soil", "Rhizosphere", "Root"), labels = c("Soil", "Rhizosphere", "Root"))


alpha_hill_plot(
  table = table,
  metadata = metadata,
  x_col = "Treatment",
 # stat = "kruskal.test",
  fill_col = "Treatment",
  facet_by = "Sample_type",
#  facet_by2 = "Time.days",
 # fill_palette = "grey",
  facet_orientation = "horizontal",
  legend_position = "right",
  free_y=FALSE
)+theme(axis.text.x = element_blank(), axis.ticks = element_blank())


randomf_lollipop_plot(
  table,
  metadata,
  variable_to_predict = "Sample_type",
  col_palette = c("red", "blue", "green"),
  top_n = 20, size = 6
)


metadata_nosoil = metadata %>% filter(!Sample_type =="Soil")
taxonomy_col <- table$taxonomy
table_filtered <- table[, match(metadata_nosoil$SAMPLEID, colnames(table))]
table_nosoil <- cbind(table_filtered, taxonomy = taxonomy_col)

aldex_anal = aldex(reads = table_nosoil %>% dplyr::select(-taxonomy), conditions = metadata_nosoil$Sample_type)

aldex_heatmap_plot(table = table_nosoil,
                   metadata= metadata_nosoil,
                   col_cond= "Sample_type", effect_threshold = 2)
