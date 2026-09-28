devtools::load_all()
library(rbiom)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

###IMPORTAR LOS DATOS ####
bacteria_table <- read.delim(
  system.file("extdata", "table_bacteria.txt", package = "MicroBioMeta"),
  row.names = 1, check.names = FALSE
)

taxonomy_bacteria <- read.delim(
  system.file("extdata", "taxonomy_bacteria.txt", package = "MicroBioMeta"),
  check.names = FALSE
) %>%
  rename(taxonomy = Taxon) %>%
  dplyr::select(-Confidence) %>%
  column_to_rownames(var = "Feature.ID")

bacteria_metadata <- read.delim(
  system.file("extdata", "metadata_bacterias.txt", package = "MicroBioMeta"),
  check.names = FALSE
) %>%
  filter(Month == "2")

samples_bac <- bacteria_metadata$SAMPLEID[
  bacteria_metadata$SAMPLEID %in% colnames(bacteria_table)]

table_bacteria <- bacteria_table[, samples_bac]

metadata_bacteria <- bacteria_metadata[
  bacteria_metadata$SAMPLEID %in% samples_bac, ]

table_bac <- merge_feature_taxonomy(table = table_bacteria,
                                    taxonomy = taxonomy_bacteria)



###CREAR ARCHIVO PHYLOSEQ ####

# Matriz de abundancias
otu_bac <- otu_table(
  as.matrix(table_bacteria),
  taxa_are_rows = TRUE
)

# Taxonomía
tax_bac <- tax_table(
  as.matrix(taxonomy_bacteria)
)

# Metadata
meta_bac <- sample_data(
  data.frame(
    metadata_bacteria,
    row.names = metadata_bacteria$SAMPLEID
  )
)

# Crear objeto phyloseq
ps_bacteria <- phyloseq(
  otu_bac,
  tax_bac,
  meta_bac
)

# Revisar el objeto
ps_bacteria

####HACER COMPATIBLE LA TAXONOMÍA ####

# 1. Convertimos la tax_table actual en un data.frame para poder manipularla
tax_df <- as.data.frame(as(tax_table(ps_bacteria), "matrix"))

#limpia los prefijos como "d__", "p__", etc.
nombre_columna <- colnames(tax_df)[1]

tax_limpia <- tax_df %>%
  separate(!!sym(nombre_columna), 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";\\s*", # Separa por punto y coma (y espacios si los hay)
           fill = "right") %>% # Si falta la especie o género, rellena con NA
  mutate(across(everything(), ~ gsub("^[a-z]__", "", .))) # Quita los prefijos "p__", "c__", etc.

# 3. Convertimos las filas vacías o "NA"
tax_limpia[tax_limpia == ""] <- NA

#Reasignamos la nueva matriz taxonómica al objeto phyloseq
rownames(tax_limpia) <- rownames(tax_df) # Mantiene los IDs de las ASVs/OTUs
tax_table(ps_bacteria) <- tax_table(as.matrix(tax_limpia))



#### ABUNDANCIA RELATIVA ####
#1)MICROBIOMETA
system.time(abundance_bar_plot(table = table_bac,
                   metadata = metadata_bacteria,
                   taxonomy_db = "silva",
                   level = "phylum",
                   top_n = 15,
                   x_col = "Type_of_soil",
                   facet_col = "Treatment",
                   add_remained = TRUE,
                   label =   "Phylum",
                   x_label_angle = 45,
                   save_table = FALSE))

#2) PHYLOSEQ
## lo hice así para que diera un gráfico relativamente similar al
#que generamos con microbiometa

system.time({
  # 1. Aglomerar a nivel de Filo
  ps_filo <- tax_glom(ps_bacteria, taxrank = "Phylum")
  
  # 2. Identificar los 15 filos más abundantes
  top15_filos <- names(sort(taxa_sums(ps_filo), decreasing = TRUE)[1:15])
  
  # 3. Agrupar los filos restantes como "Otros"
  tax_mat <- as(tax_table(ps_filo), "matrix")
  tax_mat[!(rownames(tax_mat) %in% top15_filos), "Phylum"] <- "Otros"
  tax_table(ps_filo) <- tax_table(tax_mat)
  
  # Fusionar réplicas por grupo
  sample_data(ps_filo)$Group <- paste(sample_data(ps_filo)$Type_of_soil, sample_data(ps_filo)$Treatment, sep = "_")
  ps_merged <- merge_samples(ps_filo, "Group")
  
  # Reconstruir metadatos
  df_meta <- data.frame(Group = sample_names(ps_merged))
  df_meta <- separate(df_meta, Group, into = c("Type_of_soil", "Treatment"), sep = "_", remove = FALSE)
  rownames(df_meta) <- df_meta$Group
  sample_data(ps_merged) <- sample_data(df_meta)
  
  # 4. Abundancia relativa
  ps_rel <- transform_sample_counts(ps_merged, function(x) (x / sum(x)))
  
  # 5. Generar el gráfico
  p <- plot_bar(ps_rel, x = "Type_of_soil", fill = "Phylum", facet_grid = ~Treatment)
  grafico_final <- p + geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")
  p})




## 3)RBIOM
#Convertir el objeto phyloseq a Rbiom
data <- as_rbiom(ps_bacteria)

#data2<- rarefy(data)

system.time({
  data3 <- biom_relativize(data); 
  taxa_stacked(data3, rank = 'Phylum', taxa= 10, label.by = "Type_of_soil", facet.by = "Treatment")
})



#### HEATMAPS ####
#1) MICROBIOMETA
system.time(abundance_heatmap_plot(table = table_bac,
                                   metadata = metadata_bacteria,
                                   top_n = 20,
                                   show_column_names = FALSE,
                                   condition1 = "Treatment",
                                   condition2 = "Type_of_soil",
                                   feature_prefix = "ASV"))

#2)PHYLOSEQ
# usa prune_taxa para eliminar los grupos con 0 porque si no marca
#un error

system.time({
  ps_bacteria2 <- prune_taxa(
    taxa_sums(ps_bacteria) > 0,
    ps_bacteria
  )
  
  ps_bacteria3 <- tax_glom(
    ps_bacteria2,
    taxrank = "Phylum"
  )
  
  plot_heatmap(
    ps_bacteria3,
    taxa.label = "Phylum",
    sample.label = "Type_of_soil"
  )
})


###rbiom

system.time(taxa_heatmap(
  data,
  taxa = 20,
  tracks = c('Treatment', 'Type_of_soil'),
  grid= blues9))


####DIVERSIDAD ALFA ####

#1)MICROBIOMETA
system.time(alpha_diversity_plot(table = table_bac,
                            metadata = metadata_bacteria,
                            x_col = "Treatment",
                            fill_col = "Treatment",
                            facet_by = "Type_of_soil",
                            facet_orientation = "horizontal",
                            save_table = FALSE))

#2)PHYLOSEQ
system.time({
  p <- plot_richness(
    ps_bacteria,
    x = "Type_of_soil",
    color = "Treatment",
    measures = c("Chao1", "Shannon", "Simpson")
  )
  
  p <- p + facet_grid(Type_of_soil ~ variable)
  p
})

#3)Rbiom
system.time(adiv_boxplot(
  data,
  x= "Treatment",
  facet.by = "Type_of_soil",
  adiv= c('otu', 'shannon', "simpson")))



####ORDINATION BDIV ####
#1) MICROBIOMETA
system.time(beta_ord_plot(table = table_bac,
              metadata = metadata_bacteria,
              distance = "bray",
              ordination = "PCoA",
              group_col = "Treatment",
              shape_col = "Type_of_soil"))

#2)PHYLOSEQ
system.time({GP.ord <- ordinate(ps_bacteria, "PCoA", "bray")
p1 = plot_ordination(ps_bacteria, GP.ord, type="SAMPLEID", color="Treatment",
                     title="PCoA", shape = "Type_of_soil")
print(p1)})
#Se hace en dos pasos separados, primero ordinate y luego el plot

#3)RBIOM
system.time(bdiv_ord_plot(
  data,
  bdiv = "bray",
  ord = "PCoA",
  layers = "petm",
  stat.by = 'Treatment',
  facet.by =  NULL,
  colors = TRUE,
  shapes = TRUE,
  tree = NULL,
  test = "adonis2",
  seed = 0,
  permutations = 999))

##### BDIV ESTADÍSTICA####
#1)MICROBIOMETA
system.time(beta_test_table(table = table_bac,
                metadata= metadata_bacteria,
                formula_str = "Treatment*Type_of_soil",
                method = "compositional",
                test = "permanova",
                permutations = 999))

#2)PHYLOSEQ
# phyloseq no te permite hacer los estadísticos

# 3)Rbiom
## usa estadísticos univariados
system.time(stats<- bdiv_stats(data,
                   stat.by = "Treatment"))

# 1. Genera la matriz con la clase "dist" correcta
system.time({dm <- bdiv_distmat(data, bdiv = "bray")

# 2. Extrae tu vector con nombres
sex_vector <- pull(data, "Treatment")

# 3. Ahora el PERMANOVA correrá sin inconvenientes
distmat_stats(dm, groups = sex_vector, test = "adonis2")})

#### CCA Y RDA ####

env_table_bac <- metadata_bacteria %>%
  dplyr::select(SAMPLEID, pH:Arbus_per) %>%
  remove_rownames() %>%
  column_to_rownames(var = "SAMPLEID")

#1)MICROBIOMETA
system.time(cca_rda_biplot(table = table_bac,
               env_data = env_table_bac,
               metadata = metadata_bacteria,
               env_vars = c("pH","TN","WHC", "EC", "Clay"),
               analysis = "CCA",
               show_all_env_vectors = TRUE,
               group_col = "Type_of_soil",
               scale_arrows = 3))

#2)PHYLOSEQ

system.time({cca_resultado <- ordinate(ps_bacteria, method = "CCA", formula = ~ pH+ TN+ WHC+ EC+ Clay)

p<- plot_ordination(ps_bacteria, cca_resultado, type = "SAMPLEID", color = "Type_of_soil") +
  theme_minimal() +
  geom_point(size = 3) +
  labs(title = "Análisis de Correspondencia Canónica (CCA)")

# 2. Extraer las coordenadas de las variables explicativas (vectores)
flechas <- as.data.frame(cca_resultado$CCA$biplot)
flechas$Variable <- rownames(flechas)

# 3. Superponer las flechas en el gráfico
p + geom_segment(data = flechas, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue", inherit.aes = FALSE) +
  geom_text(data = flechas, aes(x = CCA1, y = CCA2, label = Variable), 
            color = "blue", vjust = -0.5, inherit.aes = FALSE)
})
