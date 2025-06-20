#' Genera un diagrama Sankey a partir de una tabla OTU con información taxonómica
#'
#' @param table Data frame con una columna de taxonomía  y las demás columnas corresponden a muestras.
#' @param output_file Nombre del archivo HTML de salida (default: "sankey.html").
#' @param maxn Número máximo de taxones por nivel a incluir en el diagrama (default: 25).
#' @param taxRanks Niveles taxonómicos a visualizar (default: c("D","K","P","C","O","F","G","S")).
#' @param taxonomy_db Base de datos taxonómica, ej: "silva" o "kraken2" (default: "silva").

#' Genera un diagrama Sankey a partir de una tabla OTU con información taxonómica
#'
#' @param table Data frame con una columna de taxonomía y las demás columnas corresponden a muestras.
#' @param output_file Nombre del archivo HTML de salida (default: "sankey.html").
#' @param maxn Número máximo de taxones por nivel a incluir en el diagrama (default: 25).
#' @param taxRanks Niveles taxonómicos a visualizar (default: c("D","K","P","C","O","F","G","S")).
#' @param taxonomy_db Base de datos taxonómica, ej: "silva" o "kraken2" (default: "silva").

generate_sankey <- function(table, output_file = "sankey.html", maxn = 25,
                            taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S"),
                            taxonomy_db = "silva") {
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(sankeyD3)
  
  relabunda <- function(x) {
    as.data.frame(t(t(x) / colSums(x))) * 100
  }
  
  # Identificar columna de taxonomía
  tax_col <- grep("taxonomy|Taxonomy", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("No se encontró una única columna de taxonomía en table")
  
  # Reorganizar tabla: taxonomía al final
  table <- table[, c(setdiff(1:ncol(table), tax_col), tax_col)]
  
  # Matriz OTU sin taxonomía
  otu_mat <- table[, -ncol(table)]
  rownames(otu_mat) <- rownames(table)
  
  # Calcular abundancias relativas
  otu_rel <- relabunda(otu_mat)
  otu_rel$taxonomy <- table[[ncol(table)]]
  
  # Separar niveles taxonómicos
  otu_rel_parse <- otu_rel %>%
    rownames_to_column(var = "Feature.ID") %>%
    separate(taxonomy, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";", fill = "right") %>%
    mutate(across(where(is.character), ~ str_extract(., "[^_]+$")))
  
  # Modificación para Kraken2
  if (tolower(taxonomy_db) == "kraken2") {
    otu_rel_parse <- otu_rel_parse %>%
      mutate(s = ifelse(!is.na(g) & !is.na(s) & s != "NA", paste(g, s, sep = " "), s))
  }
  
  get_level_data <- function(df, level, unite_cols) {
    df %>%
      unite(col = !!level, all_of(unite_cols), remove = FALSE) %>%
      select(all_of(level), where(is.numeric)) %>%
      group_by(across(all_of(level))) %>%
      summarise(across(where(is.numeric), sum), .groups = "drop") %>%
      column_to_rownames(var = level) %>%
      t() %>%
      as.data.frame() %>%
      colMeans() %>%
      as.data.frame() %>%
      rename(abund = ".") %>%
      mutate(taxRank = toupper(substr(level, 1, 1))) %>%
      rownames_to_column(var = "Taxon") %>%
      filter(!endsWith(Taxon, "NA")) %>%
      column_to_rownames("Taxon")
  }
  
  bacterias <- get_level_data(otu_rel_parse, "kingdom", "k")
  phylum    <- get_level_data(otu_rel_parse, "phylum", c("k", "p"))
  class     <- get_level_data(otu_rel_parse, "class", c("k", "p", "c"))
  order     <- get_level_data(otu_rel_parse, "order", c("k", "p", "c", "o"))
  family    <- get_level_data(otu_rel_parse, "family", c("k", "p", "c", "o", "f"))
  genus     <- get_level_data(otu_rel_parse, "genus", c("k", "p", "c", "o", "f", "g"))
  specie    <- get_level_data(otu_rel_parse, "specie", c("k", "p", "c", "o", "f", "g", "s"))
  
  my_report <- bind_rows(bacterias, phylum, class, order, family, genus, specie) %>%
    rownames_to_column(var = "ids") %>%
    mutate(name = str_extract(ids, "[^_]+$")) %>%
    filter(!name %in% c("NA", "Bacteria", "Bacteria_Patescibacteria")) %>%
    filter(taxRank %in% taxRanks) %>%
    group_by(taxRank) %>%
    slice_max(order_by = abund, n = maxn, with_ties = FALSE) %>%
    ungroup()
  
  # Separar nodos en niveles
  splits <- strsplit(my_report$ids, "_")
  sel <- sapply(splits, length) >= 3
  splits <- splits[sel]
  
  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n = 2))),
                      stringsAsFactors = FALSE)
  colnames(links) <- c("source", "target")
  
  links$value <- sapply(seq_len(nrow(links)), function(i) {
    val <- my_report$abund[my_report$name == links$target[i]]
    if(length(val) == 1) val else NA_real_
  })
  
  links <- links[!is.na(links$value) & links$value > 0, ]
  
  # Corregir asignación de profundidad y etiquetas del eje
  valid_ranks <- taxRanks[taxRanks %in% unique(my_report$taxRank)]
  taxRank_to_depth <- setNames(seq_along(valid_ranks) - 1, valid_ranks)
  
  nodes <- data.frame(
    name = my_report$name,
    depth = taxRank_to_depth[my_report$taxRank],
    value = my_report$abund,
    stringsAsFactors = FALSE
  )
  
  names_id <- setNames(seq_len(nrow(nodes)) - 1, nodes$name)
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  
  links <- links[!is.na(links$source) & !is.na(links$target), ]
  links <- links[links$source != links$target, ]
  
  # Limpiar nombres
  nodes$name <- sub("^._", "", nodes$name)
  links$type <- sub(" .*", "", nodes[links$source + 1, "name"])
  
  sankey <- sankeyD3::sankeyNetwork(
    Links = links,
    Nodes = nodes,
    doubleclickTogglesChildren = TRUE,
    LinkGroup = "type",
    fontFamily = "Helvetica",
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    NodeGroup = "name",
    NodePosX = "depth",
    NodeValue = "value",
    dragY = TRUE,
    xAxisDomain = valid_ranks,  # Etiquetas de eje bien alineadas
    numberFormat = "pavian",
    title = NULL,
    nodeWidth = 15,
    linkGradient = TRUE,
    nodeShadow = TRUE,
    nodeCornerRadius = 5,
    units = "abund",
    fontSize = 10,
    iterations = 1000,
    align = "none",
    highlightChildLinks = TRUE,
    orderByPath = TRUE,
    scaleNodeBreadthsByString = TRUE
  )
  
  sankeyD3::saveNetwork(sankey, file = output_file)
  message("Sankey diagram saved to: ", output_file)
}
