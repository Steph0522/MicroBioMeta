#' Genera un diagrama Sankey a partir de una tabla OTU con información taxonómica
#'
#' @param table Data frame con una columna de taxonomía  y las demás columnas corresponden a muestras.
#' @param output_file Nombre del archivo HTML de salida (default: "sankey.html").
#' @param maxn Número máximo de taxones por nivel a incluir en el diagrama (default: 25).
#' @param taxRanks Niveles taxonómicos a visualizar (default: c("D","K","P","C","O","F","G","S")).
#' @param taxonomy_db Base de datos taxonómica, ej: "silva" o "kraken2" (default: "silva").


generate_sankey <- function(table, output_file = "sankey.html", maxn = 25,
                            taxRanks = c("D","K","P","C","O","F","G","S"),
                            taxonomy_db = "silva") {
  # Cargar librerías necesarias
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(sankeyD3)
  
  # Función para calcular abundancia relativa por muestra (porcentaje)
  relabunda <- function(x) {
    as.data.frame(t(t(x) / colSums(x))) * 100
  }
  
  # Identificar columna de taxonomía
  tax_col <- grep("taxonomy|Taxonomy", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("No se encontró una única columna de taxonomía en table")
  
  # Reordenar columnas: dejar taxonomía al final
  table <- table[, c(setdiff(1:ncol(table), tax_col), tax_col)]
  
  # Separar la matriz numérica de la taxonomía
  otu_mat <- table[, -ncol(table)]
  rownames(otu_mat) <- rownames(table)
  
  # Calcular abundancias relativas
  otu_rel <- relabunda(otu_mat)
  
  # Agregar columna de taxonomía nuevamente
  otu_rel$taxonomy <- table[[ncol(table)]]
  
  # Parsear niveles taxonómicos separados por “;” y extraer nombre (quitando prefijos tipo "k__")
  otu_rel_parse <- otu_rel %>%
    tibble::rownames_to_column(var = "Feature.ID") %>%
    tidyr::separate(taxonomy, into = c("k","p","c","o","f","g","s"), sep = ";", fill = "right") %>%
    dplyr::mutate(across(where(is.character), ~ stringr::str_extract(., "[^_]+$")))
  
  # Modificación para Kraken2: concatenar género + especie en la columna especie
  if(tolower(taxonomy_db) == "kraken2") {
    otu_rel_parse <- otu_rel_parse %>%
      dplyr::mutate(
        # Solo si ambas columnas existen y no son NA
        s = ifelse(!is.na(g) & !is.na(s) & s != "NA",
                   paste(g, s, sep = " "),
                   s)
      )
  }
  
  # Función para calcular abundancia media por nivel taxonómico concatenado
  get_level_data <- function(df, level, unite_cols) {
    df %>%
      tidyr::unite(col = !!level, all_of(unite_cols), remove = FALSE) %>%
      dplyr::select(all_of(level), dplyr::where(is.numeric)) %>%
      dplyr::group_by(across(all_of(level))) %>%
      dplyr::summarise(across(where(is.numeric), sum), .groups = "drop") %>%
      tibble::column_to_rownames(var = level) %>%
      t() %>%
      as.data.frame() %>%
      colMeans() %>%
      as.data.frame() %>%
      dplyr::rename(abund = ".") %>%
      dplyr::mutate(taxRank = toupper(substr(level, 1, 1))) %>%
      tibble::rownames_to_column(var = "Taxon") %>%
      dplyr::filter(!endsWith(Taxon, "NA")) %>%
      dplyr::mutate(
        Taxon = gsub("Firmicutes", "Bacillota-D", Taxon),
        Taxon = gsub("Proteobacteria", "Pseudomonadota", Taxon),
        Taxon = gsub("Actinobacteriota", "Actinomycetota", Taxon),
        Taxon = gsub("Cyanobacteria", "Cyanobacteriota", Taxon)
      ) %>%
      tibble::column_to_rownames(var = "Taxon")
  }
  
  # Obtener abundancias medias por cada nivel taxonómico
  bacterias <- get_level_data(otu_rel_parse, "k", "k")
  phylum    <- get_level_data(otu_rel_parse, "phylum", c("k", "p"))
  class     <- get_level_data(otu_rel_parse, "class", c("k", "p", "c"))
  order     <- get_level_data(otu_rel_parse, "order", c("k", "p", "c", "o"))
  family    <- get_level_data(otu_rel_parse, "family", c("k", "p", "c", "o", "f"))
  genus     <- get_level_data(otu_rel_parse, "genus", c("k", "p", "c", "o", "f", "g"))
  specie    <- get_level_data(otu_rel_parse, "specie", c("k", "p", "c", "o", "f", "g", "s"))
  
  # Combinar todos los niveles en un solo data frame
  my_report <- dplyr::bind_rows(bacterias, phylum, class, order, family, genus, specie) %>%
    tibble::rownames_to_column(var = "ids") %>%
    dplyr::mutate(name = stringr::str_extract(ids, "[^_]+$")) %>%
    dplyr::filter(!name == "NA") %>%
    dplyr::filter(!ids %in% c("Bacteria", "Bacteria_Patescibacteria")) %>%
    dplyr::filter(taxRank %in% !!taxRanks) %>%
    dplyr::group_by(taxRank) %>%
    dplyr::slice_max(order_by = abund, n = maxn, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Eliminar posibles nodos vacíos
  my_report <- my_report[!my_report$name %in% c('-_root'), ]
  
  # Separar niveles jerárquicos del ID para armar los enlaces
  splits <- strsplit(my_report$ids, "_")
  sel <- sapply(splits, length) >= 3
  splits <- splits[sel]
  
  # Crear dataframe de enlaces: fuente → destino
  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n = 2))),
                      stringsAsFactors = FALSE)
  colnames(links) <- c("source", "target")
  
  # Asignar valores de abundancia a los enlaces usando el nodo destino
  links$value <- sapply(seq_len(nrow(links)), function(i) {
    abund_val <- my_report$abund[my_report$name == links$target[i]]
    if(length(abund_val) == 1) abund_val else NA_real_
  })
  
  # Eliminar enlaces inválidos
  links <- links[!is.na(links$value) & links$value > 0, ]
  
  # Asignar profundidad por nivel taxonómico (para posicionar nodos horizontalmente)
  taxRank_to_depth <- stats::setNames(seq_along(taxRanks[taxRanks %in% my_report$taxRank]) - 1,
                                      taxRanks[taxRanks %in% my_report$taxRank])
  
  # Crear tabla de nodos
  nodes <- data.frame(
    name = my_report$name,
    depth = taxRank_to_depth[my_report$taxRank],
    value = my_report$abund,
    stringsAsFactors = FALSE
  )
  
  # Mapear nombre → índice (para la visualización con sankeyD3)
  names_id <- stats::setNames(seq_len(nrow(nodes)) - 1, nodes$name)
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  
  # Eliminar enlaces que apuntan a sí mismos o con IDs no encontrados
  links <- links[!is.na(links$source) & !is.na(links$target), ]
  links <- links[links$source != links$target, ]
  
  # Limpiar nombres de nodos (quitar prefijos)
  nodes$name <- sub("^._", "", nodes$name)
  links$type <- sub(" .*", "", nodes[links$source + 1, "name"])
  
  # Crear el diagrama Sankey interactivo
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
    xAxisDomain = taxRanks,
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
  
  # Guardar el archivo HTML
  sankeyD3::saveNetwork(sankey, file = output_file)
  message("Sankey diagram saved to: ", output_file)
}
