#' Genera un diagrama Sankey a partir de una tabla OTU con información taxonómica
#'
#' @param table Data frame con una columna de taxonomía y las demás columnas corresponden a muestras.
#' @param output_file Nombre del archivo HTML de salida (default: "sankey.html").
#' @param maxn Número máximo de taxones por nivel a incluir en el diagrama (default: 25).
#' @param taxRanks Niveles taxonómicos a visualizar (default: c("D","K","P","C","O","F","G","S")).
#' @param taxonomy_db Base de datos taxonómica, ej: "gg" o "kraken2" (default: "gg").

abundance_sankey_plot <- function(table, output_file = "sankey.html", maxn = 25,
                                  taxRanks = c("D","K","P","C","O","F","G","S"),
                                  taxonomy_db = "gg") {
  
  # Función interna para abundancia relativa
  relabunda <- function(x) as.data.frame(t(t(x) / colSums(x))) * 100
  
  # Identificar columna de taxonomía
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  table <- table[, c(setdiff(1:ncol(table), tax_col), tax_col)]
  
  otu_mat <- table[, -ncol(table)]
  rownames(otu_mat) <- rownames(table)
  otu_rel <- relabunda(otu_mat)
  otu_rel$taxonomy <- table[[ncol(table)]]
  
  otu_rel_parse <- otu_rel %>%
    tibble::rownames_to_column("Feature.ID") %>%
    tidyr::separate(taxonomy, into = c("k","p","c","o","f","g","s"), sep = ";", fill = "right")
  
  # Limpieza según base
  if(tolower(taxonomy_db) == "silva") {
    otu_rel_parse <- otu_rel_parse %>%
      dplyr::mutate(dplyr::across(c(k,p,c,o,f,g), ~ stringr::str_remove(., "^[a-zA-Z]+__"))) %>%
      dplyr::mutate(s = stringr::str_trim(s),
                    s = stringr::str_replace(s, "^\\s*[a-zA-Z]+__", ""),
                    s = stringr::str_replace_all(s, "_", " "))
  } else if(tolower(taxonomy_db) == "kraken2") {
    otu_rel_parse <- otu_rel_parse %>%
      dplyr::mutate(dplyr::across(where(is.character), ~ stringr::str_extract(., "[^_]+$")),
                    s = ifelse(!is.na(g) & !is.na(s) & s != "NA", paste(g,s,sep=" "), s))
  }
  
  # Función interna para resumir por nivel y quitar vacíos
  get_level_data <- function(df, level, unite_cols) {
    df %>%
      tidyr::unite(col = !!level, all_of(unite_cols), remove = FALSE) %>%
      dplyr::mutate(!!level := ifelse(grepl("_+$", .data[[level]]) | .data[[level]] == "", NA, .data[[level]])) %>%
      tidyr::drop_na(all_of(level)) %>%
      dplyr::select(c(all_of(level), names(df)[sapply(df,is.numeric)])) %>%
      dplyr::group_by(across(all_of(level))) %>%
      dplyr::summarise(dplyr::across(where(is.numeric), sum), .groups = "drop") %>%
      tibble::column_to_rownames(var = level) %>%
      t() %>%
      as.data.frame() %>%
      colMeans() %>%
      as.data.frame() %>%
      dplyr::rename(abund = ".") %>%
      dplyr::mutate(taxRank = toupper(substr(level,1,1))) %>%
      tibble::rownames_to_column("Taxon") %>%
      dplyr::mutate(
        Taxon = gsub("Firmicutes","Bacillota-D",Taxon),
        Taxon = gsub("Proteobacteria","Pseudomonadota",Taxon),
        Taxon = gsub("Actinobacteriota","Actinomycetota",Taxon),
        Taxon = gsub("Cyanobacteria","Cyanobacteriota",Taxon)
      ) %>%
      dplyr::filter(!endsWith(Taxon,"NA")) %>%
      tibble::column_to_rownames("Taxon")
  }
  
  # Resumir todos los niveles
  bacterias <- get_level_data(otu_rel_parse, "kingdom", "k")
  phylum    <- get_level_data(otu_rel_parse, "phylum", c("k","p"))
  class     <- get_level_data(otu_rel_parse, "class", c("k","p","c"))
  order     <- get_level_data(otu_rel_parse, "order", c("k","p","c","o"))
  family    <- get_level_data(otu_rel_parse, "family", c("k","p","c","o","f"))
  genus     <- get_level_data(otu_rel_parse, "genus", c("k","p","c","o","f","g"))
  specie    <- get_level_data(otu_rel_parse, "specie", c("k","p","c","o","f","g","s"))
  
  my_report <- dplyr::bind_rows(bacterias, phylum, class, order, family, genus, specie) %>%
    tibble::rownames_to_column("ids") %>%
    dplyr::mutate(name = stringr::str_extract(ids, "[^_]+$")) %>%
    dplyr::filter(!name %in% c("NA","Bacteria","Bacteria_Patescibacteria")) %>%
    dplyr::filter(taxRank %in% taxRanks) %>%
    dplyr::group_by(taxRank) %>%
    dplyr::slice_max(order_by = abund, n = maxn, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Construcción de nodos y links para Sankey
  splits <- strsplit(my_report$ids,"_")
  sel <- sapply(splits,length) >= 3
  splits <- splits[sel]
  
  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], 2))),
                      stringsAsFactors = FALSE)
  colnames(links) <- c("source","target")
  
  links$value <- sapply(seq_len(nrow(links)), function(i) {
    val <- my_report$abund[my_report$name == links$target[i]]
    if(length(val) == 1) val else NA_real_
  })
  
  links <- links[!is.na(links$value) & links$value>0, ]
  
  # Profundidad de los niveles
  valid_ranks <- taxRanks[taxRanks %in% unique(my_report$taxRank)]
  taxRank_to_depth <- setNames(seq_along(valid_ranks)-1, valid_ranks)
  
  nodes <- data.frame(
    name = my_report$name,
    depth = taxRank_to_depth[my_report$taxRank],
    value = my_report$abund,
    stringsAsFactors = FALSE
  )
  
  names_id <- setNames(seq_len(nrow(nodes))-1, nodes$name)
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  links <- links[!is.na(links$source) & !is.na(links$target), ]
  links <- links[links$source != links$target, ]
  
  nodes$name <- sub("^._","",nodes$name)
  links$type <- sub(" .*","",nodes[links$source + 1,"name"])
  
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
    xAxisDomain = valid_ranks,
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
