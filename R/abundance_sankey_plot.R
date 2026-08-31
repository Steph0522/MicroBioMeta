#' Generate a Sankey diagram from an OTU table with taxonomic information
#'
#' @param table Data frame with one column for taxonomy and the other columns correspond to samples.
#' @param output_file Output HTML file name (default: "sankey.html").
#' @param maxn Maximum number of taxa per level to include in the diagram (default: 25).
#' @param taxRanks Taxonomic levels to display (default: c("D","K","P","C","O","F","G","S")).
#' @param taxonomy_db Database to which the taxonomy in the table corresponds, e.g., "gg" or "kraken2" (default: "gg").
#' @param save_table Logical. If \code{TRUE}, saves a combined table of the
#'   Sankey nodes and links to disk, distinguished by a \code{table_type}
#'   column (\code{"node"} or \code{"link"}). Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"sankey_nodes_links.txt"}.
#'
#' @return Invisibly returns the Sankey diagram object and saves an HTML file.
#' @export
#'
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' abundance_sankey_plot(
#'   table        = table,
#'   output_file  = file.path(tempdir(), "sankey_output.html"),
#'   maxn         = 10,
#'   taxRanks     = c("P", "C", "G", "S"),
#'   taxonomy_db  = "silva"
#' )
#' }


abundance_sankey_plot <- function(table, output_file = "sankey.html", maxn = 25,
                                  taxRanks = c("D","K","P","C","O","F","G","S"),
                                  taxonomy_db = "gg",
                                  save_table = FALSE,
                                  table_filename = "sankey_nodes_links.txt") {

  if (!requireNamespace("networkD3", quietly = TRUE)) {
    stop(
      "Package 'networkD3' is required but not installed.\n",
      "Install it with: install.packages(\"networkD3\")",
      call. = FALSE
    )
  }

  # Internal function for relative abundance
  relabunda <- function(x) as.data.frame(t(t(x) / colSums(x))) * 100

  # Identify the taxonomy column
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  table <- table[, c(setdiff(seq_len(ncol(table)), tax_col), tax_col)]
  
  otu_mat <- table[, -ncol(table)]
  rownames(otu_mat) <- rownames(table)
  otu_rel <- relabunda(otu_mat)
  otu_rel$taxonomy <- table[[ncol(table)]]
  
  otu_rel_parse <- otu_rel %>%
    tibble::rownames_to_column("Feature.ID") %>%
    tidyr::separate(taxonomy, into = c("k","p","c","o","f","g","s"), sep = ";", fill = "right")
  
  # Cleanup according to database
  if(tolower(taxonomy_db) == "silva") {
    otu_rel_parse <- otu_rel_parse %>%
      dplyr::mutate(dplyr::across(c(k,p,c,o,f,g), ~ stringr::str_remove(., "^[a-zA-Z]+__"))) %>%
      dplyr::mutate(s = stringr::str_trim(s),
                    s = stringr::str_replace(s, "^\\s*[a-zA-Z]+__", ""),
                    s = stringr::str_replace_all(s, "_", " "))
  }else if (tolower(taxonomy_db) == "unite") {
    otu_rel_parse <- otu_rel_parse %>%
      # eliminar prefijos tipo k__, p__, c__, etc., pero conservar "incertae sedis"
      dplyr::mutate(dplyr::across(c(k,p,c,o,f,g,s), ~ ifelse(grepl("incertae sedis", .), ., stringr::str_remove(., "^[a-zA-Z]+__")))) %>%
      # reemplazar guiones bajos por espacios en todas las columnas
      dplyr::mutate(dplyr::across(c(k,p,c,o,f,g,s), ~ stringr::str_replace_all(., "_", " ")))
  
  
  } else if (tolower(taxonomy_db) %in% c("gg", "gg2", "greengenes2")) {
    
    otu_rel_parse <- otu_rel_parse %>%
      dplyr::mutate(dplyr::across(c(k,p,c,o,f,g,s), ~ stringr::str_remove(., "^[a-zA-Z]+__"))) %>%
      dplyr::mutate(dplyr::across(c(p,c,o,f,g), ~ stringr::str_replace(., "_[A-Z](_\\d+)?$", ""))) %>%
      dplyr::mutate(dplyr::across(c(k,p,c,o,f,g), ~ stringr::str_replace_all(., "_", " "))) %>%
      
      # limpiar la especie
      dplyr::mutate(
        s = stringr::str_replace(s, "^s__", ""),                   
        s = stringr::str_replace(s, "_[A-Z](_\\d+)?$", ""),       
        s = stringr::str_replace_all(s, "_", " "),               
        s = stringr::str_trim(s)                                   
      )
  
  
  }else if(tolower(taxonomy_db) == "kraken2") {
    otu_rel_parse <- otu_rel_parse %>%
      dplyr::mutate(dplyr::across(where(is.character), ~ stringr::str_extract(., "[^_]+$")),
                    s = ifelse(!is.na(g) & !is.na(s) & s != "NA", paste(g,s,sep=" "), s))
  }
  
  # Internal function to summarize by level and drop empty ones
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
  
  # Construct nodes and links for Sankey
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

  if (save_table) {
    nodes_out <- nodes
    nodes_out$table_type <- "node"
    links_out <- links
    links_out$table_type <- "link"
    combined_table <- dplyr::bind_rows(nodes_out, links_out)
    utils::write.table(combined_table, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  sankey <- networkD3::sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    NodeGroup = "name",
    LinkGroup = "type",
    units = "abund",
    fontFamily = "serif",
    fontSize = 12,
    nodeWidth = 15,
    iterations = 64
  )

  networkD3::saveNetwork(sankey, file = output_file)
  message("Sankey diagram saved to: ", output_file)

  invisible(sankey)
}
