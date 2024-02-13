
# Definir la función para convertir la taxonomía
convertir_taxonomia <- function(tax) {
  tax_levels <- strsplit(tax, ",")[[1]]
  new_tax <- character(length(tax_levels))

  for (i in seq_along(tax_levels)) {
    level <- strsplit(tax_levels[i], ":")[[1]]
    prefix <- switch(level[1],
                     "d" = "k__",
                     "p" = "p__",
                     "c" = "c__",
                     "o" = "o__",
                     "f" = "f__",
                     "g" = "g__",
                     "s" = "s__",
                     "sh" = "sh__",
                     "")
    if (prefix == "c__" && level[2] == "Saccharomycetes") {
      new_tax[i] <- "c__Leotiomycetes"
    } else if (prefix == "o__" && level[2] == "Saccharomycetales") {
      new_tax[i] <- "o__Helotiales"
    } else if (prefix == "f__" && level[2] == "Saccharomycetaceae") {
      new_tax[i] <- "f__Helotiaceae"
    } else if (prefix == "g__" && level[2] == "Issatchenkia") {
      new_tax[i] <- "g__Meliniomyces"
    } else {
      new_tax[i] <- paste0(prefix, level[2])
    }
  }

  new_tax[7] <- "s__Helotiales_sp"
  new_tax[8] <- "sh__SH1513725.09FU"

  return(paste(new_tax, collapse = "; "))
}

# Ejemplo de uso
tax <- "d:Fungi,p:Ascomycota,c:Saccharomycetes,o:Saccharomycetales,f:Saccharomycetaceae,g:Issatchenkia,s:Issatchenkia_orientalis_SH1235893.08FU"
tax_convertida <- convertir_taxonomia(tax)
print(tax_convertida)
