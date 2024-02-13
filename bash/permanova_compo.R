#' compositional PerMANOVA
#'
#' @param data_table a data frame containing the data of features counts (in rows) and samples (in columns)
#' @param metadata a data frame containing the metadata
#' @param form a formula specifying the model formula
#' @param strata a logical value indicating whether to use strata in the adonis2 function
#' @param strata_value the name of the column in the metadata to use as strata (if strata=TRUE)
#' @param write a logical value indicating whether to save the results to a file
#'
#' @return results from the perMANOVA
#' @export
#'
#' @examples
#' \dontrun{
#'   permanova_compo(data_table, metadata, form="Material+Tipo", strata_value="Individuo")
#' }
# Define a function to display usage information
print_usage <- function() {
  cat("Usage: Rscript permanova_compo.R [options] data table metadata \"Formula\" strata_value=Column_name write=TRUE\n")
  cat("Options:\n")
  cat("  -h, --help      Show this help message and exit\n")
  cat("  data_table      a txt file containing a table tabulate separated of the data\n")
  cat("                  of features counts (in rows) and samples (in columns)\n")
  cat("  metadata        a txt file with the metada with a column name sampleid\n")
  cat("  \"Formula\"     an expresion in terms of formula to use in the adonis function\n")
  cat("  strata_value    column name in the metada to use strata in adonis function\n")




}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if ("-h" %in% args || "--help" %in% args) {
  print_usage()
  quit(save = "no")
}

permanova_compo <- function(data_table, metadata, form,  strata_value = NULL, write = TRUE) {
  # Load libraries ----------------------------------------------------------
  options(repos = c(CRAN = "https://cloud.r-project.org/"))

  suppressWarnings(suppressMessages({
    if (!requireNamespace("tidyverse", quietly = TRUE))
      install.packages("tidyverse")
    library(tidyverse)
  }))

  suppressWarnings(suppressMessages({
    if (!requireNamespace("compositions", quietly = TRUE))
      install.packages("compositions")
    library(compositions)
  }))

  suppressWarnings(suppressMessages({
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    if (!requireNamespace("ALDEx2", quietly = TRUE))
      BiocManager::install("ALDEx2")
    library(ALDEx2)
  }))


  suppressWarnings(suppressMessages({
    if (!requireNamespace("vegan", quietly = TRUE))
      install.packages("vegan")
    library(vegan)
  }))

  # Format data -------------------------------------------------------------

  rows <- data_table[, 1]
  rownames(data_table) <- rows
  data_table <- data_table[, -1]
  #data_table <- data_table %>% dplyr::select(-taxonomy)

  # Transform data ----------------------------------------------------------

  suppressWarnings(suppressMessages(aldex.clr.transform <-
                                      aldex.clr(
                                        data_table,
                                        mc.samples = 2,
                                        denom = "all",
                                        verbose = FALSE,
                                        useMC = TRUE)))


  aldex.clr.transform.data <-
    t(getMonteCarloSample(aldex.clr.transform, 1))


 suppressMessages(suppressWarnings(tab<- aldex.clr.transform.data %>% as.data.frame() %>% rownames_to_column(
    var = "sampleid") %>% inner_join(metadata)))

 if (is.null(strata_value)) {
   print("No hay variable strata")
 } else {
 print(paste("La variable strata es =",strata_value)) }

 # Perform perMANOVA
 if (is.null(strata_value)) {
   formula_str <- paste("aldex.clr.transform.data", "~", form)
   formula <- as.formula(formula_str)
   perm <- how(nperm = 999)
   perma <- adonis2(
     formula = formula,
     data = tab,
     method = "euclidean",
     permutations = perm
   )
 } else {
   formula_str <- paste("aldex.clr.transform.data", "~", form)
   formula <- as.formula(formula_str)
   perm <- how(nperm = 999)
   perma <- adonis2(
     formula = formula,
     data = tab,
     method = "euclidean",
     permutations = perm,
     strata = tab[[strata_value]]
   )
 }

  # Write data or not ---------------------------------------------------
  if (write) {
    utils::write.table(
      perma,
      paste0("adonis_", format(Sys.time(), "%b_%d_%H_%M_%S"), ".txt"),
      quote = FALSE
    )
    print(perma)
  } else {
    print(perma)
  }

# Return ------------------------------------------------------------------

return(perma)
}
# Lee los argumentos desde la línea de comandos
args <- commandArgs(trailingOnly = TRUE)

# Lee los argumentos de la línea de comandos
data_table <- read.delim(args[1])  # Assuming the first argument is the data table file
metadata <- read.delim(args[2])    # Assuming the second argument is the metadata file
form <- args[3]

# Establece el valor de write según si está presente en los argumentos
write <- "write=TRUE" %in% args

# Establece el valor de strata_value si se proporciona
strata_value <- NULL
for (arg in args) {
  if (grepl("^strata_value=", arg)) {
    strata_value <- gsub("^strata_value=", "", arg)
    break
  }
}

# Realiza el análisis de adonis y guarda los resultados en un archivo si 'write' es TRUE
result <- permanova_compo(data_table, metadata, form = form, strata_value = strata_value, write = write)
