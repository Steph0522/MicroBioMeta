#' compositional PCA
#'
#' @param table a txt file containing a table tabulate separated of the data of features counts (in rows) and samples (in columns)
#' @param write a logical value indicating to save the data imported
#'
#' @return results from the PCA
#' @export
#'
#' @examples
#' \dontrun{
#' pca_compos(table, write=TRUE)
#' }

pca_compos<- function(table, write=TRUE){


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
    if (!requireNamespace("ggpubr", quietly = TRUE))
      install.packages("ggpubr")
    library(ggpubr)
  }))

  # Format data -------------------------------------------------------------

  rows <- table[, 1]
  rownames(table) <- rows
  table <- table[, -1]
  table <- table %>% dplyr::select(-taxonomy)

  # Transform data ----------------------------------------------------------

  suppressWarnings(suppressMessages(aldex.clr.transform <-
                                      aldex.clr(
                                        table,
                                        mc.samples = 2,
                                        denom = "all",
                                        verbose = FALSE,
                                        useMC = TRUE)))


  aldex.clr.transform.data <-
    t(getMonteCarloSample(aldex.clr.transform, 1))


  # PCA ---------------------------------------------------------------------


  pca <- prcomp(aldex.clr.transform.data)

  #LABELS - explanation

  PC1_comp <-
    paste("PC1", round(sum(pca$sdev[1] ^ 2) / mvar(aldex.clr.transform.data) * 100, 1), "%")
  PC2_comp <-
    paste("PC2", round(sum(pca$sdev[2] ^ 2) / mvar(aldex.clr.transform.data) * 100, 1), "%")
  PC3_comp <-
    paste("PC3", round(sum(pca$sdev[3] ^ 2) / mvar(aldex.clr.transform.data) * 100, 1), "%")
  PC4_comp <-
    paste("PC4", round(sum(pca$sdev[4] ^ 2) / mvar(aldex.clr.transform.data) * 100, 1), "%")
  explain<-paste(PC1_comp,",", PC2_comp,",", PC3_comp, ",",PC4_comp)

  # Write data or not -------------------------------------------------------

  if(isTRUE(write)){
    utils::write.table(pca$rotation, paste0("vars_pca_", format(Sys.time(), "%b_%d_%H_%M_%S"), ".txt"), quote = FALSE)
    utils::write.table(pca$x, paste0("individuals_pca_", format(Sys.time(), "%b_%d_%H_%M_%S"), ".txt"), quote = FALSE)
    utils::write.table(explain[1], paste0("explained_pca_", format(Sys.time(), "%b_%d_%H_%M_%S"), ".txt"), quote = FALSE)

  }
  else{
    return(pca)
  }


  # Return ------------------------------------------------------------------

  return(pca)
}


