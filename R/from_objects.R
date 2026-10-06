#' Convert a phyloseq object to MicroBioMeta tables
#'
#' Extracts the count table, taxonomy and sample metadata of a
#' \code{phyloseq} object into the plain data frames that every
#' MicroBioMeta function takes.
#'
#' @param physeq A \code{phyloseq} object with an OTU table and, ideally, a
#'   taxonomy table and sample data.
#' @return A list with two data frames:
#'   \describe{
#'     \item{\code{table}}{Taxa in rows, samples in columns, and a last
#'       column \code{taxonomy} with the full taxonomic string
#'       (e.g. \code{"k__Bacteria; p__Firmicutes; c__Bacilli"}).}
#'     \item{\code{metadata}}{Sample metadata whose first column,
#'       \code{SAMPLEID}, holds the sample names of \code{table}.}
#'   }
#' @seealso \code{\link{from_tse}} for \code{(Tree)SummarizedExperiment}
#'   objects.
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   counts <- matrix(c(10, 0, 5, 3, 8, 1), nrow = 3,
#'                    dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2")))
#'   tax <- matrix(c("Bacteria", "Bacteria", "Bacteria",
#'                   "Firmicutes", "Proteobacteria", NA),
#'                 nrow = 3, dimnames = list(rownames(counts), c("Kingdom", "Phylum")))
#'   samples <- data.frame(Soil = c("Rhizosphere", "Bulk soil"),
#'                         row.names = c("S1", "S2"))
#'   ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'                            phyloseq::tax_table(tax),
#'                            phyloseq::sample_data(samples))
#'   mbm <- from_phyloseq(ps)
#'   mbm$table
#'   mbm$metadata
#' }
from_phyloseq <- function(physeq) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.\n",
         "Install it with: BiocManager::install(\"phyloseq\")", call. = FALSE)
  }
  if (!methods::is(physeq, "phyloseq")) stop("`physeq` must be a phyloseq object.")

  counts <- methods::as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) counts <- t(counts)

  tax <- phyloseq::tax_table(physeq, errorIfNULL = FALSE)
  if (!is.null(tax)) tax <- methods::as(tax, "matrix")

  samples <- phyloseq::sample_data(physeq, errorIfNULL = FALSE)
  samples <- if (is.null(samples)) NULL else data.frame(samples, check.names = FALSE)

  .mbm_from_parts(counts, tax, samples)
}

#' Convert a (Tree)SummarizedExperiment object to MicroBioMeta tables
#'
#' Extracts one assay (the counts), the taxonomy (\code{rowData}) and the
#' sample metadata (\code{colData}) of a \code{TreeSummarizedExperiment} or
#' any other \code{SummarizedExperiment}, e.g. one built with the
#' \pkg{mia} package, into the plain data frames that every MicroBioMeta
#' function takes.
#'
#' @param tse A \code{TreeSummarizedExperiment} or \code{SummarizedExperiment}
#'   object with taxa in rows and samples in columns.
#' @param assay_name Character. Name of the assay with the counts. Default
#'   \code{"counts"}.
#' @return A list with two data frames:
#'   \describe{
#'     \item{\code{table}}{Taxa in rows, samples in columns, and a last
#'       column \code{taxonomy} with the full taxonomic string
#'       (e.g. \code{"k__Bacteria; p__Firmicutes; c__Bacilli"}).}
#'     \item{\code{metadata}}{Sample metadata whose first column,
#'       \code{SAMPLEID}, holds the sample names of \code{table}.}
#'   }
#'   The taxonomy is built from the \code{rowData} columns named after
#'   taxonomic ranks (Kingdom or Domain, Phylum, Class, Order, Family,
#'   Genus, Species; any case). If \code{rowData} already has a
#'   \code{taxonomy} column, it is used as is.
#' @seealso \code{\link{from_phyloseq}} for \code{phyloseq} objects.
#' @export
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   counts <- matrix(c(10, 0, 5, 3, 8, 1), nrow = 3,
#'                    dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2")))
#'   tax <- data.frame(Kingdom = "Bacteria",
#'                     Phylum = c("Firmicutes", "Proteobacteria", NA),
#'                     row.names = rownames(counts))
#'   samples <- data.frame(Soil = c("Rhizosphere", "Bulk soil"),
#'                         row.names = c("S1", "S2"))
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = counts), rowData = tax, colData = samples)
#'   mbm <- from_tse(se)
#'   mbm$table
#'   mbm$metadata
#' }
from_tse <- function(tse, assay_name = "counts") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required but not installed.\n",
         "Install it with: BiocManager::install(\"SummarizedExperiment\")",
         call. = FALSE)
  }
  if (!methods::is(tse, "SummarizedExperiment"))
    stop("`tse` must be a TreeSummarizedExperiment or SummarizedExperiment object.")
  if (!assay_name %in% SummarizedExperiment::assayNames(tse))
    stop("Assay '", assay_name, "' not found. Available assays: ",
         paste(SummarizedExperiment::assayNames(tse), collapse = ", "), call. = FALSE)

  counts <- as.matrix(SummarizedExperiment::assay(tse, assay_name))
  if (is.null(rownames(counts))) rownames(counts) <- paste0("feature", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) stop("The samples (columns) of `tse` have no names.")

  tax <- as.data.frame(SummarizedExperiment::rowData(tse))
  tax <- if (ncol(tax) == 0) NULL else tax
  samples <- as.data.frame(SummarizedExperiment::colData(tse))
  samples <- if (ncol(samples) == 0) NULL else samples

  .mbm_from_parts(counts, tax, samples)
}

# Build MicroBioMeta's table (taxonomy string as last column) and metadata
# (SAMPLEID first) from a count matrix, a taxonomy matrix/data frame with
# one column per rank, and a sample data frame with sample names as row names.
.mbm_from_parts <- function(counts, tax, samples) {
  table <- as.data.frame(counts, check.names = FALSE)
  table$taxonomy <- .mbm_taxonomy_string(tax, rownames(counts))

  if (is.null(samples)) {
    metadata <- data.frame(SAMPLEID = colnames(counts))
  } else {
    samples <- samples[colnames(counts), , drop = FALSE]
    samples$SAMPLEID <- NULL
    metadata <- data.frame(SAMPLEID = colnames(counts), samples,
                           check.names = FALSE, row.names = NULL)
  }
  list(table = table, metadata = metadata)
}

# "k__Bacteria; p__Firmicutes; ..." from one column per rank; stops at the
# first missing rank, as QIIME 2 / SILVA strings do
.mbm_taxonomy_string <- function(tax, feature_ids) {
  if (is.null(tax)) {
    warning("No taxonomy found: the `taxonomy` column is left empty.", call. = FALSE)
    return(rep("", length(feature_ids)))
  }
  tax <- as.data.frame(tax, check.names = FALSE)
  tax_col <- grep("^taxonomy$", names(tax), ignore.case = TRUE)
  if (length(tax_col) == 1) return(as.character(tax[[tax_col]]))

  ranks    <- c("kingdom|domain|superkingdom", "phylum", "class", "order",
                "family", "genus", "species")
  prefixes <- c("k", "p", "c", "o", "f", "g", "s")
  cols <- vapply(ranks, function(r) {
    hit <- grep(paste0("^(", r, ")$"), names(tax), ignore.case = TRUE)
    if (length(hit) == 0) NA_integer_ else hit[1]
  }, integer(1))
  if (all(is.na(cols))) {
    stop("No taxonomic rank columns (Kingdom/Domain, Phylum, ..., Species) ",
         "found in the taxonomy.", call. = FALSE)
  }
  prefixes <- prefixes[!is.na(cols)]
  tax <- tax[, cols[!is.na(cols)], drop = FALSE]

  vapply(seq_len(nrow(tax)), function(i) {
    values <- trimws(sub("^([kdpcofgs]__|D_[0-9]+__)", "", as.character(unlist(tax[i, ]))))
    missing <- is.na(values) | values == ""
    n <- if (any(missing)) which(missing)[1] - 1 else length(values)
    if (n == 0) return("Unassigned")
    paste0(prefixes[seq_len(n)], "__", values[seq_len(n)], collapse = "; ")
  }, character(1))
}
