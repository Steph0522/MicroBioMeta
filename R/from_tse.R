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
#'     counts <- matrix(c(10, 0, 5, 3, 8, 1),
#'         nrow = 3,
#'         dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2"))
#'     )
#'     tax <- data.frame(
#'         Kingdom = "Bacteria",
#'         Phylum = c("Firmicutes", "Proteobacteria", NA),
#'         row.names = rownames(counts)
#'     )
#'     samples <- data.frame(
#'         Soil = c("Rhizosphere", "Bulk soil"),
#'         row.names = c("S1", "S2")
#'     )
#'     se <- SummarizedExperiment::SummarizedExperiment(
#'         assays = list(counts = counts), rowData = tax, colData = samples
#'     )
#'     mbm <- from_tse(se)
#'     mbm$table
#'     mbm$metadata
#' }
from_tse <- function(tse, assay_name = "counts") {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Package 'SummarizedExperiment' is required but not installed.\n",
            "Install it with: BiocManager::install(\"SummarizedExperiment\")",
            call. = FALSE
        )
    }
    if (!methods::is(tse, "SummarizedExperiment")) {
        stop("`tse` must be a TreeSummarizedExperiment or SummarizedExperiment object.")
    }
    if (!assay_name %in% SummarizedExperiment::assayNames(tse)) {
        stop("Assay '", assay_name, "' not found. Available assays: ",
            paste(SummarizedExperiment::assayNames(tse), collapse = ", "),
            call. = FALSE
        )
    }

    counts <- as.matrix(SummarizedExperiment::assay(tse, assay_name))
    if (is.null(rownames(counts))) rownames(counts) <- paste0("feature", seq_len(nrow(counts)))
    if (is.null(colnames(counts))) stop("The samples (columns) of `tse` have no names.")

    tax <- as.data.frame(SummarizedExperiment::rowData(tse))
    tax <- if (ncol(tax) == 0) NULL else tax
    samples <- as.data.frame(SummarizedExperiment::colData(tse))
    samples <- if (ncol(samples) == 0) NULL else samples

    .mbm_from_parts(counts, tax, samples)
}
