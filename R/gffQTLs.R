#'A gff example for QTL annotation
#'
#'Data from the Animal QTLdb comprasing the bovine QTL annotation
#'
#'@docType data
#'
#'@usage data(gffQTLs)
#'
#'@format A data frame with 111742 rows and 6 variables:
#'\itemize{
#'  \item chr: Chromosome
#'  \item database: The database which the QTL information was retrieved
#'  \item QTL_type: The class of each QTL annotated in the database
#'  \item start_pos: Start position in the genome for each QTL
#'  \item end_pos: End position in the genome for each QTL
#'  \item extra_info: Additional information about the QTLs, such as QTL ID, Name, PUBMED ID, mapping type, among others
#' }
#'
#'@keywords datasets
#'
#'@examples
#'data(gffQTLs)
#' @keywords internal
"gffQTLs"
