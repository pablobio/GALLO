#'Candidate markers identified by GWAS associated with fertility traits in cattle
#'
#'bData from a systematic review which evaluated 18 articles regarding genome-wide association studies for male fertility traits in beef and dairy cattle
#'
#'@docType data
#'
#'@usage data(QTLmarkers)
#'
#'@format A data frame with 50 rows and 6 variables:
#'\itemize{
#'  \item Associated.marker: Significantly associated marker
#'  \item Trait: Trait associated
#'  \item CHR: Chromosome
#'  \item BP: Chromosomal position in base pairs (bovine reference assembly UMD3.1)
#'  \item Breed: Breed used in the study
#'  \item Reference: Study which the markers were retrieved
#' }
#'
#'@keywords datasets
#'
#'@references Fonseca et al. (2018) Journal of Animal Science, Volume 96, Issue 12, December 2018, Pages 4978–4999.
#'(\href{https://doi.org/10.1093/jas/sky382}{PubMed})
#'
#'@examples
#'data(QTLmarkers)
#'\donttest{qtl.out <- find_genes_qtls_around_markers(db_file="QTL_db.gff",
#'marker_file=QTLmarkers,method="qtl",
#'marker="snp",interval=100000)}
#'\donttest{head(qtl.out)}
#' @keywords internal
"QTLmarkers"
