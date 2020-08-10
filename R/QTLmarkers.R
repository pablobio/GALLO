#'Candidate markers identified by GWAS associated with fertility traits in cattle
#'
#'Data from a systematic review which evaluated 18 articles regarding genome-wide association studies for male fertility traits in beef and dairy cattle
#'
#'@docType data
#'
#'@usage data(QTLmarkers)
#'
#'@format A data frame with 141 rows and 7 variables:
#'\itemize{
#'  \item Associated.marker: Significantly associated marker
#'  \item SNP.reference: The rs ID when available
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
#' data(QTLmarkers)
#' @keywords internal
"QTLmarkers"
