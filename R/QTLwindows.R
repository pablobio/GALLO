#' Candidate windows identified by GWAS associated with fertility traits in cattle
#'
#' Data from a systematic review which evaluated 18 articles regarding genomw-wide association studies for male fertility traits in beef and dairy cattle
#'
#'@docType data
#'
#'@usage data(QTLwindows)
#'
#'@format A data frame with 50 rows and 8 variables:
#' \itemize{
#'  \item First.marker.in.the.window: First marker mapped in the candidate window
#'  \item Last.marker.in.the.window: Last marker mapped in the candidate window
#'  \item Trait: Trait associated
#'  \item CHR: Chromosome
#'  \item BP1: Chromosomal position in base pairs for the first marker mapped in the candidate window(bovine reference assembly UMD3.1)
#'  \item BP1: Chromosomal position in base pairs for the last marker mapped in the candidate window (bovine reference assembly UMD3.1)
#'  \item Breed: Breed used in the study
#'  \item Reference: Study which the markers were retrieved
#' }
#'
#'@keywords datasets
#'
#'@references Fonseca et al. (2018) Journal of Animal Science, Volume 96, Issue 12, December 2018, Pages 4978–4999.
#'(\href{https://doi.org/10.1093/jas/sky382}{PubMed})
#'@examples
#' data(QTLwindows)
#' @keywords internal
"QTLwindows"
