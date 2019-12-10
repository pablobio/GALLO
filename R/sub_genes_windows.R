#' Sub-function to search genes around candidate markers
#'
#' Takes a list of candidate markers and search for genes a determined interval
#' @param chr_list "Object with the chromosomes to be analyzed"
#' @param gene Data frame with the information from .gtf file
#' @param markers Data frame with the information from the candidate regions file
#' @param nThreads The number of threads to be used
#' @param int The interval in base pair
#' @name sub_genes_windows
#' @importFrom dplyr do
#' @importFrom data.table setkey
#' @importFrom data.table key
#' @importFrom data.table as.data.table
#' @importFrom parallel detectCores
#' @importFrom dynamicTreeCut printFlush
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @keywords internal
#' @return A dataframe with the genes or QTLs mapped within the specified intervals

sub_genes_windows<-function(chr_list,gene,markers,nThreads=NULL,int=0){
  nCores = detectCores()
  if (is.null(nThreads)) {
    if (nCores < 4)
      nThreads = nCores
    else nThreads = nCores - 1
    pars = list(nThreads)
    names(pars) = .threadAllowVar
    do.call(Sys.setenv, pars)
    registerDoParallel(nThreads)
    invisible(nThreads)
  }
  if (!is.numeric(nThreads) || nThreads < 2){
    stop("nThreads must be numeric and at least 2.")
  }
  if (nThreads > nCores){
    printFlush(paste("Warning in number of threads: Requested number of threads is higher than number\n",
                     "of available processors (or cores).",
                     "It is recommended that the number of threads is no more than number\n",
                     "of available processors.\n"))
    pars = list(nThreads)
    names(pars) = .threadAllowVar
    do.call(Sys.setenv, pars)
    registerDoParallel(nThreads)
    invisible(nThreads)
  }




  # tmp_search.2<-NULL
  foreach::foreach(i=1:length(chr_list),.combine="rbind")%dopar%{ # chr in 1:ncrom

    chr<-chr_list[i]
    tmp_gene<-data.table::as.data.table(gene[which(gene$chr==chr),])
    tmp_markers<-data.table::as.data.table(markers[which(markers$CHR==chr),])
    tmp_markers$tmpBP1<-tmp_markers$BP1-int
    tmp_markers$tmpBP2<-tmp_markers$BP2+int
    cat("\n")
    message(paste("Starting analysis for chromosome ",chr, sep=""))
    cat("\n")

    #selecting the genes wihthin the intervals

    # foverlap requires the second argument to be keyed
    data.table::setkey(tmp_markers, tmpBP1,tmpBP2)

    # find rows where dbh falls between dbh_min and dbh_max, and drop unnecessary
    # columns afterwards
    out<-data.table::foverlaps(tmp_gene, tmp_markers, by.x = c("start_pos","end_pos"), by.y =  data.table::key(tmp_markers),nomatch = 0)
    out[,-c("tmpBP1","tmpBP2")]
  }
}
